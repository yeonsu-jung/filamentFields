#!/usr/bin/env python3
"""
Script to generate polyscope screenshots from rod coordinate CSV files for movie creation.
Based on example.py structure.
"""

import numpy as np
import polyscope as ps
import os
import glob
from pathlib import Path
import re
import argparse
from PIL import Image, ImageDraw, ImageFont
from typing import List


def prep_for_polyscope(r_list, num_rods):
    """
    Prepare data for polyscope visualization.
    Based on the function from example.py and visualizations.py
    """
    nodes = np.vstack(r_list)
    colors = np.array([
        [76, 153, 204],   # light blue
        [204, 76, 153],   # pinkish red
        [76, 204, 153],   # mint green
        [153, 204, 76],   # light olive green
        [204, 153, 76],   # goldenrod
        [153, 76, 204],   # medium purple
        [204, 76, 102],   # crimson
        [76, 204, 204],   # cyan
        [204, 204, 76],   # sunflower yellow
        [102, 76, 204]    # indigo
    ]) / 255
    
    colors_interpolated = np.zeros((len(r_list), 3))
    for i in range(3):
        colors_interpolated[:, i] = np.interp(
            np.linspace(0, 1, num_rods),
            np.linspace(0, 1, colors.shape[0]),
            colors[:, i]
        )

    edge_colors = []
    edges = []

    starting_index = 0
    for i in range(len(r_list)):
        num_edges = len(r_list[i]) - 1
        to_add = [[starting_index + j, starting_index + j + 1] for j in range(num_edges)]
        edges.append(to_add)
        edge_colors.append(np.array([colors_interpolated[i] for j in range(num_edges)]))
        starting_index += num_edges + 1
    edges = np.vstack(edges)
    edge_colors = np.vstack(edge_colors)
    
    return nodes, edges, edge_colors


def load_rod_data(csv_path, num_nodes=9, num_rods=1536):
    """
    Load rod coordinate data from CSV file.
    """
    data = np.loadtxt(csv_path, delimiter=',')
    all_rods = data.reshape(num_rods, num_nodes, 3)
    # Flip z coordinate as done in example.py
    # all_rods[:, :, 2] *= -1
    return all_rods


def make_synthetic_rods(num_rods: int = 100, num_nodes: int = 30, seed: int = 0) -> np.ndarray:
    """Generate a set of random-walk rods (num_rods, num_nodes, 3)."""
    rng = np.random.default_rng(seed)
    rods: List[np.ndarray] = []
    for _ in range(num_rods):
        x = np.cumsum(rng.standard_normal(num_nodes))
        y = np.cumsum(rng.standard_normal(num_nodes))
        z = np.cumsum(rng.standard_normal(num_nodes))
        rods.append(np.vstack([x, y, z]).T)
    return np.asarray(rods)


def get_sorted_csv_files(positions_dir):
    """
    Get all CSV files sorted by their numerical timestamp.
    """
    csv_pattern = os.path.join(positions_dir, "rod_coordinates_*.csv")
    csv_files = glob.glob(csv_pattern)
    
    # Extract numerical part for sorting
    def extract_number(filename):
        match = re.search(r'rod_coordinates_(\d+)\.csv', filename)
        return int(match.group(1)) if match else 0
    
    # Sort by numerical timestamp
    csv_files.sort(key=extract_number)
    return csv_files


def generate_movie_screenshots(positions_dir, output_dir, skip_frames=1, rod_radius=1.21):
    """
    Generate polyscope screenshots from all CSV files for movie creation.
    
    Parameters:
    - positions_dir: Directory containing CSV files
    - output_dir: Directory to save screenshots
    - skip_frames: Only process every nth frame (for faster processing)
    - rod_radius: Radius of the rods for visualization
    """
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Get sorted CSV files
    csv_files = get_sorted_csv_files(positions_dir)
    print(f"Found {len(csv_files)} CSV files")
    if not csv_files:
        print("No CSV files found in positions directory.")
        return 0
    
    # Initialize polyscope
    ps.init()
    # ps.set_autoscale_structures(False)
    # ps.set_automatically_compute_scene_extents(False)
    ps.set_ground_plane_mode("none")
    # ps.set_length_scale(2.)
    # ps.set_bounding_box((-128., -128., -100.), (128., 128., 100.))
    ps.set_up_dir("z_up")
    
    # Load first frame to set up the visualization
    print("Loading first frame for setup...")
    first_rods = load_rod_data(csv_files[0])
    num_rods = first_rods.shape[0]
    print(f"Number of rods: {num_rods}")
    
    # Set up polyscope with first frame
    nodes, edges, edge_colors = prep_for_polyscope(first_rods, num_rods)
    ps_curves = ps.register_curve_network("filaments", nodes, edges)
    ps_curves.set_radius(rod_radius / 2, relative=False)
    ps_curves.add_color_quantity("edge_colors", edge_colors, defined_on='edges', enabled=True)
    
    # Process each frame
    frame_count = 0
    for i, csv_file in enumerate(csv_files[::skip_frames]):
    # for i, csv_file in enumerate(csv_files):
        print(f"Processing frame {frame_count} / {len(csv_files[::skip_frames])}: {os.path.basename(csv_file)}")
        
        # Load rod data
        all_rods = load_rod_data(csv_file)
        
        # Update polyscope visualization
        nodes, edges, edge_colors = prep_for_polyscope(all_rods, num_rods)
        ps_curves.update_node_positions(nodes)
        ps_curves.add_color_quantity("edge_colors", edge_colors, defined_on='edges', enabled=True)
        
        # Take screenshot
        screenshot_path = os.path.join(output_dir, f"frame_{frame_count:06d}.png")
        ps.screenshot(screenshot_path)
        
        frame_count += 1
    
    print(f"Generated {frame_count} screenshots in {output_dir}")
    print("To create a movie, run:")
    print(f"ffmpeg -framerate 30 -i {output_dir}/frame_%06d.png -c:v libx264 -pix_fmt yuv420p -y {output_dir}/output.mp4")
    # ffmpeg -framerate 30 -i /Users/yeonsu/GitHub/filamentFields/movie_frames/frame_000000.png -c:v libx264 -pix_fmt yuv420p -y /Users/yeonsu/GitHub/filamentFields/movie_frames/output.mp4
    return frame_count


def main():
    parser = argparse.ArgumentParser(description="Generate polyscope screenshots from rod coordinate CSV files.")
    parser.add_argument("--positions-dir", type=str, default="/Users/yeonsu/Downloads/nest_packing/positions",
                        help="Directory containing rod_coordinates_*.csv files")
    parser.add_argument("--output-dir", type=str, default="./movie_frames",
                        help="Directory to write screenshots")
    parser.add_argument("--skip-frames", type=int, default=3000,
                        help="Process every Nth frame (1 = process all)")
    parser.add_argument("--rod-radius", type=float, default=1.21,
                        help="Rod radius used for rendering")
    parser.add_argument("--demo", action="store_true",
                        help="If set, generate synthetic rods instead of reading CSVs")
    parser.add_argument("--entanglement", action="store_true",
                        help="Compute total entanglement for each frame and overlay on the image")
    parser.add_argument("--ent-R-omega", type=float, default=2000.0,
                        help="Neighborhood radius R_omega used when computing total entanglement")
    args = parser.parse_args()

    positions_dir = args.positions_dir
    output_dir = args.output_dir
    skip_frames = args.skip_frames
    rod_radius = args.rod_radius

    print("Starting movie generation...")
    print(f"Output directory: {output_dir}")
    print(f"Skip frames: {skip_frames}")

    def annotate(path: str, text: str):
        try:
            img = Image.open(path).convert("RGBA")
            draw = ImageDraw.Draw(img)
            # pick a basic font; default if truetype not found
            try:
                font = ImageFont.truetype("Arial.ttf", 28)
            except Exception:
                font = ImageFont.load_default()
            # draw outline for contrast (white halo), then black text on top
            x, y = 10, 10
            draw.text((x+1, y+1), text, fill=(255,255,255,200), font=font)
            draw.text((x, y), text, fill=(0,0,0,255), font=font)
            img.save(path)
        except Exception as e:
            print(f"Warning: failed to annotate image {path}: {e}")

    if args.demo:
        print("Demo mode: generating synthetic rods")
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        ps.init()
        ps.set_ground_plane_mode("none")
        ps.set_up_dir("z_up")
        rods = make_synthetic_rods(num_rods=100, num_nodes=30)
        num_rods = rods.shape[0]
        nodes, edges, edge_colors = prep_for_polyscope(list(rods), num_rods)
        ps_curves = ps.register_curve_network("filaments", nodes, edges)
        ps_curves.set_radius(rod_radius / 2, relative=False)
        ps_curves.add_color_quantity("edge_colors", edge_colors, defined_on='edges', enabled=True)
        screenshot_path = os.path.join(output_dir, f"frame_{0:06d}.png")
        ps.screenshot(screenshot_path)
        if args.entanglement:
            try:
                import filamentFields  # imported only if needed
                fil_list = [rods[i, :, :] for i in range(rods.shape[0])]
                fF = filamentFields.filamentFields(fil_list)
                fF.precompute(args.ent_R_omega)
                E = fF.return_total_entanglement()
                annotate(screenshot_path, f"Entanglement: {E:.3g}")
            except Exception as e:
                print(f"Warning: entanglement computation failed: {e}")
        print(f"Generated 1 screenshot in {output_dir}")
        return 0

    if not os.path.isdir(positions_dir):
        print(f"Error: positions directory not found: {positions_dir}")
        print("Tip: pass --positions-dir /path/to/positions or run with --demo")
        return 1

    print(f"Input directory: {positions_dir}")
    # CSV mode: loop ourselves to compute entanglement per-frame if requested
    csv_files = get_sorted_csv_files(positions_dir)[::skip_frames]
    if not csv_files:
        print("No CSV files found in positions directory.")
        print("Tip: run with --demo to generate a synthetic example.")
        return 1

    ps.init()
    ps.set_ground_plane_mode("none")
    ps.set_up_dir("z_up")
    # initialize with first frame
    first_rods = load_rod_data(csv_files[0])
    num_rods = first_rods.shape[0]
    nodes, edges, edge_colors = prep_for_polyscope(first_rods, num_rods)
    ps_curves = ps.register_curve_network("filaments", nodes, edges)
    ps_curves.set_radius(rod_radius / 2, relative=False)
    ps_curves.add_color_quantity("edge_colors", edge_colors, defined_on='edges', enabled=True)
    entanglement_over_time = []
    for frame_idx, csv_file in enumerate(csv_files):
        all_rods = load_rod_data(csv_file)
        nodes, edges, edge_colors = prep_for_polyscope(all_rods, num_rods)
        ps_curves.update_node_positions(nodes)
        ps_curves.add_color_quantity("edge_colors", edge_colors, defined_on='edges', enabled=True)

        screenshot_path = os.path.join(output_dir, f"frame_{frame_idx:06d}.png")
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        ps.screenshot(screenshot_path)

        if args.entanglement:
            try:
                real_idx = frame_idx * skip_frames
                import filamentFields
                fil_list = [all_rods[i, :, :] for i in range(all_rods.shape[0])]
                fF = filamentFields.filamentFields(fil_list)
                fF.precompute(args.ent_R_omega)
                E = fF.return_total_entanglement()
                annotate(screenshot_path, f"Frame: {real_idx}, Entanglement: {E:.3g}")
                # real index
                entanglement_over_time.append((real_idx, E))
            except Exception as e:
                print(f"Warning: entanglement computation failed for {os.path.basename(csv_file)}: {e}")

    # save entanglement over time if computed
    if args.entanglement and entanglement_over_time:
        ent_path = os.path.join(output_dir, "entanglement_over_time.csv")
        with open(ent_path, 'w') as f:
            f.write("frame_index,entanglement\n")
            for idx, E in entanglement_over_time:
                f.write(f"{idx},{E}\n")
        print(f"Saved entanglement over time to {ent_path}")

    frames = len(csv_files)
    if frames == 0:
        print("Tip: run with --demo to generate a synthetic example.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
