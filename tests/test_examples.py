import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(__file__))


def run(cmd, env=None, timeout=60):
    res = subprocess.run(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    return res.returncode, res.stdout, res.stderr


def test_use_case_basic_runs():
    env = os.environ.copy()
    env["PYTHONPATH"] = ROOT
    code, out, err = run([sys.executable, "use_cases/use_case_basic.py"], env=env)
    assert code == 0, f"basic failed: {err}\n{out}"
    assert "Entanglement:" in out


def test_use_case_visualize_runs(tmp_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = ROOT
    # use small grid internally; script already uses a small grid and writes a fixed filename
    code, out, err = run([sys.executable, "use_cases/use_case_visualize.py"], env=env)
    assert code == 0, f"visualize failed: {err}\n{out}"
    assert "entanglement_proj.png" in out or os.path.exists(os.path.join(ROOT, "entanglement_proj.png"))


def test_use_case_benchmark_runs():
    env = os.environ.copy()
    env["PYTHONPATH"] = ROOT
    code, out, err = run([sys.executable, "use_cases/use_case_benchmark.py"], env=env)
    assert code == 0, f"benchmark failed: {err}\n{out}"
    assert "time:" in out.lower()
