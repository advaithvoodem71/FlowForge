#!/usr/bin/env python3
"""
run_openfoam_batch.py

Batch-run OpenFOAM meshing + solver for a directory of STL files.

Requirements / assumptions:
- You already have a prepared OpenFOAM case template directory (see README below).
- OpenFOAM CLI tools (blockMesh, surfaceFeatureExtract, snappyHexMesh, checkMesh, simpleFoam)
  are available in PATH when this script runs (or full paths provided in CFG).
- Template case must include correctly configured:
    - system/snappyHexMeshDict referencing the STL by name (or using wildcard)
    - constant/triSurface/ directory (can be empty; this script will copy .stl files into it)
    - functionObjects in system/controlDict for forceCoeffs (or the postProcessing folder will exist)
- This script runs cases serially. For parallel runs, modify commands to use mpirun / decomposePar.

Usage:
    python run_openfoam_batch.py --stl_dir ./output_stls --template ./foam_template --work_dir ./cases

Outputs:
- cases/<case_name>/ : OpenFOAM case folder with logs
- cases/<case_name>/results_summary.csv : extracted forceCoeffs/forces summary (lift/drag)
- cases/logs/: aggregated logs
"""

import argparse
import os
import shutil
import subprocess
import glob
import csv
import time
from pathlib import Path


# Configuration

CFG = {
    "blockMesh_cmd": "blockMesh",
    "surfaceFeatureExtract_cmd": "surfaceFeatureExtract",
    "snappyHexMesh_cmd": "snappyHexMesh",
    "checkMesh_cmd": "checkMesh",
    "solver_cmd": "simpleFoam",  # change to rhoSimpleFoam or steady solver if needed
    "log_dir": "logs",
    "case_parent_dir": "cases",   # where case copies will be created
    "timeout": 60 * 60 * 2,       # 2 hours default timeout per case (seconds)
}


# Helpers


def run_cmd(cmd, cwd=None, logfile=None, timeout=None):
    """
    Run a shell command list `cmd` in working directory `cwd`.
    If logfile is provided (open file object), stdout+stderr are redirected there.
    Raises CalledProcessError on non-zero return.
    """
    start = time.time()
    with (open(logfile, "a") if logfile else open(os.devnull, "w")) as lf:
        print(f"[RUN] {' '.join(cmd)}  (cwd={cwd})")
        proc = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout)
    end = time.time()
    print(f"[DONE] {' '.join(cmd)} (elapsed {end-start:.1f}s)")
    # subprocess.run with check=True would raise; here we rely on returncode
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)

def safe_mkdir(path):
    os.makedirs(path, exist_ok=True)


# Case preparation & run


def prepare_case_from_template(stl_path, template_dir, work_parent, case_name):
    """
    Copy template_dir -> work_parent/case_name and copy stl_path into constant/triSurface/.
    Returns path to created case directory.
    """
    case_dir = os.path.join(work_parent, case_name)
    if os.path.exists(case_dir):
        print(f"[WARN] Case dir exists; removing: {case_dir}")
        shutil.rmtree(case_dir)
    print(f"[INFO] Copying template '{template_dir}' -> '{case_dir}'")
    shutil.copytree(template_dir, case_dir)
    # ensure triSurface directory exists
    tri_dir = os.path.join(case_dir, "constant", "triSurface")
    safe_mkdir(tri_dir)
    # copy stl
    stl_basename = os.path.basename(stl_path)
    dst_stl = os.path.join(tri_dir, stl_basename)
    print(f"[INFO] Copying STL '{stl_path}' -> '{dst_stl}'")
    shutil.copy2(stl_path, dst_stl)
    return case_dir, stl_basename

def run_case_pipeline(case_dir, stl_name, logfile_path, timeout=CFG["timeout"]):
    """
    Run blockMesh, surfaceFeatureExtract, snappyHexMesh, checkMesh, solver.
    Append outputs to logfile_path.
    """
    safe_mkdir(os.path.dirname(logfile_path))

    # 1) blockMesh
    run_cmd([CFG["blockMesh_cmd"]], cwd=case_dir, logfile=logfile_path, timeout=timeout)

    # 2) surfaceFeatureExtract (makes featureEdges for snappy)
    run_cmd([CFG["surfaceFeatureExtract_cmd"]], cwd=case_dir, logfile=logfile_path, timeout=timeout)

    # 3) snappyHexMesh (overwrite)
    run_cmd([CFG["snappyHexMesh_cmd"], "-overwrite"], cwd=case_dir, logfile=logfile_path, timeout=timeout)

    # 4) checkMesh
    run_cmd([CFG["checkMesh_cmd"]], cwd=case_dir, logfile=logfile_path, timeout=timeout)

    # 5) solver (steady RANS solver)
    run_cmd([CFG["solver_cmd"]], cwd=case_dir, logfile=logfile_path, timeout=timeout)

def find_force_coeffs_output(case_dir):
    """
    Attempt to locate functionObject outputs for forceCoeffs or forces.
    Common locations:
      - postProcessing/forceCoeffs/0/forceCoeffs.dat
      - postProcessing/forces/0/forces.dat
      - forces.dat (in case log)
    Returns path or None.
    """
    candidates = []
    candidates.extend(glob.glob(os.path.join(case_dir, "postProcessing", "forceCoeffs", "*", "forceCoeffs.dat")))
    candidates.extend(glob.glob(os.path.join(case_dir, "postProcessing", "forces", "*", "forces.dat")))
    candidates.extend(glob.glob(os.path.join(case_dir, "postProcessing", "*", "*forces*.dat")))
    # also check top-level
    candidates.extend(glob.glob(os.path.join(case_dir, "*.dat")))
    if candidates:
        # pick most recent
        candidates = sorted(candidates, key=os.path.getmtime, reverse=True)
        return candidates[0]
    return None

def extract_coeffs_from_file(fpath):
    """
    Extract last line of a forceCoeffs.dat style file, return dict of parsed numbers.
    Many forceCoeffs files have header then rows: time CL CD CM etc or Fx Fy Fz.
    We'll try to parse the last non-comment line.
    """
    if not fpath or not os.path.exists(fpath):
        return None
    with open(fpath, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return None
    last = lines[-1]
    parts = last.split()
    # heuristics: if first token is time (float) followed by 3+ floats -> attempt mapping
    nums = []
    for p in parts:
        try:
            nums.append(float(p))
        except:
            pass
    if len(nums) < 2:
        return None
    # Return everything as list + original line
    return {"raw_line": last, "values": nums, "file": fpath}


# Batch driver


def batch_run(stl_dir, template_dir, work_parent, out_summary_csv):
    stl_dir = os.path.abspath(stl_dir)
    template_dir = os.path.abspath(template_dir)
    work_parent = os.path.abspath(work_parent)
    safe_mkdir(work_parent)
    safe_mkdir(CFG["log_dir"])

    stl_paths = sorted([p for p in glob.glob(os.path.join(stl_dir, "*.stl"))])
    if not stl_paths:
        raise RuntimeError(f"No STL files found in {stl_dir}")

    summary_rows = []
    for stl in stl_paths:
        name = Path(stl).stem
        case_name = f"case_{name}"
        print("\n" + "="*60)
        print(f"[BATCH] Starting case for STL: {stl} -> {case_name}")

        case_dir, stl_basename = prepare_case_from_template(stl, template_dir, work_parent, case_name)
        logfile = os.path.join(CFG["log_dir"], f"{case_name}.log")
        try:
            run_case_pipeline(case_dir, stl_basename, logfile, timeout=CFG["timeout"])
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Command failed for case {case_name}: {e}. See log {logfile}")
            summary_rows.append({
                "case": case_name,
                "stl": stl,
                "status": "failed_cmd",
                "error": str(e),
                "force_file": "",
                "coeffs": ""
            })
            continue
        except subprocess.TimeoutExpired as e:
            print(f"[ERROR] Timeout for case {case_name}: {e}. See log {logfile}")
            summary_rows.append({
                "case": case_name,
                "stl": stl,
                "status": "timeout",
                "error": str(e),
                "force_file": "",
                "coeffs": ""
            })
            continue

        # locate forceCoeffs / forces output
        force_file = find_force_coeffs_output(case_dir)
        coeffs = extract_coeffs_from_file(force_file) if force_file else None

        row = {
            "case": case_name,
            "stl": stl,
            "status": "ok" if coeffs else "ok_no_coeffs",
            "force_file": force_file if force_file else "",
            "coeffs": coeffs["values"] if coeffs else ""
        }
        summary_rows.append(row)
        # write per-case CSV summary
        csvcase = os.path.join(case_dir, "results_summary.csv")
        with open(csvcase, "w", newline="") as cf:
            writer = csv.writer(cf)
            writer.writerow(["case", "stl", "status", "force_file", "coeffs_raw"])
            writer.writerow([row["case"], row["stl"], row["status"], row["force_file"], coeffs["raw_line"] if coeffs else ""])
        print(f"[BATCH] Completed case {case_name}. Summary -> {csvcase}")

    # write master summary
    with open(out_summary_csv, "w", newline="") as outf:
        writer = csv.writer(outf)
        writer.writerow(["case", "stl", "status", "force_file", "coeffs"])
        for r in summary_rows:
            writer.writerow([r["case"], r["stl"], r["status"], r["force_file"], " ".join(map(str, r["coeffs"])) if r["coeffs"] else ""])
    print(f"[DONE] Batch complete. Master summary: {out_summary_csv}")


# CLI


def main():
    parser = argparse.ArgumentParser(description="Batch-run OpenFOAM cases for STL files")
    parser.add_argument("--stl_dir", required=True, help="Directory containing STL files")
    parser.add_argument("--template", required=True, help="OpenFOAM case template directory")
    parser.add_argument("--work_dir", default=CFG["case_parent_dir"], help="Parent folder to write case copies")
    parser.add_argument("--out_summary", default="batch_summary.csv", help="CSV summary path")
    args = parser.parse_args()

    batch_run(args.stl_dir, args.template, args.work_dir, args.out_summary)

if __name__ == "__main__":
    main()