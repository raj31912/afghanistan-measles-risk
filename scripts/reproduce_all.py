#!/usr/bin/env python3
"""One-command reproducibility check: analysis -> validation -> tests -> figures."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def run(*args):
    print("+", " ".join(map(str,args)), flush=True)
    subprocess.run([str(x) for x in args],cwd=ROOT,check=True)

def main():
    run(sys.executable,"scripts/run_analysis.py")
    run(sys.executable,"scripts/validate_outputs.py")
    run(sys.executable,"scripts/audit_locked_source.py")
    run(sys.executable,"-m","pytest","-q")
    run(sys.executable,"scripts/generate_figures.py")
    print("\nREPRODUCIBILITY CHECK PASSED")

if __name__=="__main__": main()
