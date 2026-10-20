# Joint conservative measurement and chronology sensitivity; 10,000 simulations.
from pathlib import Path
import subprocess
import sys

file_path = Path(__file__).resolve().parent / "modules/Code06_joint_uncertainty/run_analysis.py"
subprocess.run([sys.executable, str(file_path), *sys.argv[1:]], check=True)
