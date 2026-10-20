# Fixed-record and joint measurement-depth rarefaction; 10,000 iterations.
from pathlib import Path
import subprocess
import sys

file_path = Path(__file__).resolve().parent / "modules/Code03_sampling_density_rarefaction/run.py"
subprocess.run([sys.executable, str(file_path), *sys.argv[1:]], check=True)
