import os
import subprocess
import sys

args = sys.argv

ntherm = int(args[1]) if len(args) > 1 else 0

temp_values = [0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0] 

for temp in temp_values:
    print(f"\n---------- TEMP = {temp} ----------")
    command = f"python continuum_limit.py {temp} {ntherm}"
    print(f"executing command: {command}")
    subprocess.run(command, shell=True)

