import os
import subprocess
import sys

L_values = [15, 16, 18, 20, 24, 30, 32, 36] # 40

for L in L_values:
    T = 2*L
    command = f"python tune_pcac_plot.py T{T}_L{L}_A/"
    print("--- executing command: %s ---" % command)
    subprocess.run(command, shell=True)

