import os
import subprocess
import sys

args = sys.argv

common_ntherm = int(args[1]) if len(args) > 1 else 0
exp_dir = args[2] if len(args) > 2 else ""
print("updating plots in ", exp_dir if exp_dir != "" else "all directories")

# Set the base_path to the current working directory
base_path = os.getcwd()

scripts = ['plot.py', 'plot_corr.py', 'pcac_mass.py', 'pion_mass.py']

# Loop through all items in the base directory
for item in os.listdir(base_path):
    # Check if the item is a directory and starts with 'T'
    if os.path.isdir(os.path.join(base_path, item)) and item.startswith('T'):
        t_dir_path = os.path.join(base_path, item)
        # Loop through all items in the T directory
        for subitem in os.listdir(t_dir_path):
            # Check if the subitem is a directory and starts with 'k'
            if os.path.isdir(os.path.join(t_dir_path, subitem)) and subitem.startswith('k'):
                # Here, instead of passing the full path, we pass the relative path: T-dir/k-dir
                rel_path = os.path.join(item, subitem)
                if exp_dir == "" or exp_dir == rel_path+"/":
                    for script in scripts:
                        # Construct the command to be executed
                        command = f"python3 {script} {rel_path}/ {common_ntherm}"
                        print("--- executing command: %s ---" % command)
                        # Execute the command
                        subprocess.run(command, shell=True)
