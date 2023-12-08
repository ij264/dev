import sys
import subprocess

# Retrieve command-line arguments
if len(sys.argv) != 3:
    print("Usage: python main.py <viscosity_kernel> <kernel_type>")
    sys.exit(1)

valid_kernels = ['surftopo', 'geoid']

visc_kernel_src = sys.argv[1]
kernel_type = sys.argv[2]


assert isinstance(visc_kernel_src, str), "Kernel source must be a string."
assert kernel_type in valid_kernels, f"Kernel type must be one of {valid_kernels}."

# Your script logic using the arguments
print(f"Argument 1: {visc_kernel_src}")

# Replace 'your_command_here' with the command you want to run
command = './MAKE_KERNEL OUTPUT_const_visc/ VISC_INPUTS/const_visc.vis'

# Run the command in the shell
try:
    result_1 = subprocess.run('ls', shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(result_1.stdout)
    result_2 = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(result_2.stdout)
    result_3 = subprocess.run('cd ..', shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(result_3.stdout)
    print("Kernels produced successfully")
except subprocess.CalledProcessError as e:
    print(f"Command failed with error: {e}")
    print(e.stderr)