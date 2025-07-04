import subprocess
import os
import sys
import re
import pandas as pd
import numpy as np

# Cell 2: Read the WHAMP output file with proper parsing
def read_whamp_output(filename):
    """
    Read WHAMP output data from text file into a pandas DataFrame.
    
    Parameters:
    filename: path to the text file containing WHAMP output
    
    Returns:
    DataFrame with columns extracted from the WHAMP output format
    """
    try:
        data = []
        
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                
                # Initialize row data
                row = {}
                
                # Extract p value - more specific pattern
                p_match = re.search(r'p=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if p_match:
                    row['p'] = float(p_match.group(1))
                
                # Extract z value - more specific pattern
                z_match = re.search(r'z=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if z_match:
                    row['z'] = float(z_match.group(1))
                
                # Extract f (frequency) - handle scientific notation properly
                # Pattern: f= number number (where numbers can be in scientific notation)
                f_match = re.search(r'f=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if f_match:
                    row['omega_r'] = float(f_match.group(1))
                    row['omega_i'] = float(f_match.group(2))
                
                # Extract EX (electric field X component) - handle spacing properly
                ex_match = re.search(r'EX=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if ex_match:
                    row['EX_real'] = float(ex_match.group(1))
                    row['EX_imag'] = float(ex_match.group(2))
                
                # Extract EY (electric field Y component)
                ey_match = re.search(r'EY=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if ey_match:
                    row['EY_real'] = float(ey_match.group(1))
                    row['EY_imag'] = float(ey_match.group(2))
                
                # Extract EZ (electric field Z component)
                ez_match = re.search(r'EZ=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if ez_match:
                    row['EZ_real'] = float(ez_match.group(1))
                    row['EZ_imag'] = float(ez_match.group(2))
                
                # Extract BETA - handle scientific notation
                beta_match = re.search(r'BETA=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if beta_match:
                    row['BETA'] = float(beta_match.group(1))
                
                # Extract A (alpha parameter) - more specific pattern to avoid matching BETA
                a_match = re.search(r'A=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$', line)
                if a_match:
                    row['A'] = float(a_match.group(1))
                
                # Only add row if we found at least p, z, and frequency data
                if 'p' in row and 'z' in row and 'omega_r' in row:
                    data.append(row)
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        return df
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def run_whamp_automated(model_file, output_file, max_iterations=1900, commands=None):
    """
    Run WHAMP with automated input commands
    
    Parameters:
    - model_file: path to the model file (e.g., '../Models/H17f3c')
    - output_file: path to the output file (e.g., '../results/parallel_firehose5.txt')
    - max_iterations: maximum number of iterations (default: 1900)
    - commands: list of commands to send to WHAMP (default: standard firehose commands)
    """
    
    if commands is None:
        # Default commands for parallel firehose instability analysis
        commands = [
            'p0z.1,1,-.1f1e-4',  # Set P and Z ranges with start frequency
            'pzfewa',            # Set output format
            'S'                  # Stop/quit
        ]
    
    # Build the command
    whamp_cmd = [
        './whamp',
        '-maxiterations', str(max_iterations),
        '-file', model_file,
        '-outfile', output_file
    ]
    
    print(f"Running WHAMP with command: {' '.join(whamp_cmd)}")
    #print(f"Input commands: {commands}")
    
    try:
        # Start the WHAMP process
        process = subprocess.Popen(
            whamp_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,  # Line buffered
            universal_newlines=True
        )
        
        # Send commands to WHAMP
        for i, command in enumerate(commands):
            print(f"Sending command {i+1}: {command}")
            process.stdin.write(command + '\n')
            process.stdin.flush()
            
            # Small delay to allow WHAMP to process the command
            #time.sleep(0.1)
        
        # Wait for the process to complete
        stdout, stderr = process.communicate(timeout=60)  # 60 second timeout
        
        print("WHAMP output:")
        print(stdout)
        
        if stderr:
            print("WHAMP errors:")
            print(stderr)
        
        if process.returncode == 0:
            print(f"WHAMP completed successfully!")
            print(f"Output saved to: {output_file}")
            return True
        else:
            print(f"WHAMP failed with return code: {process.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("WHAMP process timed out!")
        process.kill()
        return False
    except Exception as e:
        print(f"Error running WHAMP: {e}")
        return False

def run_whamp_parameter_sweep(model_file, output_base, a_values, max_iterations=1900):
    """
    Run WHAMP for multiple temperature anisotropy values
    
    Parameters:
    - model_file: path to the model file
    - output_base: base name for output files (will append _A_value.txt)
    - a_values: list of A values to test (e.g., [0.1, 0.2, 0.3, 0.4, 0.5])
    - max_iterations: maximum number of iterations
    """
    
    results = []
    
    for a_val in a_values:
        output_file = f"{output_base}_A_{a_val:.3f}.txt"
        
        commands = [
            'p0z.1,1,-.1f1e-4',  # Set P and Z ranges with start frequency
            'pzfewa',            # Set output format
            f'a{a_val}',         # Set temperature anisotropy
            'p0z.1,1,-.1f1e-4',  # Re-run with new A value
            'S'                  # Stop/quit
        ]
        
        print(f"\n{'='*60}")
        print(f"Running WHAMP for A = {a_val}")
        print(f"{'='*60}")
        
        success = run_whamp_automated(model_file, output_file, max_iterations, commands)
        results.append({
            'A_value': a_val,
            'output_file': output_file,
            'success': success
        })
    
    return results

def run_whamp_interactive_realtime(model_file, output_file, max_iterations=1900):
    """
    Run WHAMP with real-time interactive input (shows WHAMP output as it runs)
    """
    
    whamp_cmd = [
        './whamp',
        '-maxiterations', str(max_iterations),
        '-file', model_file,
        '-outfile', output_file
    ]
    
    print(f"Running WHAMP interactively: {' '.join(whamp_cmd)}")
    
    try:
        process = subprocess.Popen(
            whamp_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Commands to send
        commands = [
            'p0z.1,1,-.1f1e-4',
            'pzfewa',
            'S'
        ]
        
        command_index = 0
        
        # Real-time interaction
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            
            if output:
                print(output.strip())
                
                # Check if WHAMP is asking for input
                if '#INPUT:' in output:
                    if command_index < len(commands):
                        command = commands[command_index]
                        print(f"Sending: {command}")
                        process.stdin.write(command + '\n')
                        process.stdin.flush()
                        command_index += 1
                    else:
                        print("No more commands to send")
                        break
        
        # Wait for completion
        process.wait()
        
        return process.returncode == 0
        
    except Exception as e:
        print(f"Error in interactive mode: {e}")
        return False

# Example usage
if __name__ == "__main__":
    # Make sure we're in the right directory
    if not os.path.exists('./whamp'):
        print("Error: whamp executable not found in current directory")
        print("Please run this script from the whamp build directory")
        sys.exit(1)
    
    # Example 1: Single run
    print("Example 1: Single WHAMP run")
    success = run_whamp_automated(
        model_file='../Models/H17f3c',
        output_file='../results/parallel_firehose5.txt', max_iterations=3000
    )
    
    if success:
        print("\nSingle run completed successfully!")

    # Example 1.5: Parameter sweep for different A values saved in the same file

    a_values = np.logspace(np.log10(0.1), np.log10(1.0), num=20)  # logarithmically spaced between 0.1 and 1.0 inclusive
    c_values = np.logspace(np.log10(2*0.01986), np.log10(0.01986/4), num=40)  # logarithmically spaced cyclotron frequency
    f_values = [1e-1 for _ in range(len(c_values))]  # sweep over frequency
    commands = ['p0z0,1.5,-.1f1e-4', 'pzfewa']
    for c, f in zip(c_values, f_values):
        for a in a_values:
            commands.append(f'c{c}a{a}f{f}')
            if a > 0.8:
                commands.append('z0,.6,-.05')
            else:
                commands.append('z0,2,-.1')
            commands.append('p0z0,1,-.05')
    commands.append('S')
    print("Example 1: Single WHAMP run")
    success = run_whamp_automated(
        model_file='../Models/H17f3c',
        output_file='../results/parallel_firehose_sweep2.txt', max_iterations=3000,
        commands=commands
    )
    
    if success:
        print("\Parameter sweep runs completed successfully!")
    
#     Cell 3: Load your data
filename = '/Users/u0167590/github/whamp/results/parallel_firehose_sweep.txt'
df = read_whamp_output(filename)
# Display basic info about the data
if df is not None:
    print(f"Data loaded successfully!")
    print(f"Shape: {df.shape}")
    print(f"\nColumn names: {list(df.columns)}")
    print(f"\nUnique A values: {sorted(df['A'].unique())}")
    print(f"\nUnique BETA values: {sorted(df['BETA'].unique())}")
    print(f"\nFirst few rows:")
    print(df.head())
    
    # Group by A value to see data structure
    print(f"\nData grouped by A value:")
    for beta_val in sorted(df['BETA'].unique()):
        print(f"\nBETA = {beta_val}:")
        for a_val in sorted(df['A'].unique()):
            subset = df[(df['A'] == a_val) & (df['BETA'] == beta_val)]
            subset1e3 = df[(df['A'] == a_val) & (df['BETA'] == beta_val) & (df['omega_r'] > 1e3)]
            print(f"A = {a_val}, BETA = {beta_val}: {len(subset)} entries with {len(subset1e3)} entries above omega_r > 1e3")
    print(f" Total number of entries {len(df)}, total number of omega_r > 1e3: {len(df[df['omega_r'] > 1e3])}")
else:
    print("Failed to load data")