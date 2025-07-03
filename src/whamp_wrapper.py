import subprocess
import time
import os
import sys
from pathlib import Path



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
    print(f"Input commands: {commands}")
    
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
        output_file='../results/parallel_firehose5.txt'
    )
    
    if success:
        print("\nSingle run completed successfully!")

    # Example 1.5: Parameter sweep for different A values saved in the same file
    a_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # sweep over T_perp/T_par
    c_values = [0.01986/i for i in range(1,5)]  # sweep over cyclotron frequency
    f_values = [1e-4 for i in range(1,5)]  # sweep over frequency
    commands = ['p0z.1,1,-.1f1e-4', 'pzfewa']
    for c, f in zip(c_values, f_values):
        for a in a_values:
            commands.append(f'c{c}a{a}f{f}')
            commands.append('p0z.1,1,-.1')
    commands.append('S')
    print("Example 1: Single WHAMP run")
    success = run_whamp_automated(
        model_file='../Models/H17f3c',
        output_file='../results/parallel_firehose_sweep.txt',
        commands=commands
    )
    
    if success:
        print("\Parameter sweep runs completed successfully!")
    
    # Example 2: Parameter sweep for different A values
    #print("\n" + "="*60)
    #print("Example 2: Parameter sweep for different A values")
    
    #a_values = [0.1, 0.2, 0.3, 0.4, 0.5]
    #results = run_whamp_parameter_sweep(
    #    model_file='../Models/H17f3c',
    #    output_base='../results/parallel_firehose_sweep',
    #    a_values=a_values
    #)
    
    #print("\nParameter sweep results:")
    #for result in results:
    #    status = "✓" if result['success'] else "✗"
    #    print(f"{status} A = {result['A_value']:.3f} -> {result['output_file']}")
    
    # Example 3: Interactive mode with real-time output
    #print("\n" + "="*60)
    #print("Example 3: Interactive mode")
    
    #success = run_whamp_interactive_realtime(
    #    model_file='../Models/H17f3c',
    #    output_file='../results/parallel_firehose_interactive.txt'
    #)
    
    #if success:
    #    print("\nInteractive run completed successfully!")