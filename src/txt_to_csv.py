import pandas as pd
import sys

def txt_to_csv(input_file, output_file=None):
    """
    Convert a space-delimited text file to CSV with specified column names.
    
    Parameters:
    input_file: path to input text file
    output_file: path to output CSV file (optional, defaults to input_file.csv)
    """
    
    # Set default output filename if not provided
    if output_file is None:
        output_file = input_file.replace('.txt', '.csv')
        if not output_file.endswith('.csv'):
            output_file += '.csv'
    
    # Read the space-delimited data
    try:
        # Read data, assuming whitespace-separated values
        df = pd.read_csv(input_file, delim_whitespace=True, header=None)
        
        # Assign column names
        df.columns = ['p', 'z', 'omega_r', 'omega_i']
        
        # Save as CSV
        df.to_csv(output_file, index=False)
        
        print(f"Successfully converted {input_file} to {output_file}")
        print(f"Data shape: {df.shape}")
        print("\nFirst few rows:")
        print(df.head())
        
    except Exception as e:
        print(f"Error converting file: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python txt_to_csv.py input_file.txt [output_file.csv]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    txt_to_csv(input_file, output_file)