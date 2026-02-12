import os
import multiprocessing as mp
import re

def extract_seven_digits(filepath):
    """
    Extract the last 7-digit number from the CAF filename.
    Example:
    MiniRun6.5_1E19_RHC.caf.0000000.CAF.root -> 0000000
    """
    basename = os.path.basename(filepath.strip())
    match = re.search(r'\.caf\.(\d{7})\.CAF\.root$', basename)
    if not match:
        raise ValueError(f"Could not extract 7-digit ID from filename: {basename}")
    return match.group(1)

def process_file(filepath):
    filepath = filepath.rstrip('\r\n')

    try:
        z = extract_seven_digits(filepath)
    except ValueError as e:
        print(e)
        return filepath, 1

    cmd = (
        f"UpdateReweight2x2All "
        f"-c GENIEReWeight_100throws.fcl "
        f"-i {filepath} "
        f"-o MiniRun65Nusyst100Universes/MiniRun6.5_1E19_RHC.nuweights.{z}.nusyst.root "
        f"-N 500"
    )

    print(f"Processing file ID {z}")
    print(f"Command: {cmd}")
    result = os.system(cmd)

    return z, result

if __name__ == '__main__':
    # Read the file list
    with open('MiniRun65.txt', 'r') as file:
        lines = file.readlines()

    print(f"Found {len(lines)} files to process")

    num_processes = min(20, len(lines))
    print(f"Using {num_processes} processes")

    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_file, lines)

    print("\nProcessing complete!")
    print("Results:")
    for z, result in results:
        status = "Success" if result == 0 else f"Failed (code: {result})"
        print(f"  File ID {z}: {status}")
