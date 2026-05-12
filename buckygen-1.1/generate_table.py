import subprocess
import re
import sys
import os
import csv

# Configuration
START_ATOMS = 60
END_ATOMS = 400
STEP = 2
CSV_FILE = "ipr_results.csv"
BUCKYGEN_EXE_BASE = "buckygen"
BUCKYGEN_SRC = "buckygen.c"

def get_exe_path():
    """Returns the correct path to the buckygen executable."""
    ext = ".exe" if os.name == 'nt' else ""
    path = os.path.join(".", BUCKYGEN_EXE_BASE + ext)
    if os.path.exists(path):
        return path
    # Try without ./ just in case
    path = BUCKYGEN_EXE_BASE + ext
    if os.path.exists(path):
        return path
    return None

def compile_buckygen():
    """Compiles buckygen if the executable doesn't exist."""
    if get_exe_path() is None:
        print(f"Compiling {BUCKYGEN_SRC}...")
        # -DCPUTIME=0 bypasses the sys/times.h dependency
        cmd = ["gcc", "-O3", "-DCPUTIME=0", BUCKYGEN_SRC, "-o", BUCKYGEN_EXE_BASE]
        try:
            subprocess.run(cmd, check=True)
        except Exception as e:
            print(f"Error compiling: {e}. Ensure gcc is in your PATH.")
            sys.exit(1)

def get_buckygen_count(atoms, switches):
    """Runs buckygen for a specific atom count and returns the integer count."""
    exe = get_exe_path()
    if not exe:
        return 0

    # Removed -q to ensure we get the count from stderr/stdout
    cmd = [exe, f"{atoms}d", "-u"] + switches
    try:
        # Buckygen often writes the count to stderr
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        combined_output = result.stdout + "\n" + result.stderr
        
        # Match "Number of fullerenes generated with X vertices: Y"
        match = re.search(r"generated with \d+ vertices: (\d+)", combined_output)
        if match:
            return int(match.group(1))
        
        # Fallback: Match "Y fullerenes generated"
        match = re.search(r"(\d+) fullerenes generated", combined_output)
        if match:
            return int(match.group(1))
            
    except Exception as e:
        # print(f"Debug: Error running {atoms}d: {e}")
        pass
    return 0

def main():
    compile_buckygen()

    processed_atoms = set()
    file_exists = os.path.isfile(CSV_FILE)
    if file_exists:
        try:
            with open(CSV_FILE, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row.get('Atoms'):
                        processed_atoms.add(int(row['Atoms']))
        except Exception as e:
            print(f"Warning: Could not read existing CSV for resume: {e}")

    # Open CSV in append mode
    with open(CSV_FILE, 'a', newline='') as f:
        fieldnames = ['Atoms', 'Total_IPR', 'Asymmetric_IPR', 'Percent_Asymmetric']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        if not file_exists or os.path.getsize(CSV_FILE) == 0:
            writer.writeheader()

        print(f"{'Atoms':>6} | {'Total IPR':>18} | {'Asymmetric IPR':>18} | {'% Asym':>8}")
        print("-" * 62)

        try:
            for atoms in range(START_ATOMS, END_ATOMS + 2, STEP):
                if atoms in processed_atoms:
                    continue 

                total = get_buckygen_count(atoms, ["-I"])
                nontrivial = get_buckygen_count(atoms, ["-I", "-V"])
                asymmetric = total - nontrivial
                pct = (asymmetric / total * 100) if total > 0 else 0
                
                writer.writerow({
                    'Atoms': atoms,
                    'Total_IPR': total,
                    'Asymmetric_IPR': asymmetric,
                    'Percent_Asymmetric': round(pct, 4)
                })
                f.flush() 

                print(f"{atoms:6} | {total:18,} | {asymmetric:18,} | {pct:7.2f}%")
                sys.stdout.flush()
                
        except KeyboardInterrupt:
            print("\nStopped by user. Progress saved.")

if __name__ == "__main__":
    main()
