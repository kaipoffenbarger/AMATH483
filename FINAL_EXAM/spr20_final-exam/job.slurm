#!/bin/bash
#SBATCH -A amath
#SBATCH -N 4
#SBATCH --job-name=final_problem-4
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --partition=cpu-g2
#SBATCH --time=00:45:00
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kaipoff@uw.edu

module load gcc intel/oneAPI/2023.2.1

N=4096

processes=( 1 2 4 8 16 )
ps=(1 1 2 2 4)
qs=(1 2 2 4 4)

output_csv="output_times.csv"
echo 'N,processes,trial,p,q,wall_time,cpu_time' > "$output_csv"

source cmpl.krazeb

for i in "${!processes[@]}"; do
    current_processes=${processes[i]}
    current_p=${ps[i]}
    current_q=${qs[i]}

    # Construct the command as an array for robustness
    # This avoids 'eval' and handles arguments with spaces gracefully
    cmd_args=(mpirun -np "$current_processes" ./xkrazeb "$current_p" "$current_q" 48 48 "$N")

    for trial in {1..6}; do
        echo "Running: ${cmd_args[*]} (Trial $trial)"

        # Capture output (stdout and stderr)
        output=$("${cmd_args[@]}" 2>&1)

        # Print the program's full output with indentation (optional, for debugging/display)
        echo "$output" | sed 's/^/  /'

        # Extract wall time
        wall_time=$(echo "$output" | grep -oP 'wall_t: \K[0-9.]+')

        # Extract cpu time
        cpu_time=$(echo "$output" | grep -oP 'cpu_t: \K[0-9.]+')

        # Basic error handling for wall time extraction
        if [ -z "$wall_time" ]; then
            echo "WARNING: Wall time not found for P=$current_processes, trial=$trial. Output excerpt:" >&2
            echo "$output" | head -n 5 >&2 # Show first 5 lines of output for context
            wall_time="N/A" # Use a placeholder
        fi
        # Basic error handling for cpu time extraction
        if [ -z "$wall_time" ]; then
            echo "WARNING: CPU time not found for P=$current_processes, trial=$trial. Output excerpt:" >&2
            echo "$output" | head -n 5 >&2 # Show first 5 lines of output for context
            cpu_time="N/A" # Use a placeholder
        fi

        # Output to CSV
        echo "$N,$current_processes,$trial,$current_p,$current_q,$wall_time,$cpu_time" >> "$output_csv"
    done
    echo "---" # Separator between process configurations
done

echo "Results saved to $output_csv"
