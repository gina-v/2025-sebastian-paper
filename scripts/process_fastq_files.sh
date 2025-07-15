#!/bin/bash

# Check if at least one file is provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 threads output_file file1 [file2 [file3]]"
    exit 1
fi

threads=$1

output_file=$2

# Shift the arguments so that $@ contains only the files
shift 2

# Run the seqkit stats command and store the output
output=$(seqkit stats -j "$threads" "$@")


# Extract columns for each file
for i in $(seq 1 $#); do
    # Use awk to extract the line corresponding to the ith file
    line=$(echo "$output" | awk -v line_num=$((i + 1)) 'NR==line_num {print $1, $4, $5}')

    # Parse columns from the line
    read -r file read_count base_count <<< "$line"

    # Remove commas from numbers
    read_count=$(echo "$read_count" | tr -d ',')
    base_count=$(echo "$base_count" | tr -d ',')

    if [[ -f "$output_file" ]]; then
        if grep -q "$file" "$output_file"; then
            echo "Skipping $file as it is already in $output_file"
        else
            # Print or use the extracted information
            echo "Input file: $file"
            echo "Number of sequences: $read_count"
            echo "Number of bases or residues: $base_count"
            echo "-----"

            # Append the data to the output file
            echo -e "$file,$read_count,$base_count" >> "$output_file"
        fi
    else
     # If the file doesn't exist, create it and add the first entry
        echo "Input file: $file"
        echo "Number of sequences: $read_count"
        echo "Number of bases or residues: $base_count"
        echo "-----"

        echo -e "$file,$read_count,$base_count" > "$output_file"
    fi
done

echo "Results written to $output_file"
