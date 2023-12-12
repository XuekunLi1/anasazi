#!/bin/bash

# Specify the CSV file and model.props file
csv_file="matlab_bestlh_matrix.csv"
props_file="props/model.props"

# Check if CSV and props files exist
if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file '$csv_file' not found."
    exit 1
fi

if [ ! -f "$props_file" ]; then
    echo "Error: Properties file '$props_file' not found."
    exit 1
fi

# Process each line of the CSV file
while IFS=',' read -r corn1 corn2 corn3 corn4 corn5; do
    # Update the model.props file
    sed -i "s/^\(max\.fission\.age\s*=\s*\).*\$/\1$corn1/" $props_file
    sed -i "s/^\(max\.death\.age\s*=\s*\).*\$/\1$corn2/" $props_file
    sed -i "s/^\(annual\.variance\s*=\s*\).*\$/\1$corn3/" $props_file
    sed -i "s/^\(fertility\.prop\s*=\s*\).*\$/\1$corn4/" $props_file
    sed -i "s/^\(harvest\.adj\s*=\s*\).*\$/\1$corn5/" $props_file

    if [ $? -ne 0 ]; then
        echo "Error updating $props_file"
        exit 1
    fi

    # Run the mpirun command
    echo "Running the C++ program..."
    mpirun -n 1 bin/main.exe props/config.props props/model.props
    if [ $? -ne 0 ]; then
        echo "Error executing mpirun command"
        exit 1
    fi

    # Run the python script
    echo "Running the C++ program..."
    python3 main_group.py
    if [ $? -ne 0 ]; then
        echo "Error executing Python script"
        exit 1
    fi

    echo "Processing done for one row."
done < "$csv_file"

echo "All rows processed successfully."

