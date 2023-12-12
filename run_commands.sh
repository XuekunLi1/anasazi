#!/bin/bash

# Stop the script if any command fails
set -e

# Step 1: make clean
echo "Running 'make clean'..."
make clean

# Step 2: make all
echo "Running 'make all'..."
make all

# Step 3: Run the C++ program using mpirun
echo "Running the C++ program..."
mpirun -n 1 bin/main.exe props/config.props props/model.props

# Step 4: Run the Python script
echo "Running the Python script..."
python3 main_group.py

echo "All commands executed successfully."

