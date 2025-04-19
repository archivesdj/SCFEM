#!/bin/bash

# Exit on error
set -e

# Store the current directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Generate mesh using GMSH
echo "Generating mesh..."
cd "$SCRIPT_DIR"
gmsh -3 -format msh2 -order 1 -optimize_netgen -save_all cube.geo

# Run the solver
echo "Running solver..."
cd "$SCRIPT_DIR/.."
dotnet run --project SCFEM.CLI/SCFEM.CLI.csproj examples/cube.msh examples/materials.txt 