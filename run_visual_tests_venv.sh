#!/bin/bash
# Wrapper script to run visual tests from within the virtual environment

# Check if virtual environment is activated
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Activating virtual environment..."
    if [ -f "v_312/bin/activate" ]; then
        source v_312/bin/activate
    else
        echo "Error: v_312/bin/activate not found. Please create the virtual environment first."
        exit 1
    fi
fi

# Set PYTHONPATH
export PYTHONPATH=./src

# Run the visual test runner
python run_visual_tests.py "$@"