#!/bin/bash
# Thrust Structure FEA Demo
# =========================
#
# This script demonstrates the complete yapCAD workflow:
# 1. DSL syntax checking
# 2. Baseline geometry generation
# 3. Design comparison / optimization
# 4. Export to manufacturing formats

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=============================================="
echo "  Thrust Structure FEA Demonstration"
echo "=============================================="
echo ""

# Set up Python path
export PYTHONPATH="$SCRIPT_DIR/../../src:$PYTHONPATH"

# Check for yapCAD
if ! python3 -c "import yapcad" 2>/dev/null; then
    echo -e "${RED}Error: yapCAD not found in Python path${NC}"
    echo "Run from yapCAD environment or set PYTHONPATH"
    exit 1
fi

# Step 1: Check DSL syntax
echo -e "${YELLOW}Step 1: Checking DSL syntax...${NC}"
python3 -m yapcad.dsl check thrust_structure.dsl
echo -e "${GREEN}✓ DSL syntax OK${NC}"
echo ""

# Step 2: List available commands
echo -e "${YELLOW}Step 2: Available DSL commands:${NC}"
python3 -m yapcad.dsl list thrust_structure.dsl | grep -E "^  (MAKE|make)" || true
echo ""

# Step 3: Generate baseline geometry
echo -e "${YELLOW}Step 3: Generating baseline geometry...${NC}"
mkdir -p output
python3 -m yapcad.dsl run thrust_structure.dsl MAKE_BASELINE_PLATE \
    --output output/baseline_plate.step
echo -e "${GREEN}✓ Baseline exported to output/baseline_plate.step${NC}"
echo ""

# Step 4: Show baseline info
echo -e "${YELLOW}Step 4: Baseline design information:${NC}"
python3 visualize.py --baseline --info --no-viewer
echo ""

# Step 5: Run design comparison
echo -e "${YELLOW}Step 5: Running design comparison...${NC}"
python3 optimize.py --output output/optimization
echo ""

# Step 6: Generate optimized geometry
echo -e "${YELLOW}Step 6: Generating optimized geometry...${NC}"
python3 -m yapcad.dsl run thrust_structure.dsl MAKE_OPTIMIZED_PLATE \
    --output output/optimized_plate.step
echo -e "${GREEN}✓ Optimized design exported to output/optimized_plate.step${NC}"
echo ""

# Step 7: Show optimized info
echo -e "${YELLOW}Step 7: Optimized design information:${NC}"
python3 visualize.py --optimized --info --no-viewer
echo ""

# Summary
echo "=============================================="
echo -e "${GREEN}  Demo Complete!${NC}"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - output/baseline_plate.step"
echo "  - output/optimized_plate.step"
echo "  - output/optimization/comparison_results.json"
echo "  - output/optimization/best_thrust_plate.step"
echo ""
echo "Next steps:"
echo "  - View in CAD: Import .step files into FreeCAD, Fusion360, etc."
echo "  - Interactive view: python3 visualize.py --optimized"
echo ""
