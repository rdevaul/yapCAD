#!/bin/bash
# Thrust Structure Validation Demo
# =================================
#
# This script demonstrates the complete yapCAD validation workflow:
# 1. DSL syntax checking
# 2. Geometry generation and export
# 3. Mass budget validation (native yapCAD, no external deps)
# 4. FEA analysis (requires fenics/gmsh - optional)
#
# Usage:
#   ./run_demo.sh              # Full demo with FEA (requires fenics)
#   ./run_demo.sh --mass-only  # Mass check only (no external deps)
#   ./run_demo.sh --help       # Show usage

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
MASS_ONLY=""
DESIGN="optimized"
while [[ $# -gt 0 ]]; do
    case $1 in
        --mass-only|--skip-fea)
            MASS_ONLY="--mass-only"
            shift
            ;;
        --design)
            DESIGN="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --mass-only    Run mass check only (no FEA dependencies needed)"
            echo "  --design NAME  Design variant: baseline, light, optimized (default: optimized)"
            echo "  --help         Show this help message"
            echo ""
            echo "Requirements:"
            echo "  Mass check only: yapCAD with OCC BREP support"
            echo "  Full validation: conda install -c conda-forge fenics-dolfinx gmsh meshio"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "=============================================="
echo -e "  ${BLUE}Thrust Structure Validation Demo${NC}"
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
echo -e "${GREEN}OK${NC} DSL syntax valid"
echo ""

# Step 2: List available commands
echo -e "${YELLOW}Step 2: Available DSL commands:${NC}"
python3 -m yapcad.dsl list thrust_structure.dsl | grep -E "^  (MAKE|make)" || true
echo ""

# Step 3: Generate STEP files for CAD export
echo -e "${YELLOW}Step 3: Generating geometry exports...${NC}"
mkdir -p output

echo "   Generating baseline plate..."
python3 -m yapcad.dsl run thrust_structure.dsl MAKE_BASELINE_PLATE \
    --output output/baseline_plate.step 2>/dev/null
echo -e "   ${GREEN}OK${NC} output/baseline_plate.step"

echo "   Generating optimized plate..."
python3 -m yapcad.dsl run thrust_structure.dsl MAKE_OPTIMIZED_PLATE \
    --output output/optimized_plate.step 2>/dev/null
echo -e "   ${GREEN}OK${NC} output/optimized_plate.step"
echo ""

# Step 4: Run validation (mass check + optional FEA)
echo -e "${YELLOW}Step 4: Running validation...${NC}"
echo ""

if [ -n "$MASS_ONLY" ]; then
    echo "   (Running mass check only - use without --mass-only for FEA)"
    echo ""
fi

python3 run_validation.py --design "$DESIGN" $MASS_ONLY

# Summary
echo ""
echo "=============================================="
echo -e "  ${GREEN}Demo Complete!${NC}"
echo "=============================================="
echo ""
echo "Output files:"
echo "  Geometry:"
echo "    - output/baseline_plate.step"
echo "    - output/optimized_plate.step"
echo ""
echo "  Package:"
echo "    - thrust_structure.ycpkg/"
echo ""
if [ -z "$MASS_ONLY" ]; then
echo "  Visualization (ParaView):"
echo "    - thrust_structure.ycpkg/validation/results/thrust-fea/displacement.xdmf"
echo "    - thrust_structure.ycpkg/validation/results/thrust-fea/stress.xdmf"
echo ""
fi
echo "Next steps:"
echo "  - View STEP in CAD: FreeCAD, Fusion360, etc."
if [ -z "$MASS_ONLY" ]; then
echo "  - View FEA results: paraview thrust_structure.ycpkg/validation/results/thrust-fea/stress.xdmf"
fi
echo "  - Interactive view: python3 visualize.py --optimized"
echo ""
