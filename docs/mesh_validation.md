# Mesh Validation Pipeline

This repository now includes a small command–line helper that mirrors the
checks typically performed downstream of a CAD/boolean workflow.  The goal is
to gauge whether the generated STL will survive real fabrication tools even if
the mesh is not mathematically watertight.

## Prerequisites

The script looks for the following utilities on `PATH`:

| Tool     | Purpose                              | macOS install hint |
|----------|--------------------------------------|--------------------|
| `admesh` | STL statistics (holes, volume, etc.) | `brew install admesh` |
| `meshfix`| Mesh repair via CGAL                 | `brew install meshfix` |

If you prefer to stay in Python space, install the optional
[`pymeshfix`](https://pypi.org/project/PyMeshFix/) and
[`trimesh`](https://pypi.org/project/trimesh/) packages; the validation script
will automatically fall back to those when the stand-alone `meshfix` binary is
missing.

You can pull these extras in with a single command using the optional-deps
group defined in `pyproject.toml`:

```bash
pip install yapCAD[meshcheck]
```
| `prusa-slicer` / `prusa-slicer-console` / `slic3r` | CLI slicer check  | download from vendor; add binary to `PATH` |

If a tool is missing the script will skip the corresponding step and include a
note in the report.

## Usage

1. Export an STL using the demo driver (or your own code):

   ```bash
   PYTHONPATH=src python examples/solid_boolean_demo.py \
       --mode stl --operation difference --shapes box_hole \
       --output box_hole_difference
   ```

2. Run the validator:

   ```bash
   python tools/validate_mesh.py box_hole_difference.stl \
       --workdir build/mesh_checks
   ```

   This prints a summary to stdout and keeps any repaired STL/G-code files
   under `build/mesh_checks`.

3. (Optional) Request JSON output for CI:

   ```bash
   python tools/validate_mesh.py box_hole_difference.stl --json
   ```

The report includes highlights from `admesh`, the return code from `meshfix`
(with output file location), and the status of the slicer run.  When the slicer
succeeds you’ll find a G-code file beside the repaired meshes.

This pipeline gives a far more practical read on model quality than the strict
`issolidclosed` check, and mirrors the tools typically invoked before sending a
part to printers or CAM software.
