# Repository Guidelines

## Project Structure & Module Organization
Core library code lives in `src/yapcad`, with geometry primitives (`geom.py`, `geom3d.py`), boolean ops (`boolean/`), IO adapters (`io/`), and render helpers (`drawable.py`, `ezdxf_drawable.py`). Vendored helpers live under `src/yapcad/contrib`, with corresponding license notices in `third_party/`. PyPI metadata sits alongside in `pyproject.toml` and `setup.cfg`. Automated tests are under `tests/`, mirroring the module names (e.g., `test_geom3d.py`). Sphinx docs reside in `docs/`, while `examples/`, `dxf/`, and `images/` provide sample CAD assets for manual inspection.

## Setup, Build & Test Commands
Create a local environment with:
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .[tests]
```
Run the full suite (includes coverage) via `pytest`. Targeted runs, such as `pytest tests/test_geom3d.py`, keep iteration fast. Update geometry visual baselines with `python run_visual_tests.py --update`, then review artifacts in `build/`. Build the distribution wheel with `python -m build` and documentation with `sphinx-build docs build/sphinx`.

## Coding Style & Naming Conventions
Follow PEP 8 with 4-space indentation and descriptive, lowercase module names. Public API classes use `CamelCase`, functions and variables stay `snake_case`. Prefer explicit imports from sibling modules. Keep geometry helper names aligned with their dimensionality (e.g., `*_3d`). Run `flake8` when touching complex files to stay consistent with `setup.cfg`.

## Testing Guidelines
pytest discovers files named `test_*.py` inside `tests/`; mirror production module names to keep coverage readable. Maintain coverage at or above the default `--cov yapcad --cov-report term-missing` threshold before opening PRs. Use parametrized tests for shape variants and add integration checks for new IO handlers in `tests/test_io_*.py`. When adding CAD fixtures, store them under `examples/` and reference them relative to the repo root.

## Commit & Pull Request Guidelines
Commit messages follow an imperative summary (`Add 3D sweep helper`), optionally tagging release steps (`Prepare v0.6.0 release`). Squash noisy WIP commits locally. For pull requests, include: purpose summary, testing commands executed, links to relevant issues, and before/after renders if geometry changes affect visuals. Ensure CI (pytest + docs) passes prior to requesting review and note any follow-up tasks explicitly.
