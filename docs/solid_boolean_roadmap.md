# Boolean Engine Roadmap

## Status as of v0.5.1 (January 2025)

### Completed Work

- ✅ **Native Boolean Engine**: Production-ready implementation in `yapcad/boolean/native.py`
  - Robust normal orientation for all operations (union, intersection, difference)
  - Quality-based filtering for degenerate triangles (aspect ratio checks)
  - Containment-based filtering to eliminate interior overlap artifacts
  - All primitive boolean operations validated with correct watertight geometry

- ✅ **Code Organization**: Clean separation of boolean operations
  - Native engine in `yapcad/boolean/native.py`
  - Backward-compatible re-exports from `geom3d`
  - Modular architecture supports multiple backends

- ✅ **Engine Selector UX**: Flexible backend selection
  - `solid_boolean(..., engine='native')` API
  - Support for `trimesh:manifold` and `trimesh:blender` backends
  - Environment variables: `YAPCAD_BOOLEAN_ENGINE`, `YAPCAD_TRIMESH_BACKEND`
  - Command-line demo with `--engine` flag

- ✅ **Validation Infrastructure**: Comprehensive mesh quality checks
  - `tools/validate_mesh.py` CLI tool
  - Integration with `admesh`, `meshfix`, and slicers
  - 106 tests passing including boolean regression suite

### In Progress

- [ ] **Benchmark Harness**: Systematic comparison across engines
  - Need automated iteration across all registered engines
  - STL output comparison and quality metrics
  - Performance profiling for complex operations

- [ ] **External Engine Refinement**: Production deployment guidance
  - Document installation requirements for each backend
  - Best practices for engine selection
  - Tolerance mapping for different engines

### Near-Term Goals

- [ ] **Documentation Expansion**: Complete API reference
  - Engine selection guide with decision matrix
  - Performance characteristics and trade-offs
  - Migration guide from 0.5.0 to 0.5.1

- [ ] **Metrics Collection**: Engine-specific reporting
  - Shell count, boundary edge detection
  - Validation scores in machine-readable format (JSON)
  - Automated regression dashboards

### Open Questions

- **Engine Strategy**: Keep native engine as equal option vs fallback?
  - **Current stance**: Native engine is production-ready for most use cases
  - External engines provide alternatives for specific requirements (licensing, performance)
  - Recommend native for offline/embedded use, external for high-complexity operations

- **Quality Metrics**: What should every engine report?
  - Proposed minimum: triangle count, boundary edge count, watertightness status
  - Optional: volume calculation, self-intersection checks, performance timing

- **Hybrid Workflows**: When to combine engines?
  - Native preprocessing + external boolean + native cleanup?
  - Needs experimentation with real-world complex geometries

### Future Work (v0.6.0 and beyond)

- Automatic tolerance scaling per engine
- JSON-formatted validation output for CI/CD integration
- Performance optimization for large assemblies
- Enhanced STEP export with analytical surfaces (not just tessellated)
