# Historical Document

> **Note**: This is a historical planning document from October-December 2025. All seven
> phases of BREP integration have been completed. For current BREP implementation details,
> see `docs/yapBREP.rst`.

---

# BREP Integration Strategy

This document outlines a strategy for integrating a Boundary Representation (BREP) model into yapCAD, enabling support for STEP import and a more robust internal geometry representation. The primary goal is to achieve this while maintaining backward compatibility with existing yapCAD designs and packages.

## 1. Core Technology

We will leverage the **`pythonocc-core`** library, which provides Python bindings for the powerful OpenCascade Technology (OCC) modeling kernel.

*   **Why `pythonocc-core`?** It offers direct, low-level access to the OCCT kernel, providing the flexibility and control needed to build a robust BREP implementation within yapCAD. Its LGPL license is compatible with yapCAD's MIT license.
*   **Installation:** `pythonocc-core` will be added as a core dependency to the project.

## 2. Dependency Management

Due to the complex dependencies of `pythonocc-core`, we will use the `conda` package manager to create a dedicated environment for `yapCAD` development and usage. This approach allows us to leverage the pre-built binaries for `pythonocc-core` provided by the `conda-forge` channel, while still allowing `yapCAD` itself to be installed with `pip`.

An `environment.yml` file is provided in the root of the project to automate the creation of this environment. To create and activate the environment, run the following commands:

```bash
conda env create -f environment.yml
conda activate yapcad-brep
```

This will create a `conda` environment named `yapcad-brep` that includes `pythonocc-core` and other dependencies, and then install `yapCAD` in editable mode.

## 3. Data Model

A new set of classes will be introduced to represent the BREP topology. These classes will wrap the underlying `pythonocc-core` objects.

*   **`BrepSolid`**: Represents a solid body, composed of shells of faces.
*   **`BrepFace`**: Represents a single face of a solid, defined by an underlying analytical surface (e.g., a plane, cylinder, or NURBS surface) and bounded by a set of edges.
*   **`BrepEdge`**: Represents an edge, defined by an underlying analytical curve (e.g., a line, circle, or NURBS curve) and bounded by vertices.
*   **`BrepVertex`**: Represents a vertex, a point in 3D space.

These new classes will be defined in a new module, `src/yapcad/brep.py`.

## 4. Integration with the `Geometry` Class

The existing `yapcad.geometry.Geometry` class will be the primary mechanism for interacting with BREP objects.

*   The `Geometry` class's `__elem` attribute will be able to hold an instance of a `BrepSolid` (or other `Brep` objects).
*   A new type-checking function, `is_brep(obj)`, will be introduced to identify these new objects.
*   The `Geometry` constructor will be updated to recognize and wrap `Brep` objects.

## 5. Lazy Faceting

The current faceted representation will be generated "lazily" from the BREP model. This means that the triangular mesh will only be created when it is explicitly needed.

*   The `Geometry.surface()` method will be the trigger for lazy faceting.
*   When `surface()` is called on a `Geometry` object containing a `BrepSolid`, it will invoke a new `tessellate()` method on the `BrepSolid` object.
*   The `tessellate()` method will use `pythonocc-core`'s meshing capabilities to generate a triangular mesh of the BREP model.
*   The resulting faceted `['surface', vertices, normals, faces]` data structure will be cached within the `Geometry` object, just as it is for the current `poly2surface` workflow.

This approach ensures that the performance-intensive tessellation process is only performed when necessary, for example, for rendering or for exporting to STL.

## 6. Updating Generic Functions

The generic geometry functions in `geom.py` and `geom3d.py` will be updated to handle the new BREP types.

*   The `if/elif` chains in functions like `center()`, `bbox()`, and the transformation functions (`rotate()`, `translate()`, `scale()`) will be extended with `is_brep()` checks.
*   For transformations, the operations will be applied directly to the analytical BREP model using `pythonocc-core`'s transformation capabilities. This will ensure maximum precision.
*   For query functions like `center()` and `bbox()`, the information will be extracted from the BREP model using the appropriate `pythonocc-core` functions.

## 7. Backward Compatibility

To ensure backward compatibility, the new BREP classes will provide a mechanism to convert their representations into the legacy `geomlist` format.

*   A `to_geomlist()` method will be added to the `BrepSolid`, `BrepFace`, and `BrepEdge` classes.
*   This method will extract the boundary curves of the BREP objects and convert them into a `geomlist` of `yapCAD`'s primitive `line`, `arc`, and `spline` objects.
*   This will allow the new BREP objects to be used with existing functions that expect the old `geomlist` or `poly` formats, such as `intersectXY()`.

## 8. STEP Import

The initial STEP importer now lives in `src/yapcad/io/step_importer.py`.

*   `import_step(path)` uses `pythonocc-core` to read a STEP file (run inside `conda activate yapcad-brep`) and returns a list of `Geometry(BrepSolid)` wrappers.
*   The current implementation traverses solids inside the file and generates lazy-faceted geometries on demand. Mapping the full face/edge graph is still on the roadmap.
*   Example usage:
    ```python
    from yapcad.io.step_importer import import_step

    parts = import_step("rocket.step")
    for geom in parts:
        print(geom.bbox)
    ```

## 9. Development Roadmap

The implementation of this strategy will be broken down into the following phases:

1.  **Phase 1: Core BREP Integration:** ✅ COMPLETE
    *   Add `pythonocc-core` as a dependency. ✅
    *   Implement the core `Brep` classes in `src/yapcad/brep.py`. ✅
    *   Update the `Geometry` class to recognize and wrap `Brep` objects. ✅
    *   Implement the `tessellate()` method and integrate it with the `Geometry.surface()` method for lazy faceting. ✅
2.  **Phase 2: Generic Function Updates:** ✅ COMPLETE
    *   Update the generic geometry functions in `geom.py` and `geom3d.py` to handle the new `Brep` types. ✅
3.  **Phase 3: Backward Compatibility:** ✅ COMPLETE
    *   Implement the `to_geomlist()` method on the `Brep` classes. ✅
4.  **Phase 4: STEP Import:** ✅ COMPLETE
    *   Implement the `import_step()` function in `src/yapcad/io/step_importer.py`. ✅
    *   Serialize OCC shapes directly inside geometry JSON/`.ycpkg` metadata so analytic models round-trip without tessellation loss. ✅ *(Implemented via base64-encoded `.brep` payloads in the `metadata.brep` block.)*
    *   Provide an OCC-backed boolean engine (`--engine occ`) that consumes the stored BREP metadata and emits new analytic shapes. ✅ *(Implemented; falls back to tessellation engines when BREP data is missing.)*
5.  **Phase 5: Testing:** ✅ COMPLETE
    *   Create a comprehensive suite of tests for the new BREP functionality, including tests for STEP import, lazy faceting, and backward compatibility. ✅
6.  **Phase 6: Native BREP Representation:** ✅ COMPLETE
    *   Implement native BREP topology hierarchy in `src/yapcad/native_brep.py`. ✅
        -   `brep_vertex`, `brep_edge`, `brep_trim`, `brep_loop`, `brep_face`, `brep_shell`, `brep_solid`
        -   `TopologyGraph` class for managing entity relationships
    *   Support for curve types: line, circle/arc, B-spline/NURBS. ✅
    *   Support for surface types: plane, sphere, cylinder, cone, torus, B-spline/NURBS. ✅
    *   Serialization/deserialization to/from JSON. ✅
    *   Transformation support (translate, rotate, scale). ✅
7.  **Phase 7: Bidirectional OCC Conversion:** ✅ COMPLETE
    *   OCC → Native BREP conversion in `src/yapcad/occ_native_convert.py`. ✅
        -   `occ_surface_to_native()`: Convert OCC faces to native analytic surfaces
        -   `occ_edge_to_native()`: Convert OCC edges to native curve edges
        -   `occ_solid_to_native_brep()`: Full topology conversion preserving structure
    *   Native BREP → OCC conversion. ✅
        -   `native_vertex_to_occ()`, `native_edge_to_occ()`, `native_loop_to_occ_wire()`
        -   `native_face_to_occ()`, `native_brep_to_occ()`
    *   Round-trip fidelity with volume preservation:
        -   Box: 100% ✅
        -   Cylinder: 100% ✅
        -   Sphere: 100% ✅
        -   Cone: 100% ✅
    *   45 native BREP tests passing ✅

This phased approach has allowed for the gradual and controlled integration of BREP support into `yapCAD`, minimizing disruption to the existing codebase and ensuring that the new functionality is robust and well-tested. **All phases are now complete.**

## 10. Future Enhancements

Potential future work includes:
*   **Partial surface support**: Handle trimmed/partial surfaces in round-trip conversion
*   **IGES import**: Add support for IGES file format alongside STEP
*   **Boolean operation result tracking**: Preserve face/edge provenance through boolean operations
*   ~~**Advanced curve types**: Ellipse, hyperbola, parabola support in edge conversion~~ ✅ **COMPLETE** (November 2025) - Parabola and hyperbola primitives implemented in `native_brep.py`
