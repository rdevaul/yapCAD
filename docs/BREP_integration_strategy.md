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

A new module, `src/yapcad/io/step_importer.py`, will be created to handle STEP import.

*   This module will contain a new function, `import_step(file_path)`, which will use `pythonocc-core` to read a STEP file.
*   The `import_step()` function will traverse the imported STEP model and create corresponding `yapCAD` `BrepSolid`, `BrepFace`, `BrepEdge`, and `BrepVertex` objects.
*   The function will return a `yapCAD` `Geometry` object (or a list of `Geometry` objects) containing the imported BREP model.

## 9. Development Roadmap

The implementation of this strategy will be broken down into the following phases:

1.  **Phase 1: Core BREP Integration:**
    *   Add `pythonocc-core` as a dependency.
    *   Implement the core `Brep` classes in `src/yapcad/brep.py`.
    *   Update the `Geometry` class to recognize and wrap `Brep` objects.
    *   Implement the `tessellate()` method and integrate it with the `Geometry.surface()` method for lazy faceting.
2.  **Phase 2: Generic Function Updates:**
    *   Update the generic geometry functions in `geom.py` and `geom3d.py` to handle the new `Brep` types.
3.  **Phase 3: Backward Compatibility:**
    *   Implement the `to_geomlist()` method on the `Brep` classes.
4.  **Phase 4: STEP Import:**
    *   Implement the `import_step()` function in `src/yapcad/io/step_importer.py`.
5.  **Phase 5: Testing:**
    *   Create a comprehensive suite of tests for the new BREP functionality, including tests for STEP import, lazy faceting, and backward compatibility.

This phased approach will allow for the gradual and controlled integration of BREP support into `yapCAD`, minimizing disruption to the existing codebase and ensuring that the new functionality is robust and well-tested.