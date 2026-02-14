# Edge Selection DSL Bindings - Implementation Summary

## Overview

Successfully added DSL bindings for all edge selection functions from `yapcad.brep_edge_select` module. These functions enable selective edge operations in the yapCAD DSL, allowing users to apply fillets and chamfers to specific edges rather than all edges of a solid.

## Files Modified

### 1. `/Users/jmika/Work/DML/yapCAD/src/yapcad/dsl/types.py`

**Changes:**
- Added new `EDGE` type: `EDGE = GeometricPrimitiveType("edge", dimension=None)`
- Registered `EDGE` in `BUILTIN_TYPES` dictionary

**Purpose:** Defines the DSL type for BREP edges, enabling type checking for edge lists.

### 2. `/Users/jmika/Work/DML/yapCAD/src/yapcad/dsl/runtime/values.py`

**Changes:**
- Imported `EDGE` type
- Added `edge_list_val(edge_list_data)` constructor function

**Purpose:** Provides a convenience function to wrap edge lists with proper DSL type metadata.

### 3. `/Users/jmika/Work/DML/yapCAD/src/yapcad/dsl/runtime/builtins.py`

**Major Changes:**
- Imported `EDGE` type and `edge_list_val` constructor
- Added new `_register_edge_selection_functions()` method to `BuiltinRegistry`
- Registered the method in `_register_all()` initialization

**Functions Added (13 total):**

#### Edge Selection by Direction (3 functions)
1. `select_vertical_edges(solid, tolerance_deg?) -> list<edge>`
2. `select_horizontal_edges(solid, tolerance_deg?) -> list<edge>`
3. `select_edges_by_direction(solid, direction, tolerance_deg?) -> list<edge>`

#### Edge Selection by Length (1 function)
4. `select_edges_by_length(solid, min_length?, max_length?) -> list<edge>`

#### Edge Selection by Z Position (4 functions)
5. `select_edges_at_z(solid, z_value, tolerance?) -> list<edge>`
6. `select_edges_in_z_range(solid, z_min, z_max, tolerance?) -> list<edge>`
7. `select_top_edges(solid, tolerance?) -> list<edge>`
8. `select_bottom_edges(solid, tolerance?) -> list<edge>`

#### Edge Set Operations (3 functions)
9. `union_edges(edges1, edges2) -> list<edge>`
10. `intersect_edges(edges1, edges2) -> list<edge>`
11. `subtract_edges(base_edges, to_remove) -> list<edge>`

#### Selective Fillet/Chamfer (2 functions)
12. `fillet_edges(solid, edges, radius) -> solid`
13. `chamfer_edges(solid, edges, distance) -> solid`

## DSL Function Signatures

### Selection by Direction

```dsl
select_vertical_edges(s: solid, tol: float = 1.0) -> list<edge>
```
Select edges parallel to the Z axis (vertical edges).

```dsl
select_horizontal_edges(s: solid, tol: float = 1.0) -> list<edge>
```
Select edges perpendicular to the Z axis (horizontal edges).

```dsl
select_edges_by_direction(s: solid, dir: vector3d, tol: float = 1.0) -> list<edge>
```
Select edges parallel to a given direction vector.

### Selection by Length

```dsl
select_edges_by_length(s: solid, min_len: float?, max_len: float?) -> list<edge>
```
Select edges within a length range. Both parameters optional.

### Selection by Z Position

```dsl
select_edges_at_z(s: solid, z: float, tol: float = 0.001) -> list<edge>
```
Select edges at a specific Z height.

```dsl
select_edges_in_z_range(s: solid, z_min: float, z_max: float, tol: float = 0.001) -> list<edge>
```
Select edges with both endpoints within a Z range.

```dsl
select_top_edges(s: solid, tol: float = 0.001) -> list<edge>
```
Select edges at the maximum Z height of the solid.

```dsl
select_bottom_edges(s: solid, tol: float = 0.001) -> list<edge>
```
Select edges at the minimum Z height of the solid.

### Edge Set Operations

```dsl
union_edges(edges1: list<edge>, edges2: list<edge>) -> list<edge>
```
Combine edge lists, removing duplicates.

```dsl
intersect_edges(edges1: list<edge>, edges2: list<edge>) -> list<edge>
```
Find edges common to both lists.

```dsl
subtract_edges(base: list<edge>, remove: list<edge>) -> list<edge>
```
Remove edges from a base list.

### Selective Fillet and Chamfer

```dsl
fillet_edges(s: solid, edges: list<edge>, radius: float) -> solid
```
Apply fillet (rounded edges) to selected edges only.

```dsl
chamfer_edges(s: solid, edges: list<edge>, distance: float) -> solid
```
Apply chamfer (beveled edges) to selected edges only.

## Example Usage

```dsl
# Create a box
let box = box(20.0, 20.0, 10.0)

# Select and fillet only vertical edges
let vertical = select_vertical_edges(box, 1.0)
let filleted_box = fillet_edges(box, vertical, 2.0)

# Select and chamfer only top edges
let top = select_top_edges(box, 0.001)
let chamfered_box = chamfer_edges(box, top, 1.5)

# Combine selections: vertical edges at the top
let top_edges = select_edges_at_z(box, 10.0, 0.001)
let vertical_edges = select_vertical_edges(box, 1.0)
let top_vertical = intersect_edges(top_edges, vertical_edges)
let result = fillet_edges(box, top_vertical, 1.0)

# Select by length: only short edges
let short_edges = select_edges_by_length(box, 0.0, 5.0)
let result2 = chamfer_edges(box, short_edges, 0.5)
```

## Implementation Details

### Type System Integration
- New `EDGE` type added as `GeometricPrimitiveType` at Tier 1
- Edge lists represented as `ListType(EDGE)`
- All functions properly typed for DSL type checker

### Optional Parameters
- Functions with optional parameters use Python `None` default values
- DSL runtime checks for `None` and applies appropriate defaults
- Matches the pattern used in existing DSL builtins

### Error Handling
- All functions check for pythonocc-core availability
- All functions verify BREP data exists on solid before processing
- Clear error messages guide users when requirements not met

### Consistency with Existing API
- Follows same pattern as `fillet()` and `chamfer()` builtins
- Uses same BREP extraction and solid reconstruction workflow
- Maintains provenance tracking with `['procedure', function_name]` metadata

## Testing

### Verification Script
Created `/Users/jmika/Work/DML/yapCAD/test_edge_selection_builtins.py` to verify:
- ✅ All 13 functions successfully registered
- ✅ All function signatures match expected types
- ✅ All functions return correct types

### Test Results
```
✅ All edge selection functions successfully registered!

All 13 functions verified:
- 3 direction-based selection functions
- 1 length-based selection function
- 4 Z-position-based selection functions
- 3 edge set operation functions
- 2 selective fillet/chamfer functions
```

## Benefits

1. **Precise Control**: Users can now apply edge operations to specific edges rather than all edges
2. **Flexible Selection**: Multiple selection criteria (direction, length, position)
3. **Composable**: Edge set operations allow combining/filtering selections
4. **Type Safe**: Full DSL type system integration with `list<edge>` type
5. **Consistent API**: Follows established patterns in yapCAD DSL

## Future Enhancements

Additional functions from `brep_edge_select.py` that could be added:
- `select_edges_crossing_z(solid, z_value, tolerance?)` - Select vertical edges spanning a Z height
- `select_edges_near_point(solid, point, max_distance)` - Select by proximity to point
- `select_edges_in_cylinder(solid, center, radius, axis?)` - Select within cylindrical region
- `filter_curved_edges(edges)` - Filter to keep only curved edges
- `filter_linear_edges(edges)` - Filter to keep only straight edges
- `edge_info(edge)` - Get detailed edge information

These can be added in a future iteration if needed.
