"""STEP export utilities for yapCAD surfaces and solids."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, TextIO, Tuple

from yapcad.geometry_utils import Triangle, triangles_from_mesh
from yapcad.mesh import mesh_view

Vec3 = Tuple[float, float, float]


@dataclass(frozen=True)
class _Vertex:
    coords: Vec3
    point_id: int
    vertex_id: int


@dataclass(frozen=True)
class _EdgeData:
    edge_curve: int
    oriented_forward: int
    oriented_reverse: int
    start_vertex: int
    end_vertex: int


class _EntityWriter:
    """Helper to append STEP entities with sequential ids."""

    def __init__(self) -> None:
        self.entities: List[str] = []

    def add(self, record: str) -> int:
        idx = len(self.entities) + 1
        self.entities.append(f"#{idx} = {record};")
        return idx

    def write(self, stream: TextIO) -> None:
        for entity in self.entities:
            stream.write(entity + "\n")


def write_step(obj: Sequence,
               path_or_file,
               *,
               name: str = 'yapCAD',
               schema: str = 'AUTOMOTIVE_DESIGN_CC2') -> None:
    """Export ``obj`` (surface or solid) to a STEP file using a faceted BREP."""

    triangles = list(triangles_from_mesh(mesh_view(obj)))
    writer = _EntityWriter()

    vertex_map: Dict[Vec3, _Vertex] = {}

    for tri in triangles:
        for coords in (tri.v0, tri.v1, tri.v2):
            if coords not in vertex_map:
                pt = writer.add(
                    f"CARTESIAN_POINT('', ({coords[0]:.6f}, {coords[1]:.6f}, {coords[2]:.6f}))"
                )
                vert = writer.add(f"VERTEX_POINT('', #{pt})")
                vertex_map[coords] = _Vertex(coords, pt, vert)

    def _line_for_edge(start: _Vertex, end: _Vertex) -> Tuple[int, int, int]:
        direction_vec = _normalize((end.coords[0] - start.coords[0],
                                    end.coords[1] - start.coords[1],
                                    end.coords[2] - start.coords[2]))
        direction_id = writer.add(
            f"DIRECTION('', ({direction_vec[0]:.6f}, {direction_vec[1]:.6f}, {direction_vec[2]:.6f}))"
        )
        length = _distance(start.coords, end.coords)
        vector_id = writer.add(f"VECTOR('', #{direction_id}, {length:.6f})")
        line_id = writer.add(f"LINE('', #{start.point_id}, #{vector_id})")
        return direction_id, vector_id, line_id

    edge_store: Dict[frozenset, _EdgeData] = {}
    face_ids: List[int] = []

    for tri in triangles:
        v0 = vertex_map[tri.v0]
        v1 = vertex_map[tri.v1]
        v2 = vertex_map[tri.v2]

        loop_edges: List[int] = []
        for start, end in ((v0, v1), (v1, v2), (v2, v0)):
            key = frozenset({start.coords, end.coords})
            data = edge_store.get(key)
            if data is None:
                _, _, line_id = _line_for_edge(start, end)
                edge_curve = writer.add(
                    f"EDGE_CURVE('', #{start.vertex_id}, #{end.vertex_id}, #{line_id}, .T.)"
                )
                oriented_forward = writer.add(
                    f"ORIENTED_EDGE('', *, *, #{edge_curve}, .T.)"
                )
                oriented_reverse = writer.add(
                    f"ORIENTED_EDGE('', *, *, #{edge_curve}, .F.)"
                )
                data = _EdgeData(edge_curve, oriented_forward, oriented_reverse,
                                 start.vertex_id, end.vertex_id)
                edge_store[key] = data

            if data.start_vertex == start.vertex_id and data.end_vertex == end.vertex_id:
                loop_edges.append(data.oriented_forward)
            else:
                loop_edges.append(data.oriented_reverse)

        edge_loop = writer.add(
            "EDGE_LOOP('', (" + ", ".join(f"#{eid}" for eid in loop_edges) + "))"
        )
        face_outer = writer.add(f"FACE_OUTER_BOUND('', #{edge_loop}, .T.)")

        normal = _normalize(tri.normal)
        ref_vec = _perpendicular(normal)
        normal_id = writer.add(
            f"DIRECTION('', ({normal[0]:.6f}, {normal[1]:.6f}, {normal[2]:.6f}))"
        )
        ref_dir = writer.add(
            f"DIRECTION('', ({ref_vec[0]:.6f}, {ref_vec[1]:.6f}, {ref_vec[2]:.6f}))"
        )
        axis = writer.add(f"AXIS2_PLACEMENT_3D('', #{v0.point_id}, #{normal_id}, #{ref_dir})")
        plane = writer.add(f"PLANE('', #{axis})")
        face = writer.add(f"ADVANCED_FACE('', (#{face_outer}), #{plane}, .T.)")
        face_ids.append(face)

    shell_id = writer.add(
        "CLOSED_SHELL('', (" + ", ".join(f"#{face}" for face in face_ids) + "))"
    )
    brep_id = writer.add(f"MANIFOLD_SOLID_BREP('{name}', #{shell_id})")

    # Representation contexts and product definitions
    app_context = writer.add("APPLICATION_CONTEXT('mechanical design')")
    prod_context = writer.add(f"PRODUCT_CONTEXT('part definition','mechanical', #{app_context})")
    product = writer.add(f"PRODUCT('{name}', '{name}', '', (#{prod_context}))")
    formation = writer.add(f"PRODUCT_DEFINITION_FORMATION('', '', #{product})")
    pd_context = writer.add(f"PRODUCT_DEFINITION_CONTEXT('design', 'mechanical', #{app_context})")
    product_def = writer.add(f"PRODUCT_DEFINITION('', #{formation}, #{pd_context})")
    shape = writer.add(f"PRODUCT_DEFINITION_SHAPE('', '', #{product_def})")

    length_unit = writer.add("( LENGTH_UNIT() AND NAMED_UNIT(*) AND SI_UNIT(.MILLI., .METRE.) )")
    plane_angle_unit = writer.add("( NAMED_UNIT(*) AND PLANE_ANGLE_UNIT() AND SI_UNIT('', .RADIAN.) )")
    solid_angle_unit = writer.add("( NAMED_UNIT(*) AND SOLID_ANGLE_UNIT() AND SI_UNIT('', .STERADIAN.) )")
    uncertainty = writer.add("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06), #" + str(length_unit) + ", 'distance accuracy value', 'confusion accuracy')")
    geom_context = writer.add(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) AND GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#" + str(uncertainty) + "))"
        " AND GLOBAL_UNIT_ASSIGNED_CONTEXT((#" + str(length_unit) + ", #" + str(plane_angle_unit) + ", #" + str(solid_angle_unit) + "))"
        " AND REPRESENTATION_CONTEXT('', '') )"
    )

    brep_shape = writer.add(
        f"ADVANCED_BREP_SHAPE_REPRESENTATION('', (#{brep_id}), #{geom_context})"
    )
    writer.add(f"SHAPE_DEFINITION_REPRESENTATION(#{shape}, #{brep_shape})")

    close_stream = False
    if hasattr(path_or_file, 'write'):
        stream = path_or_file
    else:
        stream = open(path_or_file, 'w', encoding='utf-8')
        close_stream = True

    try:
        _write_header(stream, name, schema)
        stream.write("DATA;\n")
        writer.write(stream)
        stream.write("ENDSEC;\nEND-ISO-10303-21;\n")
    finally:
        if close_stream:
            stream.close()


def _write_header(stream: TextIO, name: str, schema: str) -> None:
    stream.write("ISO-10303-21;\n")
    stream.write("HEADER;\n")
    stream.write("FILE_DESCRIPTION(('yapCAD export'),'2;1');\n")
    stream.write("FILE_NAME('" + name + "','2024-01-01T00:00:00',('yapCAD'),'', 'yapCAD','', '');\n")
    stream.write("FILE_SCHEMA(('" + schema + "'));\n")
    stream.write("ENDSEC;\n")


def _normalize(vec: Vec3) -> Vec3:
    length = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5
    if length == 0.0:
        return (0.0, 0.0, 1.0)
    return (vec[0] / length, vec[1] / length, vec[2] / length)


def _distance(a: Vec3, b: Vec3) -> float:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def _perpendicular(normal: Vec3) -> Vec3:
    x, y, z = normal
    if abs(x) < abs(z):
        return _normalize((-z, 0.0, x))
    return _normalize((0.0, z, -y))


__all__ = ['write_step']
