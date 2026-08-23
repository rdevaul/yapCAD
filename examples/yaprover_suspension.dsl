module yaprover_suspension

# Tests-first, print-envelope-sized proxy for the YapRover suspension.
#
# Coordinate convention:
#   +X forward, +Y left, +Z up.  Every suspension and wheel joint uses a
#   local +Y axis.  The simple boxes and cylinders intentionally stay cheap
#   to evaluate while retaining analytic BREP in a pythonocc environment.

@meta(assembly.datums=[
    {"id": "left_rocker", "kind": "axis",
     "origin_mm": [0.0, 155.0, 0.0], "direction": [0.0, 1.0, 0.0]},
    {"id": "right_rocker", "kind": "axis",
     "origin_mm": [0.0, -155.0, 0.0], "direction": [0.0, 1.0, 0.0]}
])
command CHASSIS() -> solid:
    emit translate(box(205.0, 205.0, 50.0), 0.0, 0.0, -25.0)

@meta(assembly.datums=[
    {"id": "chassis_pivot", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 1.0, 0.0]},
    {"id": "front_axle", "kind": "axis",
     "origin_mm": [190.0, 0.0, -75.0], "direction": [0.0, 1.0, 0.0]},
    {"id": "bogie_pivot", "kind": "axis",
     "origin_mm": [-95.0, 0.0, -35.0], "direction": [0.0, 1.0, 0.0]}
])
command ROCKER() -> solid:
    # Compact kinematic proxy. The production rocker will be split at the
    # hub into front/rear printed beams; the remote datums retain exact span.
    emit translate(box(204.0, 16.0, 24.0), 0.0, 0.0, -50.0)

@meta(assembly.datums=[
    {"id": "rocker_pivot", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 1.0, 0.0]},
    {"id": "middle_axle", "kind": "axis",
     "origin_mm": [95.0, 0.0, -40.0], "direction": [0.0, 1.0, 0.0]},
    {"id": "rear_axle", "kind": "axis",
     "origin_mm": [-95.0, 0.0, -40.0], "direction": [0.0, 1.0, 0.0]}
])
command BOGIE() -> solid:
    emit translate(box(205.0, 16.0, 24.0), 0.0, 0.0, -52.0)

@meta(assembly.datums=[
    {"id": "axle", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 1.0, 0.0]}
])
command WHEEL() -> solid:
    # cylinder() is +Z; center it and turn its analytic axis onto +Y.
    let centered: solid = translate(cylinder(65.0, 36.0), 0.0, 0.0, -18.0)
    emit rotate(centered, -90.0, 0.0, 0.0)


command BUILD_SUSPENSION(rocker_angle: float = 0.0) -> solid:
    let rover: assembly = assembly("yaprover_suspension")
    add_part(rover, CHASSIS(), "chassis")
    add_part(rover, ROCKER(), "left_rocker")
    add_part(rover, ROCKER(), "right_rocker")
    add_part(rover, BOGIE(), "left_bogie")
    add_part(rover, BOGIE(), "right_bogie")
    add_part(rover, WHEEL(), "left_front_wheel")
    add_part(rover, WHEEL(), "left_middle_wheel")
    add_part(rover, WHEEL(), "left_rear_wheel")
    add_part(rover, WHEEL(), "right_front_wheel")
    add_part(rover, WHEEL(), "right_middle_wheel")
    add_part(rover, WHEEL(), "right_rear_wheel")

    add_named_mate(rover, "left_rocker_pivot", "revolute",
                   "chassis", "left_rocker", "left_rocker", "chassis_pivot")
    add_named_mate(rover, "right_rocker_pivot", "revolute",
                   "chassis", "right_rocker", "right_rocker", "chassis_pivot")
    add_named_mate(rover, "left_bogie_pivot", "revolute",
                   "left_rocker", "bogie_pivot", "left_bogie", "rocker_pivot")
    add_named_mate(rover, "right_bogie_pivot", "revolute",
                   "right_rocker", "bogie_pivot", "right_bogie", "rocker_pivot")

    add_named_mate(rover, "left_front_axle", "revolute",
                   "left_rocker", "front_axle", "left_front_wheel", "axle")
    add_named_mate(rover, "left_middle_axle", "revolute",
                   "left_bogie", "middle_axle", "left_middle_wheel", "axle")
    add_named_mate(rover, "left_rear_axle", "revolute",
                   "left_bogie", "rear_axle", "left_rear_wheel", "axle")
    add_named_mate(rover, "right_front_axle", "revolute",
                   "right_rocker", "front_axle", "right_front_wheel", "axle")
    add_named_mate(rover, "right_middle_axle", "revolute",
                   "right_bogie", "middle_axle", "right_middle_wheel", "axle")
    add_named_mate(rover, "right_rear_axle", "revolute",
                   "right_bogie", "rear_axle", "right_rear_wheel", "axle")

    set_mate_limits(rover, "left_rocker_pivot", radians(-18.0), radians(18.0))
    set_mate_limits(rover, "right_rocker_pivot", radians(-18.0), radians(18.0))
    set_mate_limits(rover, "left_bogie_pivot", radians(-35.0), radians(38.0))
    set_mate_limits(rover, "right_bogie_pivot", radians(-35.0), radians(38.0))

    # Passive differential: right rocker angle = -left rocker angle.
    add_joint_coupling(rover, "rocker_differential", "right_rocker_pivot",
                       ["left_rocker_pivot"], [-1.0], 0.0)
    solve_assembly(rover, "chassis")
    set_joint_position(rover, "left_rocker_pivot", rocker_angle)
    emit assembly_compound(rover)
