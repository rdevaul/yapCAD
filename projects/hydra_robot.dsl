module hydra_robot

# Hydra Modular Reconfigurable Robot — Starfish Configuration
# All dimensions in millimeters.
#
# Layout (X axis example):
#   Body face at x=60
#   Connector ring centered at x=60 (2.5mm into body, 2.5mm out)
#   Arm inner face at x=65 (5mm gap for connector pair)  
#   Arm outer face at x=175 (65 + 110)
#   Connector ring centered at x=175
#   Effector base at x=180 (175 + 5mm connector)

command connector_ring() -> solid:
    # Magnetic alignment ring + power coil housing
    # 5mm tall, centered at z=0 (extends -2.5 to +2.5)
    let outer: solid = cylinder(23.0, 5.0)
    let inner: solid = cylinder(10.0, 6.0)
    let ring: solid = difference(outer, inner)
    let centered: solid = translate(ring, 0.0, 0.0, -2.5)
    emit centered

command arm_segment() -> solid:
    # 56x56x110 beam, bottom at z=0, top at z=110
    let beam: solid = box(56.0, 56.0, 110.0)
    let bc: solid = translate(beam, -28.0, -28.0, 0.0)
    # Spherical joint cap at top
    let js: solid = sphere(24.0)
    let jp: solid = translate(js, 0.0, 0.0, 110.0)
    let arm: solid = union(bc, jp)
    emit arm

command gripper_unit() -> solid:
    # Base at z=0, fingers extend upward
    let base: solid = box(50.0, 50.0, 10.0)
    let bc: solid = translate(base, -25.0, -25.0, 0.0)
    let fl: solid = box(10.0, 10.0, 70.0)
    let flp: solid = translate(fl, -25.0, -5.0, 10.0)
    let fr: solid = box(10.0, 10.0, 70.0)
    let frp: solid = translate(fr, 15.0, -5.0, 10.0)
    let g1: solid = union(bc, flp)
    let g2: solid = union(g1, frp)
    emit g2

command wheel_unit() -> solid:
    let hub: solid = cylinder(30.0, 8.0)
    let rim_o: solid = cylinder(35.0, 18.0)
    let rim_i: solid = cylinder(28.0, 20.0)
    let rim: solid = difference(rim_o, rim_i)
    let rim_up: solid = translate(rim, 0.0, 0.0, 8.0)
    let w: solid = union(hub, rim_up)
    emit w

# === PER-LAYER COMMANDS ===

command BODY() -> solid:
    let body: solid = box(120.0, 120.0, 80.0)
    let bc: solid = translate(body, -60.0, -60.0, -40.0)
    emit bc, name="body", layer="body", material="aluminum"

command CONNECTORS() -> solid:
    let cr: solid = connector_ring()

    # Body-to-arm connectors at each face (x=±60, y=±60)
    let c1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c1p: solid = translate(c1, 62.5, 0.0, 0.0)

    let c2: solid = rotate(cr, 0.0, -90.0, 0.0)
    let c2p: solid = translate(c2, -62.5, 0.0, 0.0)

    let c3: solid = rotate(cr, -90.0, 0.0, 0.0)
    let c3p: solid = translate(c3, 0.0, 62.5, 0.0)

    let c4: solid = rotate(cr, 90.0, 0.0, 0.0)
    let c4p: solid = translate(c4, 0.0, -62.5, 0.0)

    # Arm-to-effector connectors at outer arm faces
    # Arm outer face at body_face(60) + gap(5) + arm(110) = 175
    let ao1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ao1p: solid = translate(ao1, 177.5, 0.0, 0.0)

    let ao2: solid = rotate(cr, 0.0, -90.0, 0.0)
    let ao2p: solid = translate(ao2, -177.5, 0.0, 0.0)

    let ao3: solid = rotate(cr, -90.0, 0.0, 0.0)
    let ao3p: solid = translate(ao3, 0.0, 177.5, 0.0)

    let ao4: solid = rotate(cr, 90.0, 0.0, 0.0)
    let ao4p: solid = translate(ao4, 0.0, -177.5, 0.0)

    let u1: solid = union(c1p, c2p)
    let u2: solid = union(u1, c3p)
    let u3: solid = union(u2, c4p)
    let u4: solid = union(u3, ao1p)
    let u5: solid = union(u4, ao2p)
    let u6: solid = union(u5, ao3p)
    let u7: solid = union(u6, ao4p)

    emit u7, name="connectors", layer="connector", material="copper"

command ARMS() -> solid:
    let seg: solid = arm_segment()

    # Arm inner face at body_face(60) + connector(5) = 65
    # Arm extends from 65 to 175 along each axis
    # +X arm: rotate so arm Z axis points along +X, base at x=65
    let a1r: solid = rotate(seg, 0.0, -90.0, 0.0)
    let a1p: solid = translate(a1r, 65.0, 0.0, 0.0)

    # -X arm: rotate so arm Z axis points along -X, base at x=-65
    let a2r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a2p: solid = translate(a2r, -65.0, 0.0, 0.0)

    # +Y arm
    let a3r: solid = rotate(seg, 90.0, 0.0, 0.0)
    let a3p: solid = translate(a3r, 0.0, 65.0, 0.0)

    # -Y arm
    let a4r: solid = rotate(seg, -90.0, 0.0, 0.0)
    let a4p: solid = translate(a4r, 0.0, -65.0, 0.0)

    let u1: solid = union(a1p, a2p)
    let u2: solid = union(u1, a3p)
    let u3: solid = union(u2, a4p)

    emit u3, name="arms", layer="arm", material="steel"

command EFFECTORS() -> solid:
    # Effector base at arm_outer(175) + connector(5) = 180
    # +X gripper
    let g1: solid = gripper_unit()
    let g1r: solid = rotate(g1, 0.0, -90.0, 0.0)
    let g1p: solid = translate(g1r, 180.0, 0.0, 0.0)

    # -X gripper
    let g2: solid = gripper_unit()
    let g2r: solid = rotate(g2, 0.0, 90.0, 0.0)
    let g2p: solid = translate(g2r, -180.0, 0.0, 0.0)

    # +Y wheel
    let w1: solid = wheel_unit()
    let w1r: solid = rotate(w1, 90.0, 0.0, 0.0)
    let w1p: solid = translate(w1r, 0.0, 180.0, 0.0)

    let u1: solid = union(g1p, g2p)
    let u2: solid = union(u1, w1p)

    emit u2, name="effectors", layer="effector", material="plastic"

# === FULL ASSEMBLY ===

command ASSEMBLY() -> solid:
    let body: solid = BODY()
    let conn: solid = CONNECTORS()
    let arms: solid = ARMS()
    let eff: solid = EFFECTORS()

    let u1: solid = union(body, conn)
    let u2: solid = union(u1, arms)
    let u3: solid = union(u2, eff)
    emit u3, name="hydra_assembly", layer="assembly"
