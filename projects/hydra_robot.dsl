module hydra_robot

# Hydra Modular Reconfigurable Robot — Assembly Visualization
# All dimensions in millimeters.

# === SHARED PRIMITIVES ===

command connector_ring() -> solid:
    let outer: solid = cylinder(23.0, 5.0)
    let inner: solid = cylinder(10.0, 6.0)
    let ring: solid = difference(outer, inner)
    emit ring

command arm_segment() -> solid:
    let beam: solid = box(56.0, 56.0, 110.0)
    let bc: solid = translate(beam, -28.0, -28.0, -55.0)
    let js: solid = sphere(24.0)
    let jp: solid = translate(js, 0.0, 0.0, 55.0)
    let arm: solid = union(bc, jp)
    emit arm

command gripper_unit() -> solid:
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

# === PER-LAYER COMMANDS (for multi-layer packaging) ===

command BODY() -> solid:
    let body: solid = box(120.0, 120.0, 80.0)
    let bc: solid = translate(body, -60.0, -60.0, -40.0)
    emit bc, name="body", layer="body", material="aluminum"

command CONNECTORS() -> solid:
    let cr: solid = connector_ring()

    # Body face connectors
    let c1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c1p: solid = translate(c1, 60.0, 0.0, 0.0)
    let c2: solid = rotate(cr, 0.0, -90.0, 0.0)
    let c2p: solid = translate(c2, -60.0, 0.0, 0.0)
    let c3: solid = rotate(cr, -90.0, 0.0, 0.0)
    let c3p: solid = translate(c3, 0.0, 60.0, 0.0)
    let c4: solid = rotate(cr, 90.0, 0.0, 0.0)
    let c4p: solid = translate(c4, 0.0, -60.0, 0.0)

    # Arm inner connectors (body-side)
    let ai1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ai1p: solid = translate(ai1, 70.0, 0.0, 0.0)
    let ai2: solid = rotate(cr, 0.0, -90.0, 0.0)
    let ai2p: solid = translate(ai2, -70.0, 0.0, 0.0)
    let ai3: solid = rotate(cr, -90.0, 0.0, 0.0)
    let ai3p: solid = translate(ai3, 0.0, 70.0, 0.0)
    let ai4: solid = rotate(cr, 90.0, 0.0, 0.0)
    let ai4p: solid = translate(ai4, 0.0, -70.0, 0.0)

    # Arm outer connectors (effector-side)
    let ao1: solid = rotate(cr, 0.0, -90.0, 0.0)
    let ao1p: solid = translate(ao1, 180.0, 0.0, 0.0)
    let ao2: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ao2p: solid = translate(ao2, -180.0, 0.0, 0.0)
    let ao3: solid = rotate(cr, 90.0, 0.0, 0.0)
    let ao3p: solid = translate(ao3, 0.0, 180.0, 0.0)
    let ao4: solid = rotate(cr, -90.0, 0.0, 0.0)
    let ao4p: solid = translate(ao4, 0.0, -186.0, 0.0)

    # Union all connectors
    let u1: solid = union(c1p, c2p)
    let u2: solid = union(u1, c3p)
    let u3: solid = union(u2, c4p)
    let u4: solid = union(u3, ai1p)
    let u5: solid = union(u4, ai2p)
    let u6: solid = union(u5, ai3p)
    let u7: solid = union(u6, ai4p)
    let u8: solid = union(u7, ao1p)
    let u9: solid = union(u8, ao2p)
    let u10: solid = union(u9, ao3p)
    let u11: solid = union(u10, ao4p)

    emit u11, name="connectors", layer="connector", material="copper"

command ARMS() -> solid:
    let seg: solid = arm_segment()

    let a1r: solid = rotate(seg, 0.0, -90.0, 0.0)
    let a1p: solid = translate(a1r, 125.0, 0.0, 0.0)

    let a2r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a2p: solid = translate(a2r, -125.0, 0.0, 0.0)

    let a3r: solid = rotate(seg, 90.0, 0.0, 0.0)
    let a3p: solid = translate(a3r, 0.0, 125.0, 0.0)

    let a4r: solid = rotate(seg, -90.0, 0.0, 0.0)
    let a4p: solid = translate(a4r, 0.0, -125.0, 0.0)

    let u1: solid = union(a1p, a2p)
    let u2: solid = union(u1, a3p)
    let u3: solid = union(u2, a4p)

    emit u3, name="arms", layer="arm", material="steel"

command EFFECTORS() -> solid:
    # Gripper +X
    let g1: solid = gripper_unit()
    let g1r: solid = rotate(g1, 0.0, -90.0, 0.0)
    let g1p: solid = translate(g1r, 190.0, 0.0, 0.0)

    # Gripper -X
    let g2: solid = gripper_unit()
    let g2r: solid = rotate(g2, 0.0, 90.0, 0.0)
    let g2p: solid = translate(g2r, -190.0, 0.0, 0.0)

    # Wheel +Y
    let w1: solid = wheel_unit()
    let w1r: solid = rotate(w1, 90.0, 0.0, 0.0)
    let w1p: solid = translate(w1r, 0.0, 190.0, 0.0)

    let u1: solid = union(g1p, g2p)
    let u2: solid = union(u1, w1p)

    emit u2, name="effectors", layer="effector", material="plastic"

# === FULL ASSEMBLY (all layers unioned for single-solid export) ===

command ASSEMBLY() -> solid:
    let body: solid = BODY()
    let conn: solid = CONNECTORS()
    let arms: solid = ARMS()
    let eff: solid = EFFECTORS()

    let u1: solid = union(body, conn)
    let u2: solid = union(u1, arms)
    let u3: solid = union(u2, eff)

    emit u3, name="hydra_assembly", layer="assembly"
