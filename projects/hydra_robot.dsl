module hydra_robot

# Hydra Modular Reconfigurable Robot — Starfish Configuration
# All dimensions in millimeters.
#
# PRIMITIVE ORIGINS:
#   box(w,d,h)      → centered at origin (±w/2, ±d/2, ±h/2)
#   cylinder(r,h)   → base at z=0, top at z=h, centered on XY
#   sphere(r)        → centered at origin
#
# DESIGN:
#   Hexagonal body (6 faces) with 5 radial arms (asymmetric starfish).
#   Arms are cylindrical tubes. Face 6 left open for docking.
#   Body height: 80mm. Arm diameter: 56mm, length: 120mm.
#   Hexagon circumradius: 70mm (flat-to-flat ~121mm).
#
# LAYOUT:
#   5 arms at 72° intervals (360/5) around Z axis.
#   Arm angles: 0°, 72°, 144°, 216°, 288°
#   Each arm: connector flush on body face → arm cylinder → connector → effector

command connector_ring() -> solid:
    # 46mm OD ring, 5mm tall, base at z=0
    let outer: solid = cylinder(23.0, 5.0)
    let inner: solid = cylinder(10.0, 6.0)
    let ring: solid = difference(outer, inner)
    emit ring

command hexagonal_body() -> solid:
    # Hexagonal prism via regular_polygon extrusion
    # Circumradius 70mm, height 80mm
    let hex_profile: region2d = regular_polygon(6, 70.0)
    let body: solid = extrude(hex_profile, 80.0)
    # Body base at z=0, top at z=80. Center vertically.
    let centered: solid = translate(body, 0.0, 0.0, -40.0)
    emit centered

command arm_cylinder() -> solid:
    # Cylindrical arm segment, r=28mm, length=120mm
    # Base at z=0, extends to z=120
    let arm: solid = cylinder(28.0, 120.0)
    emit arm

command gripper_unit() -> solid:
    # Base plate + two fingers. Base at z=0.
    let base: solid = cylinder(25.0, 10.0)
    let fl: solid = box(10.0, 10.0, 70.0)
    let flp: solid = translate(fl, -15.0, 0.0, 10.0)
    let fr: solid = box(10.0, 10.0, 70.0)
    let frp: solid = translate(fr, 15.0, 0.0, 10.0)
    let g1: solid = union(base, flp)
    let g2: solid = union(g1, frp)
    emit g2

command wheel_unit() -> solid:
    # Hub + rim. Base at z=0.
    let hub: solid = cylinder(25.0, 8.0)
    let rim_o: solid = cylinder(35.0, 20.0)
    let rim_i: solid = cylinder(28.0, 22.0)
    let rim: solid = difference(rim_o, rim_i)
    let rim_up: solid = translate(rim, 0.0, 0.0, 8.0)
    let w: solid = union(hub, rim_up)
    emit w

# === PER-LAYER COMMANDS ===

command BODY() -> solid:
    let body: solid = hexagonal_body()
    emit body, name="body", layer="body", material="aluminum"

command CONNECTORS() -> solid:
    let cr: solid = connector_ring()

    # 5 connector pairs at 72° intervals
    # Body hex circumradius = 70mm. Connector at face.
    # Arm 0 at 0° (+X direction)
    let c0a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c0b: solid = translate(c0a, 70.0, 0.0, 0.0)

    # Arm 1 at 72°
    let c1a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c1t: solid = translate(c1a, 70.0, 0.0, 0.0)
    let c1b: solid = rotate(c1t, 0.0, 0.0, 72.0)

    # Arm 2 at 144°
    let c2a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c2t: solid = translate(c2a, 70.0, 0.0, 0.0)
    let c2b: solid = rotate(c2t, 0.0, 0.0, 144.0)

    # Arm 3 at 216°
    let c3a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c3t: solid = translate(c3a, 70.0, 0.0, 0.0)
    let c3b: solid = rotate(c3t, 0.0, 0.0, 216.0)

    # Arm 4 at 288°
    let c4a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let c4t: solid = translate(c4a, 70.0, 0.0, 0.0)
    let c4b: solid = rotate(c4t, 0.0, 0.0, 288.0)

    # Outer connectors at arm tips (70 + 5 + 120 = 195)
    let co0a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co0b: solid = translate(co0a, 195.0, 0.0, 0.0)

    let co1a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co1t: solid = translate(co1a, 195.0, 0.0, 0.0)
    let co1b: solid = rotate(co1t, 0.0, 0.0, 72.0)

    let co2a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co2t: solid = translate(co2a, 195.0, 0.0, 0.0)
    let co2b: solid = rotate(co2t, 0.0, 0.0, 144.0)

    let co3a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co3t: solid = translate(co3a, 195.0, 0.0, 0.0)
    let co3b: solid = rotate(co3t, 0.0, 0.0, 216.0)

    let co4a: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co4t: solid = translate(co4a, 195.0, 0.0, 0.0)
    let co4b: solid = rotate(co4t, 0.0, 0.0, 288.0)

    let u1: solid = union(c0b, c1b)
    let u2: solid = union(u1, c2b)
    let u3: solid = union(u2, c3b)
    let u4: solid = union(u3, c4b)
    let u5: solid = union(u4, co0b)
    let u6: solid = union(u5, co1b)
    let u7: solid = union(u6, co2b)
    let u8: solid = union(u7, co3b)
    let u9: solid = union(u8, co4b)

    emit u9, name="connectors", layer="connector", material="copper"

command ARMS() -> solid:
    let seg: solid = arm_cylinder()

    # Arm base at body face (70) + connector (5) = 75
    # Rotate cylinder so Z axis → radial outward, then translate out
    # Arm 0 at 0°
    let a0r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a0p: solid = translate(a0r, 75.0, 0.0, 0.0)

    # Arm 1 at 72°
    let a1r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a1t: solid = translate(a1r, 75.0, 0.0, 0.0)
    let a1p: solid = rotate(a1t, 0.0, 0.0, 72.0)

    # Arm 2 at 144°
    let a2r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a2t: solid = translate(a2r, 75.0, 0.0, 0.0)
    let a2p: solid = rotate(a2t, 0.0, 0.0, 144.0)

    # Arm 3 at 216°
    let a3r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a3t: solid = translate(a3r, 75.0, 0.0, 0.0)
    let a3p: solid = rotate(a3t, 0.0, 0.0, 216.0)

    # Arm 4 at 288°
    let a4r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a4t: solid = translate(a4r, 75.0, 0.0, 0.0)
    let a4p: solid = rotate(a4t, 0.0, 0.0, 288.0)

    let u1: solid = union(a0p, a1p)
    let u2: solid = union(u1, a2p)
    let u3: solid = union(u2, a3p)
    let u4: solid = union(u3, a4p)

    emit u4, name="arms", layer="arm", material="steel"

command EFFECTORS() -> solid:
    # Effector base at arm tip (75 + 120) + connector (5) = 200
    # Arms 0,1: grippers. Arm 2: wheel. Arms 3,4: open (no effector).

    # Gripper at 0°
    let g0: solid = gripper_unit()
    let g0r: solid = rotate(g0, 0.0, 90.0, 0.0)
    let g0t: solid = translate(g0r, 200.0, 0.0, 0.0)

    # Gripper at 72°
    let g1: solid = gripper_unit()
    let g1r: solid = rotate(g1, 0.0, 90.0, 0.0)
    let g1t: solid = translate(g1r, 200.0, 0.0, 0.0)
    let g1p: solid = rotate(g1t, 0.0, 0.0, 72.0)

    # Wheel at 144°
    let w2: solid = wheel_unit()
    let w2r: solid = rotate(w2, 0.0, 90.0, 0.0)
    let w2t: solid = translate(w2r, 200.0, 0.0, 0.0)
    let w2p: solid = rotate(w2t, 0.0, 0.0, 144.0)

    let u1: solid = union(g0t, g1p)
    let u2: solid = union(u1, w2p)

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
