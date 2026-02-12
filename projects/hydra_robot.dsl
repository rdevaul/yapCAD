module hydra_robot

# Hydra Modular Reconfigurable Robot — Starfish Configuration
# All dimensions in millimeters.
#
# PRIMITIVE ORIGINS:
#   box(w,d,h)      → centered at origin
#   cylinder(r,h)   → base at z=0, top at z=h, centered on XY
#   sphere(r)       → centered at origin
#   extrude(p, h)   → base at z=0, top at z=h
#
# DESIGN:
#   Hexagonal body with 6 faces at 60° intervals.
#   5 arms on faces 0-4, face 5 (300°) left open.
#   Face normals at: 0°, 60°, 120°, 180°, 240°, 300°
#   Hex circumradius: 70mm. Body height: 80mm.
#   Arm: cylinder r=28mm, length=120mm.

command connector_ring() -> solid:
    let outer: solid = cylinder(23.0, 5.0)
    let inner: solid = cylinder(10.0, 6.0)
    let ring: solid = difference(outer, inner)
    emit ring

command hexagonal_body() -> solid:
    let hex_profile: region2d = regular_polygon(6, 70.0)
    let body: solid = extrude(hex_profile, 80.0)
    let centered: solid = translate(body, 0.0, 0.0, -40.0)
    emit centered

command arm_cylinder() -> solid:
    let arm: solid = cylinder(28.0, 120.0)
    emit arm

command gripper_unit() -> solid:
    let base: solid = cylinder(25.0, 10.0)
    let fl: solid = box(10.0, 10.0, 70.0)
    let flp: solid = translate(fl, -15.0, 0.0, 10.0)
    let fr: solid = box(10.0, 10.0, 70.0)
    let frp: solid = translate(fr, 15.0, 0.0, 10.0)
    let g1: solid = union(base, flp)
    let g2: solid = union(g1, frp)
    emit g2

command wheel_unit() -> solid:
    let hub: solid = cylinder(25.0, 8.0)
    let rim_o: solid = cylinder(35.0, 20.0)
    let rim_i: solid = cylinder(28.0, 22.0)
    let rim: solid = difference(rim_o, rim_i)
    let rim_up: solid = translate(rim, 0.0, 0.0, 8.0)
    let w: solid = union(hub, rim_up)
    emit w

# Helper: place object at hex face
# Pattern: rotate to face outward (+Y to +X via 90° Y rot),
#          translate to face radius, rotate around Z to face angle.
# Hex face distance (apothem) = circumradius * cos(30°) ≈ 70 * 0.866 = 60.6mm
# But we use circumradius for connector placement since faces are flat.
# Actually: hex apothem (center to face midpoint) = R * cos(30°) = 60.6mm

# === PER-LAYER COMMANDS ===

command BODY() -> solid:
    let body: solid = hexagonal_body()
    emit body, name="body", layer="body", material="aluminum"

command CONNECTORS() -> solid:
    let cr: solid = connector_ring()

    # Hex apothem = 70 * cos(30°) ≈ 60.6mm (distance from center to face midpoint)
    # Connector placed at apothem distance on each face normal
    # Face angles: 0°, 60°, 120°, 180°, 240°, 300°

    # Inner connectors (body side) — all 6 faces
    let ci0: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci0t: solid = translate(ci0, 60.6, 0.0, 0.0)

    let ci1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci1t: solid = translate(ci1, 60.6, 0.0, 0.0)
    let ci1p: solid = rotate(ci1t, 0.0, 0.0, 60.0)

    let ci2: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci2t: solid = translate(ci2, 60.6, 0.0, 0.0)
    let ci2p: solid = rotate(ci2t, 0.0, 0.0, 120.0)

    let ci3: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci3t: solid = translate(ci3, 60.6, 0.0, 0.0)
    let ci3p: solid = rotate(ci3t, 0.0, 0.0, 180.0)

    let ci4: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci4t: solid = translate(ci4, 60.6, 0.0, 0.0)
    let ci4p: solid = rotate(ci4t, 0.0, 0.0, 240.0)

    let ci5: solid = rotate(cr, 0.0, 90.0, 0.0)
    let ci5t: solid = translate(ci5, 60.6, 0.0, 0.0)
    let ci5p: solid = rotate(ci5t, 0.0, 0.0, 300.0)

    # Outer connectors (arm tip side) — faces 0-4 only
    # Arm tip at apothem(60.6) + connector(5) + arm(120) = 185.6
    let co0: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co0t: solid = translate(co0, 185.6, 0.0, 0.0)

    let co1: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co1t: solid = translate(co1, 185.6, 0.0, 0.0)
    let co1p: solid = rotate(co1t, 0.0, 0.0, 60.0)

    let co2: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co2t: solid = translate(co2, 185.6, 0.0, 0.0)
    let co2p: solid = rotate(co2t, 0.0, 0.0, 120.0)

    let co3: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co3t: solid = translate(co3, 185.6, 0.0, 0.0)
    let co3p: solid = rotate(co3t, 0.0, 0.0, 180.0)

    let co4: solid = rotate(cr, 0.0, 90.0, 0.0)
    let co4t: solid = translate(co4, 185.6, 0.0, 0.0)
    let co4p: solid = rotate(co4t, 0.0, 0.0, 240.0)

    # Union all
    let u1: solid = union(ci0t, ci1p)
    let u2: solid = union(u1, ci2p)
    let u3: solid = union(u2, ci3p)
    let u4: solid = union(u3, ci4p)
    let u5: solid = union(u4, ci5p)
    let u6: solid = union(u5, co0t)
    let u7: solid = union(u6, co1p)
    let u8: solid = union(u7, co2p)
    let u9: solid = union(u8, co3p)
    let u10: solid = union(u9, co4p)

    emit u10, name="connectors", layer="connector", material="copper"

command ARMS() -> solid:
    let seg: solid = arm_cylinder()

    # Arm base at apothem(60.6) + connector(5) = 65.6
    # Faces 0-4 get arms, face 5 (300°) is blank

    let a0r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a0p: solid = translate(a0r, 65.6, 0.0, 0.0)

    let a1r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a1t: solid = translate(a1r, 65.6, 0.0, 0.0)
    let a1p: solid = rotate(a1t, 0.0, 0.0, 60.0)

    let a2r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a2t: solid = translate(a2r, 65.6, 0.0, 0.0)
    let a2p: solid = rotate(a2t, 0.0, 0.0, 120.0)

    let a3r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a3t: solid = translate(a3r, 65.6, 0.0, 0.0)
    let a3p: solid = rotate(a3t, 0.0, 0.0, 180.0)

    let a4r: solid = rotate(seg, 0.0, 90.0, 0.0)
    let a4t: solid = translate(a4r, 65.6, 0.0, 0.0)
    let a4p: solid = rotate(a4t, 0.0, 0.0, 240.0)

    let u1: solid = union(a0p, a1p)
    let u2: solid = union(u1, a2p)
    let u3: solid = union(u2, a3p)
    let u4: solid = union(u3, a4p)

    emit u4, name="arms", layer="arm", material="steel"

command EFFECTORS() -> solid:
    # Effector base at apothem(60.6) + connector(5) + arm(120) + connector(5) = 190.6
    # Face 0: gripper, Face 1: gripper, Face 2: wheel
    # Faces 3,4: open (no effector)

    let g0: solid = gripper_unit()
    let g0r: solid = rotate(g0, 0.0, 90.0, 0.0)
    let g0p: solid = translate(g0r, 190.6, 0.0, 0.0)

    let g1: solid = gripper_unit()
    let g1r: solid = rotate(g1, 0.0, 90.0, 0.0)
    let g1t: solid = translate(g1r, 190.6, 0.0, 0.0)
    let g1p: solid = rotate(g1t, 0.0, 0.0, 60.0)

    let w2: solid = wheel_unit()
    let w2r: solid = rotate(w2, 0.0, 90.0, 0.0)
    let w2t: solid = translate(w2r, 190.6, 0.0, 0.0)
    let w2p: solid = rotate(w2t, 0.0, 0.0, 120.0)

    let u1: solid = union(g0p, g1p)
    let u2: solid = union(u1, w2p)

    emit u2, name="effectors", layer="effector", material="plastic"

command ASSEMBLY() -> solid:
    let body: solid = BODY()
    let conn: solid = CONNECTORS()
    let arms: solid = ARMS()
    let eff: solid = EFFECTORS()
    let u1: solid = union(body, conn)
    let u2: solid = union(u1, arms)
    let u3: solid = union(u2, eff)
    emit u3, name="hydra_assembly", layer="assembly"
