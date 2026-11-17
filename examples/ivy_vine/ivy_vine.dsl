module ivy_vine

# =============================================================================
#  IVY'S VINE  --  a print-in-place articulated ivy vine
#
#  A name leaf with IVY cut through it, and a vine that grows out of its
#  stalk.  The vine is a recursion: every link carries a leaf, a hinge to the
#  next link, and -- while the recursion has depth left -- a bud from which a
#  smaller copy of the whole vine grows.  The leaves' veins branch by the
#  same kind of rule.
#
#  Everything is a signed distance field (build with --representation sdf):
#    * leaf outlines are smooth unions of pointed lobes, so the sinuses
#      between lobes are the blend itself;
#    * a hinge's running clearance is the *offset* of the mating part's
#      field, so the gap is the same everywhere, bevels included;
#    * the joints are revolute mates on declared datums, which is what lets
#      yapcad.rom sweep them for clearance without meshing anything.
#
#  Prints flat on the bed, in one piece, with no supports: every overhang is
#  a 45 degree cone or bevel.
#
#  Units are millimetres.  Angles are degrees unless a name says otherwise.
# =============================================================================


# ----------------------------------------------------------------------------
#  Shared dimensions
# ----------------------------------------------------------------------------

command gap() -> float:
    # Running clearance between a hinge's two parts, measured normal to the
    # surfaces.  0.35 suits a 0.4 mm nozzle at 0.2 mm layers.
    emit 0.35

command ring_h() -> float:
    # Height of a hinge ring, and of the stems.
    emit 3.6

command knuckle_h() -> float:
    # Height of a hinge pin and the arch that holds its ring down.
    emit 6.4

command flare() -> float:
    # How far a pin's foot spreads: what stops a ring sliding off downwards.
    emit 1.0

command leaf_t() -> float:
    emit 2.8

# The per-level tables.  A link's hinge datums are literals in its @meta, so
# these must agree with the datum table on VINE_LINK below.

command pitch(level: int) -> float:
    if level <= 0:
        emit 24.0
    elif level == 1:
        emit 20.0
    emit 16.5

command bud_x(level: int) -> float:
    if level <= 0:
        emit 11.0
    emit 9.5

command bud_y(level: int) -> float:
    if level <= 0:
        emit 11.0
    emit 9.5

command pin_r(level: int) -> float:
    if level <= 0:
        emit 3.0
    elif level == 1:
        emit 2.6
    emit 2.3

command wall(level: int) -> float:
    if level <= 0:
        emit 2.8
    elif level == 1:
        emit 2.4
    emit 2.2

command stem_w(level: int) -> float:
    if level <= 0:
        emit 4.6
    elif level == 1:
        emit 4.0
    emit 3.6

command leaf_size(level: int) -> float:
    if level <= 0:
        emit 16.0
    elif level == 1:
        emit 13.5
    emit 11.5

command ring_ro(level: int) -> float:
    emit pin_r(level) + gap() + wall(level)


# ----------------------------------------------------------------------------
#  Small helpers
# ----------------------------------------------------------------------------

# A flat bar from (ax, ay) to (bx, by), standing on z = 0, edges bevelled.
# Its ends are square, so run it a little way into whatever it joins.
command bar(ax: float, ay: float, bx: float, by: float,
            width: float, height: float) -> solid:
    let dx: float = bx - ax
    let dy: float = by - ay
    let length: float = sqrt(dx * dx + dy * dy)
    let blank: solid = chamfer(box(length, width, height), 0.6)
    let turned: solid = rotate(blank, 0.0, 0.0, degrees(atan2(dy, dx)))
    emit translate(turned, (ax + bx) / 2.0, (ay + by) / 2.0, height / 2.0)

# Everything above z = 0, out to +/- half_size.
command above_bed(half_size: float) -> solid:
    emit translate(box(2.0 * half_size, 2.0 * half_size, 40.0), 0.0, 0.0, 20.0)


# ----------------------------------------------------------------------------
#  The hinge
#
#  Parent side: a pin whose foot flares at 45 degrees, an arm reaching it
#  from behind, and under the arm an arch -- the offset of the child's ring.
#  Child side: a ring whose bore is the offset of the pin and whose top is a
#  45 degree ridge, so the arch above it is a ridge too and needs no bridge.
#  The flare holds the ring up, the arch holds it down, nothing touches.
# ----------------------------------------------------------------------------

command hinge_pin(level: int) -> solid:
    let r: float = pin_r(level)
    emit union(cylinder(r, knuckle_h()), cone(r + flare(), r, flare()))

# The ring without its bore: a cylinder whose top is a 45 degree ridge, an
# outer bevel meeting an inner countersink.
command ring_blank(level: int) -> solid:
    let ro: float = ring_ro(level)
    let c: float = wall(level) / 2.0
    let h: float = ring_h()
    let rb: float = pin_r(level) + gap()
    # The countersink is a whole cone, apex below the bed, so that the arch
    # offset from it slopes at 45 degrees all the way in to the pin.
    let sink: solid = translate(cone(0.0, rb + c + 1.0, rb + c + 1.0),
                                0.0, 0.0, h - c - rb)
    emit difference(intersection(cylinder(ro, h), cone(ro + h - c, ro - c, h)),
                    sink)

# The bore of a ring: the pin grown by the running clearance.  A link cuts
# this last, so nothing it is built from can reach into its own hinge.
command ring_bore(level: int) -> solid:
    emit offset_solid(hinge_pin(level), gap())

# The space a hinge's ring needs, whichever way it is turned: a ring is a
# solid of revolution about the joint, so its sweep is itself.  Anything
# else on the parent is trimmed back to this.
command ring_keepout(level: int) -> solid:
    emit translate(cylinder(ring_ro(level) + gap() + 0.7, 30.0),
                   0.0, 0.0, -10.0)

# How far behind the joint centre a knuckle's arm reaches.
command knuckle_reach(level: int) -> float:
    emit ring_ro(level) + gap() + 3.0

# The parent's end of a hinge, centred on the origin, its arm reaching back
# along the direction back_deg.
command hinge_knuckle(level: int, back_deg: float) -> solid:
    let reach: float = knuckle_reach(level)
    let width: float = 2.0 * pin_r(level)
    let h: float = knuckle_h()
    let arm: solid = translate(box(reach, width, h), reach / 2.0, 0.0, h / 2.0)
    let arch: solid = offset_solid(ring_blank(level), gap())
    let body: solid = union(difference(arm, arch), hinge_pin(level))
    emit rotate(body, 0.0, 0.0, back_deg)


# ----------------------------------------------------------------------------
#  The ivy leaf
#
#  Five pointed lobes fanned from the stalk, each the lens where two circles
#  overlap, smooth-unioned so the sinuses between them come from the blend.
#  The leaf lies on z = 0 with its stalk end at the origin and its tip
#  along +x.
# ----------------------------------------------------------------------------

# One lobe: a lens of the given length and width, tall enough to be cut to
# thickness later, its near point a little behind the origin.
command lobe(length: float, width: float, angle: float) -> solid:
    let r: float = (length * length / (2.0 * width) + width / 2.0) / 2.0
    let d: float = r - width / 2.0
    let cx: float = length * 0.38
    let a: solid = translate(cylinder(r, 20.0), cx, d, -10.0)
    let b: solid = translate(cylinder(r, 20.0), cx, 0.0 - d, -10.0)
    emit rotate(intersection(a, b), 0.0, 0.0, angle)

# The veins of one lobe, as grooves to cut: a vein, and from two points
# along it a smaller pair of veins, and so on down.
command veins(x: float, y: float, angle: float, length: float,
              z: float, r0: float, depth: int) -> list<solid>:
    let ex: float = x + length * cos(radians(angle))
    let ey: float = y + length * sin(radians(angle))
    let r: float = r0 * (1.0 + 0.22 * depth)
    let own: list<solid> = [capsule(point(x, y, z), point(ex, ey, z), r)]
    if depth <= 0 or length < 3.0:
        emit own
    let mx: float = x + 0.42 * length * cos(radians(angle))
    let my: float = y + 0.42 * length * sin(radians(angle))
    let nx: float = x + 0.70 * length * cos(radians(angle))
    let ny: float = y + 0.70 * length * sin(radians(angle))
    let near: list<solid> = concat(
        veins(mx, my, angle + 42.0, 0.34 * length, z, r0, depth - 1),
        veins(mx, my, angle - 42.0, 0.34 * length, z, r0, depth - 1))
    let far: list<solid> = concat(
        veins(nx, ny, angle + 38.0, 0.22 * length, z, r0, depth - 1),
        veins(nx, ny, angle - 38.0, 0.22 * length, z, r0, depth - 1))
    emit concat(own, concat(near, far))

# The blade: lobes blended, cut to thickness with a rounded top edge.
# `fullness` 1.0 gives the sharp five-lobed leaf of a climbing shoot; larger
# values fatten the lobes and fill the sinuses, toward the broad, barely
# lobed leaf ivy grows once it is mature.
command leaf_blade(size: float, thickness: float, fullness: float) -> solid:
    let blend: float = (0.06 * size + 0.5) * fullness * fullness
    let middle: solid = lobe(size, 0.40 * size * fullness, 0.0)
    let sides: solid = smooth_union(
        lobe(0.78 * size, 0.33 * size * fullness, 58.0),
        lobe(0.78 * size, 0.33 * size * fullness, -58.0), blend)
    let basal: solid = smooth_union(
        lobe(0.50 * size, 0.30 * size * fullness, 122.0),
        lobe(0.50 * size, 0.30 * size * fullness, -122.0), blend)
    let outline: solid = smooth_union(smooth_union(middle, sides, blend),
                                      basal, blend)
    let lid: solid = translate(box(4.0 * size, 4.0 * size, 20.0),
                               0.0, 0.0, thickness - 10.0)
    emit intersection(smooth_intersection(outline, lid, 0.9),
                      above_bed(2.0 * size))

# A whole leaf: blade, veins cut into its top face, and a stalk of the given
# length running back along -x.
command leaf(size: float, thickness: float, stalk: float,
             vein_depth: int) -> solid:
    let z: float = thickness + 0.1
    let grooves: list<solid> = concat(
        veins(0.04 * size, 0.0, 0.0, 0.80 * size, z, 0.40, vein_depth),
        concat(
            concat(veins(0.0, 0.0, 58.0, 0.58 * size, z, 0.40, vein_depth),
                   veins(0.0, 0.0, -58.0, 0.58 * size, z, 0.40, vein_depth)),
            concat(veins(0.0, 0.0, 122.0, 0.33 * size, z, 0.40, vein_depth - 1),
                   veins(0.0, 0.0, -122.0, 0.33 * size, z, 0.40, vein_depth - 1))))
    let blade: solid = difference_all(leaf_blade(size, thickness, 1.0), grooves)
    emit union(blade, bar(0.0 - stalk, 0.0, 0.10 * size, 0.0, 2.6, thickness))


# ----------------------------------------------------------------------------
#  A vine link
#
#  Root ring at the origin, stem along +x.  It may end in a pin for the next
#  link of the same stem ("tip"), and may carry a pin for a side shoot one
#  level smaller ("bud") on the side given by `side` (+1 left, -1 right).
#  Its leaf grows on the other side, or straight ahead if the stem ends here.
#
#  The datum table lists every hinge axis a link of any level can have; a
#  link is built with the ones its level uses.
# ----------------------------------------------------------------------------

@meta(assembly.datums=[
    {"id": "root", "kind": "axis",
     "origin_mm": [0.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "tip0", "kind": "axis",
     "origin_mm": [24.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "tip1", "kind": "axis",
     "origin_mm": [20.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "tip2", "kind": "axis",
     "origin_mm": [16.5, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "budL0", "kind": "axis",
     "origin_mm": [11.0, 11.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "budR0", "kind": "axis",
     "origin_mm": [11.0, -11.0, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "budL1", "kind": "axis",
     "origin_mm": [9.5, 9.5, 0.0], "direction": [0.0, 0.0, 1.0]},
    {"id": "budR1", "kind": "axis",
     "origin_mm": [9.5, -9.5, 0.0], "direction": [0.0, 0.0, 1.0]}
])
@meta(assembly.joint_kind="axial", material="PLA")
command VINE_LINK(
    level: int @ui(label="Level", min=0, max=2) = 0,
    has_tip: bool @ui(label="Continues the stem") = true,
    has_bud: bool @ui(label="Carries a side shoot") = true,
    side: float @ui(label="Bud side (+1 left, -1 right)") = 1.0
) -> solid:
    require level >= 0 and level <= 2
    let p: float = pitch(level)
    let ws: float = stem_w(level)
    let size: float = leaf_size(level)
    let parts: list<solid> = [ring_blank(level)]

    if has_tip:
        # Stem to the tip knuckle, and a leaf on the side away from the bud,
        # trimmed clear of the next link's ring.
        let blade: solid = translate(
            rotate(leaf(size, leaf_t(), 3.0, 1), 0.0, 0.0, 0.0 - side * 98.0),
            0.44 * p, 0.0 - side * (ws / 2.0 + 2.2), 0.0)
        parts = concat(parts, [
            bar(ring_ro(level) - 0.8, 0.0,
                p - knuckle_reach(level) + 1.5, 0.0, ws, ring_h()),
            translate(hinge_knuckle(level, 180.0), p, 0.0, 0.0),
            difference(blade, translate(ring_keepout(level), p, 0.0, 0.0))])
    else:
        # The stem ends in a leaf, a little larger, pointing ahead.
        parts = concat(parts, [
            translate(leaf(1.15 * size, leaf_t(), 3.0, 1),
                      ring_ro(level) + 2.4, 0.0, 0.0)])

    if has_bud:
        # The bud's arm reaches back to the stem; its shoot leaves at about
        # 60 degrees on the bud's side.
        let bx: float = bud_x(level)
        let by: float = side * bud_y(level)
        let back: float = 0.0 - side * 120.0
        let reach: float = knuckle_reach(level + 1) - 1.5
        parts = concat(parts, [
            translate(hinge_knuckle(level + 1, back), bx, by, 0.0),
            bar(bx + reach * cos(radians(back)), by + reach * sin(radians(back)),
                ring_ro(level) + pin_r(level + 1) + 1.0, 0.0,
                2.0 * pin_r(level + 1), ring_h())])

    emit difference(union_all(parts), ring_bore(level))


# ----------------------------------------------------------------------------
#  The name leaf
#
#  A large leaf with I, V and Y cut through it, each letter a blend of
#  capsules.  Held by the leaf with the vine hanging down, the name reads
#  the right way up.  Its stalk ends in the pin the vine grows from.
# ----------------------------------------------------------------------------

# One stroke of a letter, as a cutter with near-vertical walls.
command stroke(ax: float, ay: float, bx: float, by: float, z: float) -> solid:
    let pen: solid = capsule(point(ax, ay, 0.0), point(bx, by, 0.0), 1.4)
    emit translate(scale(pen, 1.0, 1.0, 6.0), 0.0, 0.0, z)

# IVY, in letter coordinates: u to the right, v upwards, origin at the
# middle of the baseline.  (u, v) maps to leaf (x, y) = (x0 - v, u).
command name_cutters(x0: float, height: float, z: float) -> list<solid>:
    let top: float = x0 - height
    let mid: float = x0 - 0.46 * height
    emit [
        # I
        stroke(x0, -12.5, top, -12.5, z),
        # V
        smooth_union(stroke(top, -8.2, x0, -4.0, z),
                     stroke(x0, -4.0, top, 0.2, z), 0.5),
        # Y
        smooth_union(
            smooth_union(stroke(top, 5.0, mid, 8.7, z),
                         stroke(top, 12.4, mid, 8.7, z), 0.5),
            stroke(mid, 8.7, x0, 8.7, z), 0.5)
    ]

@meta(assembly.datums=[
    {"id": "stem", "kind": "axis",
     "origin_mm": [23.0, 0.0, 0.0], "direction": [0.0, 0.0, 1.0]}
])
@meta(assembly.joint_kind="axial", material="PLA", tags=["name_leaf"])
command NAME_LEAF(
    size: float @ui(widget="slider", label="Leaf size", min=34.0, max=44.0) = 38.0
) -> solid:
    let t: float = 3.2
    # Leaf tip along -x, stalk toward the vine on +x.
    let blade: solid = rotate(leaf_blade(size, t, 1.55), 0.0, 0.0, 180.0)
    let grooves: list<solid> = concat(
        veins(0.0 - 0.52 * size, 0.0, 180.0, 0.30 * size, t, 0.62, 1),
        concat(veins(0.08 * size, 0.24 * size, 96.0, 0.18 * size, t, 0.62, 1),
               veins(0.08 * size, 0.0 - 0.24 * size, -96.0, 0.18 * size, t, 0.62, 1)))
    let cut: solid = difference_all(
        blade, concat(concat(grooves,
                             name_cutters(0.0 - 0.07 * size, 0.27 * size, t)),
                      [translate(ring_keepout(0), 23.0, 0.0, 0.0)]))
    emit union_all([
        cut,
        bar(1.0, 0.0, 23.0 - knuckle_reach(0) + 1.5, 0.0,
            stem_w(0), ring_h()),
        translate(hinge_knuckle(0, 180.0), 23.0, 0.0, 0.0)])


# ----------------------------------------------------------------------------
#  Growing the vine
#
#  grow() adds one link to `parent` at `datum`, then calls itself twice: once
#  to continue the stem, and once -- from the bud -- to grow a whole smaller
#  vine.  That second call is what makes the vine self-similar.  Parts are
#  named by the path that reached them: "vsb" is the vine's first link (v),
#  the next along its stem (s), then that link's side shoot (b).
# ----------------------------------------------------------------------------

# The range a side shoot's joint can turn through, measured from its
# parent's stem direction, before its first link meets the parent.  These
# and the stem swing were set from yapcad.rom sweeps of the joints.
command bud_lo() -> float:
    emit 55.0

command bud_hi() -> float:
    emit 135.0

command shoot_length(level: int, remaining: int) -> int:
    # A side shoot is one link shorter than what is left of its stem.
    emit remaining - 1

command grow(a: assembly, parent: string, datum: string, name: string,
             level: int, remaining: int, max_level: int,
             side: float) -> assembly:
    if remaining <= 0:
        emit a
    let has_tip: bool = remaining > 1
    let has_bud: bool = level < max_level and shoot_length(level, remaining) >= 1

    add_part(a, VINE_LINK(level, has_tip, has_bud, side), name)
    add_named_mate(a, str("j_", name), "revolute", parent, datum, name, "root")

    # The stem carries on, its next bud on the other side ...
    grow(a, name, str("tip", level), str(name, "s"),
         level, remaining - 1, max_level, 0.0 - side)
    if has_bud:
        # ... and from this link's bud, a whole smaller vine.
        let bud: string = str("budL", level) if side > 0.0 else str("budR", level)
        grow(a, name, bud, str(name, "b"),
             level + 1, shoot_length(level, remaining), max_level, side)
    emit a

# Walk the same tree again once it is solved: put each joint at the angle it
# is printed at, and declare how far it can turn.  A stem joint swings
# either way about straight ahead.  A side shoot leaves on its bud's side
# and stays there, between bud_lo() and bud_hi() from its parent's stem.
command settle(a: assembly, name: string, level: int, remaining: int,
               max_level: int, side: float, is_bud: bool,
               wander_deg: float, swing_deg: float,
               bud_deg: float) -> assembly:
    if remaining <= 0:
        emit a
    let joint: string = str("j_", name)
    if is_bud:
        set_joint_position(a, joint, radians(side * bud_deg))
        let lo: float = bud_lo() if side > 0.0 else 0.0 - bud_hi()
        let hi: float = bud_hi() if side > 0.0 else 0.0 - bud_lo()
        set_mate_limits(a, joint, radians(lo), radians(hi))
    else:
        set_joint_position(a, joint, radians(side * wander_deg))
        set_mate_limits(a, joint, radians(0.0 - swing_deg), radians(swing_deg))

    settle(a, str(name, "s"), level, remaining - 1, max_level, 0.0 - side,
           false, wander_deg, swing_deg, bud_deg)
    if level < max_level and shoot_length(level, remaining) >= 1:
        settle(a, str(name, "b"), level + 1, shoot_length(level, remaining),
               max_level, side, true, wander_deg, swing_deg, bud_deg)
    emit a

@meta(material="PLA", tags=["print_in_place"])
command IVY_VINE(
    links: int @ui(label="Links in the main stem", min=1, max=6) = 4,
    max_level: int @ui(label="Levels of side shoots", min=0, max=2) = 2,
    swing_deg: float @ui(widget="slider", label="Stem joint swing either way",
                         min=10.0, max=60.0) = 48.0,
    wander_deg: float @ui(widget="slider", label="Stem wander as printed",
                          min=0.0, max=30.0) = 8.0,
    bud_deg: float @ui(widget="slider", label="Side shoot angle as printed",
                       min=55.0, max=120.0) = 74.0
) -> solid:
    require links >= 1 and links <= 6
    require max_level >= 0 and max_level <= 2
    require wander_deg <= swing_deg
    require bud_deg >= bud_lo() and bud_deg <= bud_hi()
    let a: assembly = assembly("ivy_vine")
    add_part(a, NAME_LEAF(), "ivy")
    grow(a, "ivy", "stem", "v", 0, links, max_level, 1.0)
    solve_assembly(a, "ivy")
    settle(a, "v", 0, links, max_level, 1.0, false,
           wander_deg, swing_deg, bud_deg)
    emit assembly_compound(a)
