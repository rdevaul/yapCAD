module enclosure_shell_cutouts_v1

# =============================================================================
# Enclosure Shell Cutouts — generic test fixture (yapCAD metadata v1.1)
# =============================================================================
#
# Generic, project-neutral fixture exercising the @meta(operation.*) cutter
# workflow: three radial box cutters declared as separate commands, each
# decorated so the assembly resolver (BasicResolver) applies it as a boolean
# subtract against a cylindrical `enclosure_shell` target part.
#
# This mirrors the public "process-aware subtract cutter" capability without
# referencing any specific design. Parts and dimensions are invented.
#
# Geometry convention:
#   - Cutter local frame: box centered at origin (width=tangential,
#     depth=radial, height=axial). Cutter is shifted so its near (-X) face
#     sits at R_OUTER, rotated about Z to the clock angle, then translated
#     along Z to its axial position.
#   - Priorities 10/20/30 drive resolver application order.
# =============================================================================

command SHELL_PARAMS() -> float:
    # Enclosure outer skin radius (mm).
    R_OUTER : float = 157.0
    emit R_OUTER

# Helper: build a rectangular radial cutter positioned on the shell outer skin.
#
# Args:
#   tang_w   : tangential width, mm
#   axial_h  : axial height, mm
#   depth    : radial depth of cut, mm
#   clock_deg: clock angle, degrees
#   z_mm     : axial position of cutter center, mm
command RADIAL_BOX_CUTTER(
    tang_w    : float,
    axial_h   : float,
    depth     : float,
    clock_deg : float,
    z_mm      : float
) -> solid:
    R_OUTER : float = 157.0

    raw : solid = box(tang_w, depth, axial_h)
    shifted : solid = translate(raw, R_OUTER + depth / 2.0, 0.0, 0.0)
    rotated : solid = rotate(shifted, 0.0, 0.0, clock_deg)
    placed : solid = translate(rotated, 0.0, 0.0, z_mm)
    emit placed


# =============================================================================
# CUT_VENT_PANEL  (generic access panel A)
# =============================================================================

@meta(operation.kind="subtract", operation.target_filter=["enclosure_shell"],
      operation.priority=10, operation.through=true, operation.policy="strict",
      operation.feature_id="vent_panel", operation.feature_kind="vent")
command CUT_VENT_PANEL() -> solid:
    cutter : solid = RADIAL_BOX_CUTTER(
        tang_w    = 101.6,
        axial_h   = 127.0,
        depth     = 14.0,
        clock_deg = 247.5,
        z_mm      = 162.5
    )
    emit cutter


# =============================================================================
# CUT_ACCESS_DOOR  (generic access panel B)
# =============================================================================

@meta(operation.kind="subtract", operation.target_filter=["enclosure_shell"],
      operation.priority=20, operation.through=true, operation.policy="strict",
      operation.feature_id="access_door", operation.feature_kind="access_panel")
command CUT_ACCESS_DOOR() -> solid:
    cutter : solid = RADIAL_BOX_CUTTER(
        tang_w    = 165.1,
        axial_h   = 203.2,
        depth     = 30.0,
        clock_deg = 180.0,
        z_mm      = 162.5
    )
    emit cutter


# =============================================================================
# CUT_ACCESSORY_PORT  (generic access panel C)
# =============================================================================

@meta(operation.kind="subtract", operation.target_filter=["enclosure_shell"],
      operation.priority=30, operation.through=true, operation.policy="strict",
      operation.feature_id="accessory_port", operation.feature_kind="access_panel")
command CUT_ACCESSORY_PORT() -> solid:
    cutter : solid = RADIAL_BOX_CUTTER(
        tang_w    = 101.6,
        axial_h   = 127.0,
        depth     = 14.0,
        clock_deg = 90.0,
        z_mm      = 162.5
    )
    emit cutter
