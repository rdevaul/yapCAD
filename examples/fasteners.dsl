module fasteners

# Fastener generation examples using the yapCAD DSL
# Demonstrates metric (ISO) and unified (ASME) hex bolts and nuts
#
# Usage examples:
#   # Default M10x30 bolt
#   python -m yapcad.dsl run examples/fasteners.dsl MAKE_METRIC_BOLT
#
#   # Custom size bolt
#   python -m yapcad.dsl run examples/fasteners.dsl MAKE_METRIC_BOLT --param size=M8 --param length=25.0
#
#   # Bolt and nut pair
#   python -m yapcad.dsl run examples/fasteners.dsl MAKE_METRIC_BOLT_NUT_PAIR --param size=M10 --param length=30.0
#
#   # Export to STEP
#   python -m yapcad.dsl run examples/fasteners.dsl MAKE_METRIC_BOLT --output bolt.step
#
#   # Export to package
#   python -m yapcad.dsl run examples/fasteners.dsl MAKE_METRIC_BOLT_NUT_PAIR --package fasteners.ycpkg


# =============================================================================
# Metric Fasteners (ISO)
# =============================================================================

# Create a metric hex bolt (ISO 4014/4017)
# Available sizes: M3, M4, M5, M6, M8, M10, M12, M14, M16, M20, M24
command MAKE_METRIC_BOLT(
    size: string = "M10",
    length: float = 30.0
) -> solid:
    require length > 0.0
    bolt: solid = metric_hex_bolt(size, length)
    emit bolt


# Create a metric hex nut (ISO 4032)
command MAKE_METRIC_NUT(
    size: string = "M10"
) -> solid:
    nut: solid = metric_hex_nut(size)
    emit nut


# Create a matched bolt and nut pair, with nut positioned at thread end
command MAKE_METRIC_BOLT_NUT_PAIR(
    size: string = "M10",
    length: float = 30.0
) -> solid:
    require length > 0.0

    # Create bolt (head up, threads at Z=0)
    bolt: solid = metric_hex_bolt(size, length)

    # Create nut (base at Z=0)
    nut: solid = metric_hex_nut(size)

    # Position nut at the bottom of the thread
    # Nut sits near Z=0 where the thread tip is
    nut_positioned: solid = translate(nut, 0.0, 0.0, 2.0)

    # Combine into assembly using compound (not union, keeps separate bodies)
    assembly: solid = compound(bolt, nut_positioned)

    emit assembly


# =============================================================================
# Unified Fasteners (ASME - UNC)
# =============================================================================

# Create a unified (UNC) hex bolt (ASME B18.2.1)
# Available sizes: #4-40, #6-32, #8-32, #10-24, #12-24,
#                  1/4-20, 5/16-18, 3/8-16, 1/2-13, 5/8-11, 3/4-10, 1-8
# Note: length is in INCHES
command MAKE_UNIFIED_BOLT(
    size: string = "1/4-20",
    length: float = 1.0
) -> solid:
    require length > 0.0
    bolt: solid = unified_hex_bolt(size, length)
    emit bolt


# Create a unified hex nut (ASME B18.2.2)
command MAKE_UNIFIED_NUT(
    size: string = "1/4-20"
) -> solid:
    nut: solid = unified_hex_nut(size)
    emit nut


# Create a matched unified bolt and nut pair
command MAKE_UNIFIED_BOLT_NUT_PAIR(
    size: string = "1/4-20",
    length: float = 1.0
) -> solid:
    require length > 0.0

    # Create bolt (head up, threads at Z=0)
    bolt: solid = unified_hex_bolt(size, length)

    # Create nut (base at Z=0)
    nut: solid = unified_hex_nut(size)

    # Position nut at the bottom of the thread
    nut_positioned: solid = translate(nut, 0.0, 0.0, 2.0)

    # Combine into assembly
    assembly: solid = compound(bolt, nut_positioned)

    emit assembly


# =============================================================================
# Demo Commands
# =============================================================================

# Demo: Default M10x30 metric bolt
command DEMO_METRIC_BOLT() -> solid:
    bolt: solid = metric_hex_bolt("M10", 30.0)
    emit bolt


# Demo: Default 1/4-20 x 1" unified bolt
command DEMO_UNIFIED_BOLT() -> solid:
    bolt: solid = unified_hex_bolt("1/4-20", 1.0)
    emit bolt


# Demo: M8 bolt with M8 nut assembly
command DEMO_BOLT_NUT_ASSEMBLY() -> solid:
    # M8x25 bolt
    bolt: solid = metric_hex_bolt("M8", 25.0)

    # M8 nut positioned near thread end
    nut: solid = metric_hex_nut("M8")
    nut_positioned: solid = translate(nut, 0.0, 0.0, 3.0)

    # Create compound assembly
    assembly: solid = compound(bolt, nut_positioned)

    emit assembly


# Demo: Multiple fasteners arranged in a row
command DEMO_FASTENER_ARRAY() -> solid:
    # Create three M6x20 bolts spaced apart
    bolt1: solid = metric_hex_bolt("M3", 20.0)

    bolt2: solid = metric_hex_bolt("M6", 20.0)
    bolt2_pos: solid = translate(bolt2, 20.0, 0.0, 0.0)

    bolt3: solid = metric_hex_bolt("M8", 20.0)
    bolt3_pos: solid = translate(bolt3, 40.0, 0.0, 0.0)

    # Combine all
    assembly: solid = compound(bolt1, bolt2_pos, bolt3_pos)

    emit assembly
