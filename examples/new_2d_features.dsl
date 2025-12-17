# New 2D Features Example
# Demonstrates ellipses, splines, 2D booleans, and curve sampling

module new_2d_features

# =============================================================================
# Example 1: Ellipse Curves
# =============================================================================

# Create an ellipse and sample points along it
command DEMO_ELLIPSE() -> ellipse:
    # Create an ellipse centered at origin
    # Semi-major axis 50mm along X, semi-minor 30mm along Y
    center: point = point(0.0, 0.0)
    e: ellipse = ellipse(center, 50.0, 30.0)

    # Sample some points along the ellipse
    p0: point = sample_curve(e, 0.0)
    p1: point = sample_curve(e, 0.25)
    p2: point = sample_curve(e, 0.5)
    p3: point = sample_curve(e, 0.75)

    print("Ellipse sampled at t=0:", p0)
    print("Ellipse sampled at t=0.25:", p1)
    print("Ellipse sampled at t=0.5:", p2)
    print("Ellipse sampled at t=0.75:", p3)

    emit e

# Rotated ellipse with partial arc
command DEMO_ELLIPSE_ARC() -> ellipse:
    center: point = point(0.0, 0.0)
    # Ellipse rotated 45 degrees, arc from 0 to 180 degrees
    e: ellipse = ellipse(center, 40.0, 20.0, 45.0, 0.0, 180.0)
    emit e

# =============================================================================
# Example 2: Catmull-Rom Splines
# =============================================================================

# Create a smooth curve through control points
command DEMO_CATMULLROM() -> catmullrom:
    # Define control points for an organic curve
    p0: point = point(0.0, 0.0, 0.0)
    p1: point = point(20.0, 30.0, 0.0)
    p2: point = point(50.0, 25.0, 0.0)
    p3: point = point(80.0, 40.0, 0.0)
    p4: point = point(100.0, 10.0, 0.0)

    pts: list<point> = [p0, p1, p2, p3, p4]

    # Create centripetal Catmull-Rom spline (alpha=0.5)
    curve: catmullrom = catmullrom(pts, false, 0.5)

    # Sample multiple points along the curve
    sample0: point = sample_curve(curve, 0.0)
    sample1: point = sample_curve(curve, 0.5)
    sample2: point = sample_curve(curve, 1.0)

    print("Spline at t=0:", sample0)
    print("Spline at t=0.5:", sample1)
    print("Spline at t=1:", sample2)

    # Get the curve length
    length: float = curve_length(curve)
    print("Curve length:", length)

    emit curve

# Closed Catmull-Rom spline (loop)
command DEMO_CLOSED_SPLINE() -> catmullrom:
    # Create a closed smooth curve
    p0: point = point(0.0, 0.0, 0.0)
    p1: point = point(30.0, 50.0, 0.0)
    p2: point = point(60.0, 50.0, 0.0)
    p3: point = point(90.0, 0.0, 0.0)
    p4: point = point(60.0, -30.0, 0.0)
    p5: point = point(30.0, -30.0, 0.0)

    pts: list<point> = [p0, p1, p2, p3, p4, p5]

    # Closed spline with uniform parameterization (alpha=0)
    curve: catmullrom = catmullrom(pts, true, 0.0)
    emit curve

# =============================================================================
# Example 3: NURBS Curves
# =============================================================================

# Create a NURBS curve
command DEMO_NURBS() -> nurbs:
    # Control points for a NURBS curve
    p0: point = point(0.0, 0.0, 0.0)
    p1: point = point(25.0, 50.0, 0.0)
    p2: point = point(75.0, 50.0, 0.0)
    p3: point = point(100.0, 0.0, 0.0)

    pts: list<point> = [p0, p1, p2, p3]

    # Uniform weights for NURBS (all 1.0)
    weights: list<float> = [1.0, 1.0, 1.0, 1.0]

    # Degree 3 (cubic) NURBS with uniform weights
    curve: nurbs = nurbs(pts, weights, 3)

    # Sample points
    start_pt: point = sample_curve(curve, 0.0)
    mid_pt: point = sample_curve(curve, 0.5)
    end_pt: point = sample_curve(curve, 1.0)

    print("NURBS start:", start_pt)
    print("NURBS mid:", mid_pt)
    print("NURBS end:", end_pt)

    emit curve

# =============================================================================
# Example 4: Polygon from Points
# =============================================================================

# Create a polygon from arbitrary points
command DEMO_POLYGON() -> region2d:
    # Define vertices of a custom polygon (star shape)
    outer_r: float = 50.0
    inner_r: float = 20.0

    # 5-pointed star vertices (alternating outer/inner)
    angle0: float = radians(90.0)
    angle1: float = radians(126.0)
    angle2: float = radians(162.0)
    angle3: float = radians(198.0)
    angle4: float = radians(234.0)
    angle5: float = radians(270.0)
    angle6: float = radians(306.0)
    angle7: float = radians(342.0)
    angle8: float = radians(18.0)
    angle9: float = radians(54.0)

    # Outer points
    p0: point = point(outer_r * cos(angle0), outer_r * sin(angle0))
    p2: point = point(outer_r * cos(angle2), outer_r * sin(angle2))
    p4: point = point(outer_r * cos(angle4), outer_r * sin(angle4))
    p6: point = point(outer_r * cos(angle6), outer_r * sin(angle6))
    p8: point = point(outer_r * cos(angle8), outer_r * sin(angle8))

    # Inner points
    p1: point = point(inner_r * cos(angle1), inner_r * sin(angle1))
    p3: point = point(inner_r * cos(angle3), inner_r * sin(angle3))
    p5: point = point(inner_r * cos(angle5), inner_r * sin(angle5))
    p7: point = point(inner_r * cos(angle7), inner_r * sin(angle7))
    p9: point = point(inner_r * cos(angle9), inner_r * sin(angle9))

    pts: list<point> = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9]
    star: region2d = polygon(pts)
    emit star

# =============================================================================
# Example 5: Disk Region
# =============================================================================

# Create a filled circular region
command DEMO_DISK() -> region2d:
    center: point = point(0.0, 0.0)
    # 64-sided approximation for smooth circle
    d: region2d = disk(center, 40.0, 64)

    area_val: float = area(d)
    print("Disk area:", area_val)
    print("Expected (pi*r^2):", pi() * 40.0 * 40.0)

    emit d

# =============================================================================
# Example 6: 2D Boolean Operations
# =============================================================================

# Difference: plate with hole
command DEMO_DIFFERENCE2D() -> region2d:
    # Outer plate
    plate: region2d = rectangle(100.0, 60.0)

    # Center hole
    hole: region2d = disk(point(0.0, 0.0), 20.0, 32)

    # Subtract hole from plate
    result: region2d = difference2d(plate, hole)
    emit result

# Union: combine two overlapping shapes
command DEMO_UNION2D() -> region2d:
    # Two overlapping rectangles forming an L-shape
    rect1: region2d = rectangle(60.0, 20.0)
    rect2: region2d = rectangle(20.0, 60.0)

    # Position rect2 to form L (translate not available for region2d in DSL)
    # For now, use overlapping centered rectangles
    combined: region2d = union2d(rect1, rect2)
    emit combined

# Intersection: overlapping area only
command DEMO_INTERSECTION2D() -> region2d:
    # Two overlapping rectangles
    rect1: region2d = rectangle(80.0, 40.0)
    rect2: region2d = rectangle(40.0, 80.0)

    # Keep only the overlapping region (40x40 square)
    overlap: region2d = intersection2d(rect1, rect2)
    emit overlap

# =============================================================================
# Example 7: Complex Profile with Multiple Operations
# =============================================================================

# Create a mounting plate profile with multiple features
# Demonstrates chained difference2d operations with proper hole accumulation
command MAKE_MOUNTING_PLATE_PROFILE() -> region2d:
    # Main plate outline
    plate_width: float = 120.0
    plate_height: float = 80.0
    plate: region2d = rectangle(plate_width, plate_height)

    # Center mounting hole
    center_hole: region2d = disk(point(0.0, 0.0), 15.0, 32)

    # Corner mounting holes
    corner_offset_x: float = 45.0
    corner_offset_y: float = 30.0
    corner_radius: float = 5.0

    hole_tl: region2d = disk(point(0.0 - corner_offset_x, corner_offset_y), corner_radius, 24)
    hole_tr: region2d = disk(point(corner_offset_x, corner_offset_y), corner_radius, 24)
    hole_bl: region2d = disk(point(0.0 - corner_offset_x, 0.0 - corner_offset_y), corner_radius, 24)
    hole_br: region2d = disk(point(corner_offset_x, 0.0 - corner_offset_y), corner_radius, 24)

    # Chain difference operations to create all 5 holes
    # difference2d now properly accumulates holes from previous operations
    step1: region2d = difference2d(plate, center_hole)
    step2: region2d = difference2d(step1, hole_tl)
    step3: region2d = difference2d(step2, hole_tr)
    step4: region2d = difference2d(step3, hole_bl)
    result: region2d = difference2d(step4, hole_br)

    print("Mounting plate profile with 5 holes (chained difference2d)")
    emit result

# =============================================================================
# Example 8: Curve Sampling for Analysis
# =============================================================================

# Sample an ellipse at multiple points and print coordinates
command ANALYZE_ELLIPSE(semi_major: float = 50.0, semi_minor: float = 30.0, num_samples: int = 8) -> ellipse:
    center: point = point(0.0, 0.0)
    e: ellipse = ellipse(center, semi_major, semi_minor)

    # Calculate expected perimeter (Ramanujan approximation)
    h: float = pow(semi_major - semi_minor, 2.0) / pow(semi_major + semi_minor, 2.0)
    approx_perim: float = pi() * (semi_major + semi_minor) * (1.0 + 3.0 * h / (10.0 + sqrt(4.0 - 3.0 * h)))
    print("Approximate perimeter (Ramanujan):", approx_perim)

    # Sample points around ellipse
    print("Sampling", num_samples, "points around ellipse:")

    # Manual unrolled loop for 8 samples
    t0: float = 0.0
    t1: float = 0.125
    t2: float = 0.25
    t3: float = 0.375
    t4: float = 0.5
    t5: float = 0.625
    t6: float = 0.75
    t7: float = 0.875

    pt0: point = sample_curve(e, t0)
    pt1: point = sample_curve(e, t1)
    pt2: point = sample_curve(e, t2)
    pt3: point = sample_curve(e, t3)
    pt4: point = sample_curve(e, t4)
    pt5: point = sample_curve(e, t5)
    pt6: point = sample_curve(e, t6)
    pt7: point = sample_curve(e, t7)

    print("  t=0.000:", pt0)
    print("  t=0.125:", pt1)
    print("  t=0.250:", pt2)
    print("  t=0.375:", pt3)
    print("  t=0.500:", pt4)
    print("  t=0.625:", pt5)
    print("  t=0.750:", pt6)
    print("  t=0.875:", pt7)

    emit e
