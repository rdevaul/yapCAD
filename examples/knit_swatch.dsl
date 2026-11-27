module knit_swatch

# =============================================================================
#  A knitted swatch
#
#  A patch of stockinette with a heart worked into it in purl stitches and a
#  curved outline.  Nothing is modelled stitch by stitch: one stitch is two
#  capsules, and the fabric is that stitch repeated without end and cut to
#  the outline.
#
#      python -m yapcad.render examples/knit_swatch.dsl SWATCH --view top
#      python -m yapcad.dsl run examples/knit_swatch.dsl SWATCH \
#          --representation sdf --sdf-cell 0.12 --output swatch.stl
#
#  The cell has to be finer than the gaps between the yarns, or the mesher
#  finds two surfaces in one cell; 0.12 mm suits this gauge.
# =============================================================================

command stitch_w() -> float:
    emit 3.0

command stitch_h() -> float:
    emit 2.3

command yarn_r() -> float:
    emit 0.58

# A knit stitch: the two legs of a V, centred on its cell.  The legs reach
# into the rows above and below, which is what makes the columns of Vs nest.
command knit_stitch(z: float) -> solid:
    let w: float = stitch_w()
    let h: float = stitch_h()
    let left: solid = capsule(point(0.0 - 0.10 * w, 0.0 - 0.56 * h, z - 0.20),
                              point(0.0 - 0.30 * w, 0.62 * h, z + 0.12),
                              yarn_r())
    let right: solid = capsule(point(0.10 * w, 0.0 - 0.56 * h, z - 0.20),
                               point(0.30 * w, 0.62 * h, z + 0.12), yarn_r())
    emit smooth_union(left, right, 0.25)

# A purl stitch: the head of the loop, a bump lying across its cell.
command purl_stitch(z: float) -> solid:
    let w: float = stitch_w()
    let h: float = stitch_h()
    let left: solid = capsule(point(0.0 - 0.28 * w, 0.0 - 0.10 * h, z),
                              point(0.0, 0.10 * h, z + 0.28), yarn_r())
    let right: solid = capsule(point(0.28 * w, 0.0 - 0.10 * h, z),
                               point(0.0, 0.10 * h, z + 0.28), yarn_r())
    emit smooth_union(left, right, 0.2)

# The motif as a knitting chart, top row first: 1 is purl, 0 is knit.
command heart() -> list<list<int>>:
    emit [
        [0, 1, 1, 0, 1, 1, 0],
        [1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 1, 1, 0],
        [0, 0, 1, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0]
    ]

# The cells a chart marks, as tall boxes; the chart's top-left stitch is the
# cell in column `col`, row `row`.
command chart_cells(chart: list<list<int>>, col: int, row: int) -> solid:
    let cells: list<solid> = []
    for r in range(len(chart)):
        for c in range(len(chart[r])):
            if chart[r][c] == 1:
                cells = append(cells, translate(
                    box(stitch_w(), stitch_h(), 12.0),
                    (col + c) * stitch_w(), (row - r) * stitch_h(), 0.0))
    emit union_all(cells)

command SWATCH(
    backing: float @ui(widget="slider", label="Backing thickness",
                       min=1.0, max=3.0) = 1.6
) -> solid:
    let edge: catmullrom = catmullrom([
        point(-24.0, -17.0), point(0.0, -20.0), point(24.0, -17.0),
        point(27.0, 0.0), point(24.0, 17.0), point(0.0, 20.0),
        point(-24.0, 17.0), point(-27.0, 0.0)], true, 0.5)
    let profile: region2d = region_from_spline(edge, 240)
    let outline: solid = translate(extrude_field(profile, 20.0), 0.0, 0.0, -10.0)
    let motif: solid = chart_cells(heart(), -3, 3)
    let z: float = backing + 0.1
    # Knit everywhere but the motif; the blend rounds the stitches it cuts.
    let knit: solid = smooth_difference(
        repeat_fill(knit_stitch(z), outline, stitch_w(), stitch_h(), 0.0),
        motif, 0.3)
    let purl: solid = repeat_fill(purl_stitch(z), motif,
                                  stitch_w(), stitch_h(), 0.0)
    emit union_all([extrude_field(profile, backing), knit, purl])
