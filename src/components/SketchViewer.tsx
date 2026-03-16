/**
 * SketchViewer — Canvas2D renderer for yapCAD sketch entities.
 * Renders polylines (sampled curves) and control point handles for splines.
 * Supports interactive dragging of control points when interactive=true.
 */

import { useRef, useEffect, useCallback } from 'react';
import type { LoadedSketch } from '../yapcad-loader';

// Re-export as SketchEntity for App.tsx compatibility
export type SketchEntity = LoadedSketch;

interface SketchViewerProps {
  sketch: SketchEntity;
  /** Optional override for control point positions (from live param values) */
  controlPoints?: number[][];
  width?: number;
  height?: number;
  interactive?: boolean;
  onControlPointsChanged?: (points: number[][]) => void;
}

const PAD = 40;
const CTRL_RADIUS = 7;
const HIT_RADIUS = 12;
const SAMPLES_PER_SEGMENT = 16;

// ── Client-side Catmull-Rom sampler ──────────────────────────────────────────
// Centripetal (alpha=0.5) parametrisation — matches yapCAD default.
// Uses the Barry-Goldman algorithm with correct knot intervals.

function lerp(a: number, b: number, t: number, t0: number, t1: number): number {
  const d = t1 - t0;
  if (Math.abs(d) < 1e-10) return a;
  return a + (b - a) * (t - t0) / d;
}

function catmullRomSample(pts: number[][], alpha = 0.5, closed = false): [number, number][] {
  if (pts.length < 2) return pts.map(p => [p[0], p[1]]);

  // Duplicate endpoints so every original point is interior to a 4-point window
  const P: number[][] = closed
    ? [pts[pts.length - 1], ...pts, pts[0], pts[1]]
    : [pts[0], ...pts, pts[pts.length - 1]];

  const result: [number, number][] = [];

  for (let seg = 0; seg < P.length - 3; seg++) {
    const p0 = P[seg], p1 = P[seg + 1], p2 = P[seg + 2], p3 = P[seg + 3];

    // Centripetal knot intervals: dt = dist^alpha
    const dt01 = Math.pow(Math.hypot(p1[0] - p0[0], p1[1] - p0[1]), alpha) || 1e-4;
    const dt12 = Math.pow(Math.hypot(p2[0] - p1[0], p2[1] - p1[1]), alpha) || 1e-4;
    const dt23 = Math.pow(Math.hypot(p3[0] - p2[0], p3[1] - p2[1]), alpha) || 1e-4;

    // Knot values (accumulate from 0 at p0)
    const t0 = 0;
    const t1 = t0 + dt01;
    const t2 = t1 + dt12;
    const t3 = t2 + dt23;

    const n = SAMPLES_PER_SEGMENT;
    // Parametric range for this segment is [t1, t2]
    for (let i = (seg === 0 ? 0 : 1); i <= n; i++) {
      const t = t1 + (t2 - t1) * (i / n);

      // Level 1: lerp adjacent pairs
      const A1x = lerp(p0[0], p1[0], t, t0, t1);
      const A1y = lerp(p0[1], p1[1], t, t0, t1);
      const A2x = lerp(p1[0], p2[0], t, t1, t2);
      const A2y = lerp(p1[1], p2[1], t, t1, t2);
      const A3x = lerp(p2[0], p3[0], t, t2, t3);
      const A3y = lerp(p2[1], p3[1], t, t2, t3);

      // Level 2: lerp A values using outer intervals
      const B1x = lerp(A1x, A2x, t, t0, t2);
      const B1y = lerp(A1y, A2y, t, t0, t2);
      const B2x = lerp(A2x, A3x, t, t1, t3);
      const B2y = lerp(A2y, A3y, t, t1, t3);

      // Level 3: final lerp over the segment interval
      const Cx = lerp(B1x, B2x, t, t1, t2);
      const Cy = lerp(B1y, B2y, t, t1, t2);

      result.push([Cx, Cy]);
    }
  }

  return result;
}

// ── coordinate transform ────────────────────────────────────────────────────

function makeTransform(bbox: number[] | undefined, w: number, h: number) {
  const xmin = bbox?.[0] ?? 0;
  const ymin = bbox?.[1] ?? 0;
  const xmax = bbox?.[2] ?? 1;
  const ymax = bbox?.[3] ?? 1;
  const dw = xmax - xmin || 1;
  const dh = ymax - ymin || 1;
  const drawW = w - PAD * 2;
  const drawH = h - PAD * 2;
  const scale = Math.min(drawW / dw, drawH / dh);
  // Centre the drawing in the available space
  const ox = PAD + (drawW - dw * scale) / 2;
  const oy = PAD + (drawH - dh * scale) / 2;

  // Canvas: Y grows DOWN. Model: Y grows UP. Flip around canvas centre.
  const toCanvas = (x: number, y: number): [number, number] => [
    ox + (x - xmin) * scale,
    h - oy - (y - ymin) * scale,
  ];
  const fromCanvas = (cx: number, cy: number): [number, number] => [
    (cx - ox) / scale + xmin,
    (h - oy - cy) / scale + ymin,
  ];

  return { toCanvas, fromCanvas, scale, xmin, ymin, xmax, ymax };
}

// ── component ───────────────────────────────────────────────────────────────

export function SketchViewer({
  sketch,
  controlPoints,
  width = 600,
  height = 360,
  interactive = false,
  onControlPointsChanged,
}: SketchViewerProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const dragRef = useRef<{
    active: boolean;
    idx: number;
    pts: number[][];
  } | null>(null);
  // Track current canvas buffer size independently of props
  const sizeRef = useRef({ w: width, h: height });

  // ── draw ─────────────────────────────────────────────────────────────────

  const draw = useCallback((ctrlOverride?: number[][]) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const w = canvas.width;
    const h = canvas.height;

    if (w === 0 || h === 0) return;

    // When live control points are active, compute bbox from them so the
    // canvas auto-scales to fit — avoids clipping when points move outside
    // the stored bbox from the last eval.
    const activePtsForBbox = ctrlOverride ?? controlPoints;
    let liveBbox = sketch.boundingBox;
    if (activePtsForBbox && activePtsForBbox.length > 0) {
      const xs = activePtsForBbox.map(p => p[0]);
      const ys = activePtsForBbox.map(p => p[1]);
      // Also include the sampled curve to account for spline overshoot
      const sampled = catmullRomSample(activePtsForBbox);
      const cxs = sampled.map(p => p[0]);
      const cys = sampled.map(p => p[1]);
      liveBbox = [
        Math.min(...xs, ...cxs),
        Math.min(...ys, ...cys),
        Math.max(...xs, ...cxs),
        Math.max(...ys, ...cys),
      ];
    }

    const { toCanvas, xmin, ymin, xmax, ymax, scale } = makeTransform(liveBbox, w, h);

    // ── background
    ctx.fillStyle = '#0d0d1a';
    ctx.fillRect(0, 0, w, h);

    // ── subtle grid
    ctx.strokeStyle = '#1e1e38';
    ctx.lineWidth = 1;
    const rawStep = (xmax - xmin) / 5;
    const gridStep = Math.pow(10, Math.floor(Math.log10(rawStep)));
    for (let x = Math.ceil(xmin / gridStep) * gridStep; x <= xmax + gridStep * 0.001; x += gridStep) {
      const [cx] = toCanvas(x, ymin);
      ctx.beginPath(); ctx.moveTo(cx, PAD); ctx.lineTo(cx, h - PAD); ctx.stroke();
    }
    for (let y = Math.ceil(ymin / gridStep) * gridStep; y <= ymax + gridStep * 0.001; y += gridStep) {
      const [, cy] = toCanvas(xmin, y);
      ctx.beginPath(); ctx.moveTo(PAD, cy); ctx.lineTo(w - PAD, cy); ctx.stroke();
    }

    // ── axes
    ctx.strokeStyle = '#2a2a50';
    ctx.lineWidth = 1.5;
    const [ax0] = toCanvas(xmin, 0);
    const [ax1] = toCanvas(xmax, 0);
    const [, ay] = toCanvas(0, 0);
    if (ay >= PAD && ay <= h - PAD) {
      ctx.beginPath(); ctx.moveTo(PAD, ay); ctx.lineTo(w - PAD, ay); ctx.stroke();
    }
    void ax0; void ax1;

    // ── polylines (sampled curve) ──────────────────────────────────────────
    // If we have live control points, re-sample client-side so the curve
    // updates immediately during drag and after param edits (before re-eval).
    const activePts = ctrlOverride ?? controlPoints;
    const polylinesToDraw: [number, number][][] =
      activePts && activePts.length >= 2
        ? [catmullRomSample(activePts)]
        : sketch.polylines.map(pl => pl.map(p => [p[0], p[1]] as [number, number]));

    for (const polyline of polylinesToDraw) {
      if (polyline.length < 2) continue;
      ctx.beginPath();
      ctx.strokeStyle = '#00d4ff';
      ctx.lineWidth = 2.5;
      ctx.lineJoin = 'round';
      ctx.lineCap = 'round';
      const [sx, sy] = toCanvas(polyline[0][0], polyline[0][1]);
      ctx.moveTo(sx, sy);
      for (let i = 1; i < polyline.length; i++) {
        const [px, py] = toCanvas(polyline[i][0], polyline[i][1]);
        ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // ── control points ────────────────────────────────────────────────────
    for (const prim of sketch.primitives ?? []) {
      if (!['catmullrom', 'nurbs', 'bezier'].includes(prim.kind)) continue;
      // Priority: drag override > external controlPoints prop > eval primitives
      const rawPts = ctrlOverride ?? controlPoints ?? (prim.points as number[][]);
      if (!rawPts?.length) continue;

      // Dashed hull line
      ctx.setLineDash([5, 4]);
      ctx.strokeStyle = 'rgba(255, 140, 0, 0.4)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      const [hx0, hy0] = toCanvas(rawPts[0][0], rawPts[0][1]);
      ctx.moveTo(hx0, hy0);
      for (let i = 1; i < rawPts.length; i++) {
        const [hx, hy] = toCanvas(rawPts[i][0], rawPts[i][1]);
        ctx.lineTo(hx, hy);
      }
      ctx.stroke();
      ctx.setLineDash([]);

      // Handles
      for (let i = 0; i < rawPts.length; i++) {
        const [px, py] = toCanvas(rawPts[i][0], rawPts[i][1]);
        ctx.beginPath();
        ctx.arc(px, py, CTRL_RADIUS, 0, Math.PI * 2);
        ctx.fillStyle = '#ff8c00';
        ctx.fill();
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 1.5;
        ctx.stroke();

        // Index label
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 10px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(String(i), px, py);
      }
    }

    // ── info label ─────────────────────────────────────────────────────────
    const dw = (xmax - xmin).toFixed(1);
    const dh_val = (ymax - ymin).toFixed(1);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#445';
    ctx.font = '11px monospace';
    ctx.fillText(`${dw} × ${dh_val} mm  |  scale: ${scale.toFixed(2)}px/mm`, PAD, h - 10);

    // kind label
    const kind = sketch.primitives?.[0]?.kind ?? 'sketch';
    ctx.fillStyle = '#7c7cf8';
    ctx.font = 'bold 12px monospace';
    ctx.fillText(kind.toUpperCase(), PAD, 20);

  }, [sketch]);

  useEffect(() => { draw(); }, [draw]);

  // ── resize observer — update canvas buffer size and redraw ───────────────
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ro = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width: w, height: h } = entry.contentRect;
        if (w > 0 && h > 0) {
          canvas.width = Math.round(w);
          canvas.height = Math.round(h);
          sizeRef.current = { w: canvas.width, h: canvas.height };
          draw();
        }
      }
    });
    ro.observe(canvas);
    // Set initial size from actual rendered dimensions
    const { clientWidth, clientHeight } = canvas;
    if (clientWidth > 0 && clientHeight > 0) {
      canvas.width = clientWidth;
      canvas.height = clientHeight;
      sizeRef.current = { w: clientWidth, h: clientHeight };
    }
    draw();
    return () => ro.disconnect();
  }, [draw]);

  // ── interaction ───────────────────────────────────────────────────────────

  const getCtrlPts = (): number[][] | null => {
    // Prefer externally provided control points (from live param values)
    if (controlPoints && controlPoints.length > 0) {
      return controlPoints.map(p => [...p]);
    }
    for (const prim of sketch.primitives ?? []) {
      if (['catmullrom', 'nurbs', 'bezier'].includes(prim.kind) && prim.points?.length) {
        return (prim.points as number[][]).map(p => [...p]);
      }
    }
    return null;
  };

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!interactive) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (canvas.width / rect.width);
    const my = (e.clientY - rect.top) * (canvas.height / rect.height);
    const { toCanvas } = makeTransform(sketch.boundingBox, canvas.width, canvas.height);
    const pts = getCtrlPts();
    if (!pts) return;
    for (let i = 0; i < pts.length; i++) {
      const [cx, cy] = toCanvas(pts[i][0], pts[i][1]);
      if (Math.hypot(mx - cx, my - cy) <= HIT_RADIUS) {
        dragRef.current = { active: true, idx: i, pts };
        e.preventDefault();
        return;
      }
    }
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!dragRef.current?.active) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (canvas.width / rect.width);
    const my = (e.clientY - rect.top) * (canvas.height / rect.height);
    const { fromCanvas } = makeTransform(sketch.boundingBox, canvas.width, canvas.height);
    const [nx, ny] = fromCanvas(mx, my);
    const { idx, pts } = dragRef.current;
    const updated = pts.map((p, i) =>
      i === idx ? [nx, ny, p[2] ?? 0, p[3] ?? 1] : p
    );
    dragRef.current.pts = updated;
    draw(updated);
  };

  const handleMouseUp = () => {
    if (!dragRef.current?.active) return;
    const pts = dragRef.current.pts;
    dragRef.current = null;
    onControlPointsChanged?.(pts);
  };

  return (
    <div style={{
      width: '100%', height: '100%',
      display: 'flex', alignItems: 'stretch', justifyContent: 'stretch',
      backgroundColor: '#0d0d1a',
    }}>
      <canvas
        ref={canvasRef}
        style={{
          display: 'block',
          width: '100%',
          height: '100%',
          cursor: interactive ? 'crosshair' : 'default',
        }}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
      />
    </div>
  );
}
