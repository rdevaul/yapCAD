/**
 * SketchViewer — Canvas2D renderer for yapCAD sketch entities.
 * Renders polylines (sampled curves) and control point handles for splines.
 * Supports interactive dragging of control points when interactive=true.
 * Also supports circle radius drag handles for circle_r params.
 */

import { useRef, useEffect, useCallback } from 'react';
import type { LoadedSketch } from '../yapcad-loader';

// Re-export as SketchEntity for App.tsx compatibility
export type SketchEntity = LoadedSketch;

export interface CircleParamInfo {
  paramName: string;
  snapMode?: string;
  label?: string;
  /** "radius" (default) or "diameter" — controls how the drag value maps to the param */
  unit?: 'radius' | 'diameter';
}

interface SketchViewerProps {
  sketch: SketchEntity;
  /** Optional override for control point positions (from live param values) */
  controlPoints?: number[][];
  width?: number;
  height?: number;
  interactive?: boolean;
  onControlPointsChanged?: (points: number[][]) => void;
  /** circle_r param declarations in declaration order */
  circleParams?: CircleParamInfo[];
  /** Called on mouseup after a circle radius drag */
  onCircleRadiusChanged?: (paramName: string, radius: number) => void;
}

const PAD = 40;
const CTRL_RADIUS = 7;
const HIT_RADIUS = 12;
const CIRCLE_HANDLE_SIZE = 9;  // half-diagonal of the diamond
const CIRCLE_HANDLE_HIT = 14;
const SAMPLES_PER_SEGMENT = 16;

// ── Snapping ─────────────────────────────────────────────────────────────────

/** M3-M10 tap drill diameters in mm */
const METRIC_TAP_DIAMETERS = [2.5, 3.3, 4.2, 5.0, 6.75, 8.5];
/** Common unified clearance hole diameters in mm (1/4-20, 5/16-18, 3/8-16, 1/2-13) */
const UNIFIED_CLEARANCE_DIAMETERS = [6.8, 8.4, 10.0, 13.5];

function snapRadius(r: number, snapMode: string | undefined): number {
  const mode = snapMode ?? 'none';
  if (mode === 'mm') {
    // Round to nearest 0.5mm on diameter (0.25mm on radius) — 1mm diameter steps
    const dia = r * 2;
    return Math.round(dia) / 2;
  }
  if (mode === '0.5mm') {
    // Round to nearest 1mm on radius = 2mm on diameter
    return Math.round(r);
  }
  if (mode === 'metric_tap') {
    const dia = r * 2;
    let closest = METRIC_TAP_DIAMETERS[0];
    let bestDist = Math.abs(dia - closest);
    for (const d of METRIC_TAP_DIAMETERS) {
      const dist = Math.abs(dia - d);
      if (dist < bestDist) { bestDist = dist; closest = d; }
    }
    return closest / 2;
  }
  if (mode === 'unified_clearance') {
    const dia = r * 2;
    let closest = UNIFIED_CLEARANCE_DIAMETERS[0];
    let bestDist = Math.abs(dia - closest);
    for (const d of UNIFIED_CLEARANCE_DIAMETERS) {
      const dist = Math.abs(dia - d);
      if (dist < bestDist) { bestDist = dist; closest = d; }
    }
    return closest / 2;
  }
  return r;
}

// ── Native circle primitive extraction ───────────────────────────────────────

interface NativeCircle {
  cx: number;
  cy: number;
  r: number;
  paramName: string;   // from primitive.meta.param
  label?: string;      // from primitive.meta.label
}

/**
 * Extract tagged circle primitives from a sketch.
 * Replaces the old polyline-fitting detectCircles + matchCirclesToParams approach.
 *
 * geometry_json serialises full circles as:
 *   { kind: "circle", center: [x, y, z, w], radius: r, orientation: n, meta: { param: "...", label: "..." } }
 * Note: center and radius are top-level fields (NOT nested under params).
 */
function extractCirclePrimitives(sketch: LoadedSketch): NativeCircle[] {
  if (!sketch.primitives) return [];
  const result: NativeCircle[] = [];
  for (const prim of sketch.primitives) {
    if (prim.kind !== 'circle') continue;
    const paramName = (prim.meta as any)?.param as string | undefined;
    if (!paramName) continue;  // skip untagged circles
    // center and radius are top-level fields in the JSON (not under params)
    const center = (prim as any).center as number[] | undefined;
    const radius  = (prim as any).radius  as number  | undefined;
    if (center == null || radius == null) continue;
    result.push({
      cx: center[0],
      cy: center[1],
      r: radius,
      paramName,
      label: (prim.meta as any)?.label as string | undefined,
    });
  }
  return result;
}

// ── Client-side Catmull-Rom sampler ──────────────────────────────────────────

function lerp(a: number, b: number, t: number, t0: number, t1: number): number {
  const d = t1 - t0;
  if (Math.abs(d) < 1e-10) return a;
  return a + (b - a) * (t - t0) / d;
}

function catmullRomSample(pts: number[][], alpha = 0.5, closed = false): [number, number][] {
  if (pts.length < 2) return pts.map(p => [p[0], p[1]]);

  const P: number[][] = closed
    ? [pts[pts.length - 1], ...pts, pts[0], pts[1]]
    : [pts[0], ...pts, pts[pts.length - 1]];

  const result: [number, number][] = [];

  for (let seg = 0; seg < P.length - 3; seg++) {
    const p0 = P[seg], p1 = P[seg + 1], p2 = P[seg + 2], p3 = P[seg + 3];

    const dt01 = Math.pow(Math.hypot(p1[0] - p0[0], p1[1] - p0[1]), alpha) || 1e-4;
    const dt12 = Math.pow(Math.hypot(p2[0] - p1[0], p2[1] - p1[1]), alpha) || 1e-4;
    const dt23 = Math.pow(Math.hypot(p3[0] - p2[0], p3[1] - p2[1]), alpha) || 1e-4;

    const t0 = 0;
    const t1 = t0 + dt01;
    const t2 = t1 + dt12;
    const t3 = t2 + dt23;

    const n = SAMPLES_PER_SEGMENT;
    for (let i = (seg === 0 ? 0 : 1); i <= n; i++) {
      const t = t1 + (t2 - t1) * (i / n);

      const A1x = lerp(p0[0], p1[0], t, t0, t1);
      const A1y = lerp(p0[1], p1[1], t, t0, t1);
      const A2x = lerp(p1[0], p2[0], t, t1, t2);
      const A2y = lerp(p1[1], p2[1], t, t1, t2);
      const A3x = lerp(p2[0], p3[0], t, t2, t3);
      const A3y = lerp(p2[1], p3[1], t, t2, t3);

      const B1x = lerp(A1x, A2x, t, t0, t2);
      const B1y = lerp(A1y, A2y, t, t0, t2);
      const B2x = lerp(A2x, A3x, t, t1, t3);
      const B2y = lerp(A2y, A3y, t, t1, t3);

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
  const ox = PAD + (drawW - dw * scale) / 2;
  const oy = PAD + (drawH - dh * scale) / 2;

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

// ── drag state types ──────────────────────────────────────────────────────────

type ControlPointDrag = {
  kind: 'controlPoint';
  active: boolean;
  idx: number;
  pts: number[][];
};

type CircleDrag = {
  kind: 'circleDrag';
  active: boolean;
  cx: number;       // model coords
  cy: number;
  originalR: number;
  currentR: number;
  paramName: string;
  snapMode?: string;
};

type DragState = ControlPointDrag | CircleDrag | null;

// ── component ───────────────────────────────────────────────────────────────

export function SketchViewer({
  sketch,
  controlPoints,
  width = 600,
  height = 360,
  interactive = false,
  onControlPointsChanged,
  circleParams,
  onCircleRadiusChanged,
}: SketchViewerProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const dragRef = useRef<DragState>(null);
  const sizeRef = useRef({ w: width, h: height });
  const bboxRef = useRef<number[] | undefined>(undefined);

  // Live circle drag radius — updated every mousemove
  const circleDragLiveRef = useRef<{
    cx: number; cy: number; r: number; paramName: string; snapMode?: string;
  } | null>(null);

  // ── draw ─────────────────────────────────────────────────────────────────

  const draw = useCallback((ctrlOverride?: number[][], circleDragPreview?: {
    cx: number; cy: number; r: number;
  }) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const w = canvas.width;
    const h = canvas.height;
    if (w === 0 || h === 0) return;

    // Compute bbox
    const activePtsForBbox = ctrlOverride ?? controlPoints;
    if (activePtsForBbox && activePtsForBbox.length > 0) {
      const xs = activePtsForBbox.map(p => p[0]);
      const ys = activePtsForBbox.map(p => p[1]);
      const sampled = catmullRomSample(activePtsForBbox);
      const cxs = sampled.map(p => p[0]);
      const cys = sampled.map(p => p[1]);
      bboxRef.current = [
        Math.min(...xs, ...cxs),
        Math.min(...ys, ...cys),
        Math.max(...xs, ...cxs),
        Math.max(...ys, ...cys),
      ];
    } else {
      bboxRef.current = sketch.boundingBox;
    }

    const { toCanvas, xmin, ymin, xmax, ymax, scale } = makeTransform(bboxRef.current, w, h);

    // ── background
    ctx.fillStyle = '#0d0d1a';
    ctx.fillRect(0, 0, w, h);

    // ── grid
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

    // ── extract native circle primitives (replaces polyline-fitting) ─────
    const nativeCircles = extractCirclePrimitives(sketch);
    // Build a Set of paramNames being dragged for quick lookup
    const draggingParam = dragRef.current?.kind === 'circleDrag' ? dragRef.current.paramName : null;

    // ── polylines ──────────────────────────────────────────────────────────
    const activePts = ctrlOverride ?? controlPoints;
    const polylinesToDraw: [number, number][][] =
      activePts && activePts.length >= 2
        ? [catmullRomSample(activePts)]
        : sketch.polylines.map(pl => pl.map(p => [p[0], p[1]] as [number, number]));

    for (let plIdx = 0; plIdx < polylinesToDraw.length; plIdx++) {
      const polyline = polylinesToDraw[plIdx];
      if (polyline.length < 2) continue;

      ctx.beginPath();
      ctx.strokeStyle = '#00d4ff';
      ctx.lineWidth = 2.5;
      ctx.lineJoin = 'round';
      ctx.lineCap = 'round';
      ctx.setLineDash([]);
      const [sx, sy] = toCanvas(polyline[0][0], polyline[0][1]);
      ctx.moveTo(sx, sy);
      for (let i = 1; i < polyline.length; i++) {
        const [px, py] = toCanvas(polyline[i][0], polyline[i][1]);
        ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // ── circle drag preview ring ──────────────────────────────────────────
    if (circleDragPreview) {
      const [pcx, pcy] = toCanvas(circleDragPreview.cx, circleDragPreview.cy);
      const pr = circleDragPreview.r * scale;

      // Dashed preview ring in cyan
      ctx.beginPath();
      ctx.arc(pcx, pcy, pr, 0, Math.PI * 2);
      ctx.strokeStyle = '#00ffcc';
      ctx.lineWidth = 2;
      ctx.setLineDash([6, 4]);
      ctx.stroke();
      ctx.setLineDash([]);

      // Label new radius
      const label = `r=${circleDragPreview.r.toFixed(2)}mm  ⌀${(circleDragPreview.r * 2).toFixed(2)}mm`;
      ctx.fillStyle = '#00ffcc';
      ctx.font = 'bold 11px monospace';
      ctx.textAlign = 'center';
      ctx.fillText(label, pcx, pcy - pr - 8);
    }

    // ── spline control points ─────────────────────────────────────────────
    for (const prim of sketch.primitives ?? []) {
      if (!['catmullrom', 'nurbs', 'bezier'].includes(prim.kind)) continue;
      const rawPts = ctrlOverride ?? controlPoints ?? (prim.points as number[][]);
      if (!rawPts?.length) continue;

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

      for (let i = 0; i < rawPts.length; i++) {
        const [px, py] = toCanvas(rawPts[i][0], rawPts[i][1]);
        ctx.beginPath();
        ctx.arc(px, py, CTRL_RADIUS, 0, Math.PI * 2);
        ctx.fillStyle = '#ff8c00';
        ctx.fill();
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 1.5;
        ctx.stroke();

        ctx.fillStyle = '#fff';
        ctx.font = 'bold 10px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(String(i), px, py);
      }
    }

    // ── native circle primitives — draw outlines + drag handles ──────────
    for (const nc of nativeCircles) {
      const [scx, scy] = toCanvas(nc.cx, nc.cy);
      const sr = nc.r * scale;

      // Draw the circle outline (dashed for bolt-circle param, solid for others)
      ctx.beginPath();
      ctx.arc(scx, scy, sr, 0, Math.PI * 2);
      if (nc.paramName === 'bolt_circle_mm') {
        ctx.setLineDash([8, 5]);
        ctx.strokeStyle = 'rgba(150, 150, 255, 0.7)';
      } else {
        ctx.setLineDash([]);
        ctx.strokeStyle = 'rgba(80, 180, 255, 0.85)';
      }
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.setLineDash([]);

      // Diamond drag handle at 3-o'clock
      const isDragging = draggingParam === nc.paramName;
      const [hx, hy] = toCanvas(nc.cx + nc.r, nc.cy);
      const s = CIRCLE_HANDLE_SIZE;
      ctx.beginPath();
      ctx.moveTo(hx,     hy - s);
      ctx.lineTo(hx + s, hy);
      ctx.lineTo(hx,     hy + s);
      ctx.lineTo(hx - s, hy);
      ctx.closePath();
      ctx.fillStyle = isDragging ? '#00ffcc' : '#ff8c00';
      ctx.fill();
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 1.5;
      ctx.stroke();

      // "R" label on handle
      ctx.fillStyle = '#fff';
      ctx.font = 'bold 9px monospace';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('R', hx, hy);
    }

    // ── info label ─────────────────────────────────────────────────────────
    const dw = (xmax - xmin).toFixed(1);
    const dh_val = (ymax - ymin).toFixed(1);
    ctx.setLineDash([]);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#445';
    ctx.font = '11px monospace';
    ctx.fillText(`${dw} × ${dh_val} mm  |  scale: ${scale.toFixed(2)}px/mm`, PAD, h - 10);

    const kind = sketch.primitives?.[0]?.kind ?? 'sketch';
    ctx.fillStyle = '#7c7cf8';
    ctx.font = 'bold 12px monospace';
    ctx.fillText(kind.toUpperCase(), PAD, 20);

  }, [sketch, circleParams]);

  useEffect(() => { draw(); }, [draw]);

  // ── resize observer ───────────────────────────────────────────────────────
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
    const { toCanvas, fromCanvas } = makeTransform(bboxRef.current, canvas.width, canvas.height);

    // ── Check native circle drag handles first ────────────────────────────
    if (onCircleRadiusChanged) {
      const nativeCircles = extractCirclePrimitives(sketch);
      for (const nc of nativeCircles) {
        const [hx, hy] = toCanvas(nc.cx + nc.r, nc.cy);
        if (Math.hypot(mx - hx, my - hy) <= CIRCLE_HANDLE_HIT) {
          // Look up snapMode from circleParams prop (if provided)
          const cp = circleParams?.find(p => p.paramName === nc.paramName);
          dragRef.current = {
            kind: 'circleDrag',
            active: true,
            cx: nc.cx,
            cy: nc.cy,
            originalR: nc.r,
            currentR: nc.r,
            paramName: nc.paramName,
            snapMode: cp?.snapMode,
          };
          circleDragLiveRef.current = {
            cx: nc.cx, cy: nc.cy, r: nc.r,
            paramName: nc.paramName, snapMode: cp?.snapMode,
          };
          e.preventDefault();
          return;
        }
      }
    }

    // ── Check spline control points ───────────────────────────────────────
    const pts = getCtrlPts();
    if (!pts) return;
    for (let i = 0; i < pts.length; i++) {
      const [cx, cy] = toCanvas(pts[i][0], pts[i][1]);
      if (Math.hypot(mx - cx, my - cy) <= HIT_RADIUS) {
        dragRef.current = { kind: 'controlPoint', active: true, idx: i, pts };
        e.preventDefault();
        return;
      }
    }
    void fromCanvas;
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!dragRef.current?.active) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (canvas.width / rect.width);
    const my = (e.clientY - rect.top) * (canvas.height / rect.height);
    const { fromCanvas } = makeTransform(bboxRef.current, canvas.width, canvas.height);

    if (dragRef.current.kind === 'circleDrag') {
      const drag = dragRef.current;
      const [nx, ny] = fromCanvas(mx, my);
      let newR = Math.hypot(nx - drag.cx, ny - drag.cy);
      newR = Math.max(0.1, newR);
      newR = snapRadius(newR, drag.snapMode);

      drag.currentR = newR;
      circleDragLiveRef.current = { cx: drag.cx, cy: drag.cy, r: newR, paramName: drag.paramName, snapMode: drag.snapMode };

      // Redraw with preview
      draw(undefined, { cx: drag.cx, cy: drag.cy, r: newR });
      return;
    }

    if (dragRef.current.kind === 'controlPoint') {
      const [nx, ny] = fromCanvas(mx, my);
      const { idx, pts } = dragRef.current;
      const updated = pts.map((p, i) =>
        i === idx ? [nx, ny, p[2] ?? 0, p[3] ?? 1] : p
      );
      dragRef.current.pts = updated;
      draw(updated);
    }
  };

  const handleMouseUp = () => {
    if (!dragRef.current?.active) return;

    if (dragRef.current.kind === 'circleDrag') {
      const live = circleDragLiveRef.current;
      if (live) {
        onCircleRadiusChanged?.(live.paramName, live.r);
      }
      dragRef.current = null;
      circleDragLiveRef.current = null;
      draw(); // redraw without preview
      return;
    }

    if (dragRef.current.kind === 'controlPoint') {
      const pts = dragRef.current.pts;
      dragRef.current = null;
      onControlPointsChanged?.(pts);
      return;
    }

    dragRef.current = null;
  };

  // ── cursor style ──────────────────────────────────────────────────────────
  const getCursor = () => {
    if (!interactive) return 'default';
    if (dragRef.current?.kind === 'circleDrag') return 'ew-resize';
    return 'crosshair';
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
          cursor: getCursor(),
        }}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
      />
    </div>
  );
}
