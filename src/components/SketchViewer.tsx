/**
 * SketchViewer — Canvas2D renderer for yapCAD sketch entities.
 * Renders polylines (sampled curves) and control points for splines.
 * Supports interactive dragging of control points when interactive=true.
 */

import { useRef, useEffect, useCallback } from 'react';
import type { LoadedSketch } from '../yapcad-loader';

// Re-export as SketchEntity for App.tsx compatibility
export type SketchEntity = LoadedSketch;

interface SketchViewerProps {
  sketch: SketchEntity;
  width?: number;
  height?: number;
  interactive?: boolean;
  onControlPointsChanged?: (points: number[][]) => void;
}

const PAD = 32;           // canvas padding in px
const CTRL_RADIUS = 6;    // control point handle radius
const HIT_RADIUS = 10;    // hit test radius for dragging

export function SketchViewer({
  sketch,
  width = 600,
  height = 320,
  interactive = false,
  onControlPointsChanged,
}: SketchViewerProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  // Mutable state for drag — kept in a ref to avoid re-render on every mousemove
  const dragState = useRef<{
    dragging: boolean;
    pointIdx: number;
    ctrlPts: number[][];
    toCanvas: (p: number[]) => [number, number];
    fromCanvas: (x: number, y: number) => number[];
  } | null>(null);

  // ── coordinate helpers ────────────────────────────────────────────────────

  function buildTransform(bbox: number[] | undefined, w: number, h: number) {
    // bbox = [xmin, ymin, xmax, ymax]
    const xmin = bbox?.[0] ?? 0;
    const ymin = bbox?.[1] ?? 0;
    const xmax = bbox?.[2] ?? 1;
    const ymax = bbox?.[3] ?? 1;
    const dw = xmax - xmin || 1;
    const dh = ymax - ymin || 1;
    const drawW = w - PAD * 2;
    const drawH = h - PAD * 2;
    const scale = Math.min(drawW / dw, drawH / dh);
    // centre the drawing
    const offX = PAD + (drawW - dw * scale) / 2;
    const offY = PAD + (drawH - dh * scale) / 2;

    const toCanvas = (p: number[]): [number, number] => [
      offX + (p[0] - xmin) * scale,
      // flip Y: canvas Y grows down, model Y grows up
      h - (offY + (p[1] - ymin) * scale),
    ];
    const fromCanvas = (cx: number, cy: number): number[] => [
      (cx - offX) / scale + xmin,
      ((h - cy) - offY) / scale + ymin,
    ];
    return { toCanvas, fromCanvas, scale, xmin, ymin, xmax, ymax };
  }

  // ── draw ─────────────────────────────────────────────────────────────────

  const draw = useCallback((ctrlOverride?: number[][]) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const w = canvas.width;
    const h = canvas.height;
    const { toCanvas, xmin, ymin, xmax, ymax } = buildTransform(sketch.boundingBox, w, h);

    // Background
    ctx.fillStyle = '#0d0d1a';
    ctx.fillRect(0, 0, w, h);

    // Subtle grid
    ctx.strokeStyle = '#1e1e32';
    ctx.lineWidth = 1;
    const gridStep = 10;
    const { toCanvas: tc0 } = buildTransform(sketch.boundingBox, w, h);
    for (let x = Math.ceil(xmin / gridStep) * gridStep; x <= xmax; x += gridStep) {
      const [cx] = tc0([x, ymin]);
      ctx.beginPath(); ctx.moveTo(cx, PAD); ctx.lineTo(cx, h - PAD); ctx.stroke();
    }
    for (let y = Math.ceil(ymin / gridStep) * gridStep; y <= ymax; y += gridStep) {
      const [, cy] = tc0([xmin, y]);
      ctx.beginPath(); ctx.moveTo(PAD, cy); ctx.lineTo(w - PAD, cy); ctx.stroke();
    }

    // Polylines (sampled curve)
    for (const polyline of sketch.polylines) {
      if (!polyline.length) continue;
      ctx.beginPath();
      ctx.strokeStyle = '#00d4ff';
      ctx.lineWidth = 2;
      ctx.lineJoin = 'round';
      const [sx, sy] = toCanvas(polyline[0]);
      ctx.moveTo(sx, sy);
      for (let i = 1; i < polyline.length; i++) {
        const [px, py] = toCanvas(polyline[i]);
        ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // Control points from primitives
    const primitives = sketch.primitives ?? [];
    for (const prim of primitives) {
      if (!['catmullrom', 'nurbs', 'bezier'].includes(prim.kind)) continue;
      const ctrlPts = ctrlOverride ?? (prim.points as number[][]);
      if (!ctrlPts?.length) continue;

      // Dashed polyline connecting control points
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = '#ff8c0066';
      ctx.lineWidth = 1;
      ctx.beginPath();
      const [sx, sy] = toCanvas(ctrlPts[0]);
      ctx.moveTo(sx, sy);
      for (let i = 1; i < ctrlPts.length; i++) {
        const [px, py] = toCanvas(ctrlPts[i]);
        ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.setLineDash([]);

      // Control point handles
      for (let i = 0; i < ctrlPts.length; i++) {
        const [px, py] = toCanvas(ctrlPts[i]);
        ctx.beginPath();
        ctx.arc(px, py, CTRL_RADIUS, 0, Math.PI * 2);
        ctx.fillStyle = '#ff8c00';
        ctx.fill();
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.5;
        ctx.stroke();
      }
    }

    // Bounding box label
    const dw = (xmax - xmin).toFixed(1);
    const dh = (ymax - ymin).toFixed(1);
    ctx.fillStyle = '#556';
    ctx.font = '11px monospace';
    ctx.fillText(`${dw} × ${dh} mm`, PAD, h - 8);

    // Command name
    ctx.fillStyle = '#7c7cf8';
    ctx.font = 'bold 11px monospace';
    ctx.fillText(sketch.name ?? 'sketch', PAD, 18);
  }, [sketch]);

  useEffect(() => { draw(); }, [draw]);

  // ── interaction ───────────────────────────────────────────────────────────

  const getPrimCtrlPts = (): number[][] | null => {
    for (const prim of sketch.primitives ?? []) {
      if (['catmullrom', 'nurbs', 'bezier'].includes(prim.kind) && prim.points?.length) {
        return (prim.points as number[][]).map(p => [...p]);
      }
    }
    return null;
  };

  const handleMouseDown = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!interactive) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const { toCanvas, fromCanvas } = buildTransform(sketch.boundingBox, canvas.width, canvas.height);
    const ctrlPts = getPrimCtrlPts();
    if (!ctrlPts) return;

    for (let i = 0; i < ctrlPts.length; i++) {
      const [cx, cy] = toCanvas(ctrlPts[i]);
      if (Math.hypot(mx - cx, my - cy) <= HIT_RADIUS) {
        dragState.current = { dragging: true, pointIdx: i, ctrlPts, toCanvas, fromCanvas };
        break;
      }
    }
  }, [interactive, sketch]);

  const handleMouseMove = useCallback((e: React.MouseEvent<HTMLCanvasElement>) => {
    if (!dragState.current?.dragging) return;
    const canvas = canvasRef.current;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const { fromCanvas, ctrlPts, pointIdx } = dragState.current;
    const newPos = fromCanvas(mx, my);
    const updated = ctrlPts.map((p, i) => i === pointIdx ? [newPos[0], newPos[1], p[2] ?? 0, p[3] ?? 1] : p);
    dragState.current.ctrlPts = updated;
    draw(updated);
  }, [draw]);

  const handleMouseUp = useCallback(() => {
    if (!dragState.current?.dragging) return;
    const pts = dragState.current.ctrlPts;
    dragState.current = null;
    onControlPointsChanged?.(pts);
  }, [onControlPointsChanged]);

  return (
    <div style={styles.wrapper}>
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        style={{ ...styles.canvas, cursor: interactive ? 'crosshair' : 'default' }}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
      />
    </div>
  );
}

const styles = {
  wrapper: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: '#0d0d1a',
    borderRadius: '6px',
    overflow: 'hidden',
    border: '1px solid #2a2a48',
  } as React.CSSProperties,
  canvas: {
    display: 'block',
    maxWidth: '100%',
  } as React.CSSProperties,
};
