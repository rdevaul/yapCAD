/**
 * yapCAD Web Viewer - Main Application
 * Three-column layout: Editor | Viewer | Packages
 */

import { useState, useCallback, useRef, useEffect } from 'react';
import { PackageSelector } from './components/PackageSelector';
import { LightingSelector } from './components/LightingSelector';
import { RenderModeSelector } from './components/RenderModeSelector';
import { ClippingControls } from './components/ClippingControls';
import { ResizeHandle } from './components/ResizeHandle';
import { LayerPanel } from './components/LayerPanel';
import { DSLEditor } from './components/DSLEditor';
import { TabBar } from './components/TabBar';
import { ParameterPanel } from './components/ParameterPanel';
import { ChatPanel } from './components/ChatPanel';
import { SessionBar } from './components/SessionBar';
import { FileBar } from './components/FileBar';
import { SketchViewer, type SketchEntity } from './components/SketchViewer';
import { useTabState } from './hooks/useTabState';
import { useWsEval, type EvalResult } from './hooks/useWsEval';
import { ResizableSplit } from './components/ResizableSplit';
import { YapCADViewer, type GeometryEntity, type Measurement, type BoundingBoxInfo, type MeasureUnit } from './viewer';
import { loadYapCADGeometry, type LoadedSketch } from './yapcad-loader';
import type { PackageEntry } from './types/package';

// Helper function to set up different view angles
const setupViewCamera = (viewer: YapCADViewer, viewType: string) => {
  const camera = viewer.getCamera();
  const distance = 500;
  
  switch (viewType) {
    case 'front':
      camera.position.set(0, 0, distance);
      camera.up.set(0, 1, 0);
      camera.lookAt(0, 0, 0);
      break;
    case 'top':
      camera.position.set(0, distance, 0);
      camera.up.set(0, 0, -1);
      camera.lookAt(0, 0, 0);
      break;
    case 'right':
      camera.position.set(distance, 0, 0);
      camera.up.set(0, 1, 0);
      camera.lookAt(0, 0, 0);
      break;
    default:
      break;
  }
  
  (camera as THREE.PerspectiveCamera | THREE.OrthographicCamera).updateProjectionMatrix();
};

/**
 * Wrapper that measures its container and renders SketchViewer at the right size.
 * Also handles interactive control-point editing → parameter update → re-eval.
 */
function SketchViewerContainer({
  sketchEntities,
  activeTab,
  updateParam,
  handleCommandEvaluate,
}: {
  sketchEntities: SketchEntity[];
  activeTab: ReturnType<typeof useTabState>['activeTab'];
  updateParam: ReturnType<typeof useTabState>['updateParam'];
  handleCommandEvaluate: (command: string, parameters: Record<string, unknown>, opts?: { silent?: boolean }) => Promise<void>;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState({ width: 600, height: 400 });

  useEffect(() => {
    if (!containerRef.current) return;
    const ro = new ResizeObserver((entries) => {
      for (const entry of entries) {
        const { width, height } = entry.contentRect;
        if (width > 0 && height > 0) setSize({ width, height });
      }
    });
    ro.observe(containerRef.current);
    // Initial read
    const { clientWidth, clientHeight } = containerRef.current;
    if (clientWidth > 0 && clientHeight > 0) setSize({ width: clientWidth, height: clientHeight });
    return () => ro.disconnect();
  }, []);

  // Handle control point drag → map each point back to its point2d param by index,
  // then re-eval the active command AND any sibling commands that share the same params.
  const handleControlPointsChanged = useCallback((newPoints: number[][]) => {
    if (!activeTab) return;
    const cmd = activeTab.selectedCommand;
    if (!cmd) return;

    const cmdDef = activeTab.commands.find(c => c.name === cmd);
    if (!cmdDef) return;

    // Find all point2d params in declaration order
    const point2dParams = cmdDef.params.filter(
      p => p.type === 'point2d' || p.ui_hint?.widget === 'point2d'
    );

    if (point2dParams.length === 0) return;

    // Map each dragged control point back to its DSL param by index
    const updatedParams: Record<string, unknown> = { ...(activeTab.paramValues[cmd] ?? {}) };
    newPoints.forEach((pt, i) => {
      if (i < point2dParams.length) {
        const paramName = point2dParams[i].name;
        const homogeneous = [pt[0], pt[1], pt[2] ?? 0, pt[3] ?? 1];
        updatedParams[paramName] = homogeneous;
        updateParam(activeTab.id, cmd, paramName, homogeneous);
      }
    });

    // Re-eval the active command (updates 2D sketch view)
    handleCommandEvaluate(cmd, updatedParams);

    // Also re-eval any sibling commands that share these param names
    // (e.g. lathe_solid shares p0..p4 with spline_profile → update 3D view)
    const sharedParamNames = new Set(point2dParams.map(p => p.name));
    for (const sibling of activeTab.commands) {
      if (sibling.name === cmd) continue;
      const siblingSharesParams = sibling.params.some(p => sharedParamNames.has(p.name));
      if (siblingSharesParams) {
        const siblingParams: Record<string, unknown> = { ...(activeTab.paramValues[sibling.name] ?? {}) };
        // Inject updated point values
        newPoints.forEach((pt, i) => {
          if (i < point2dParams.length) {
            const paramName = point2dParams[i].name;
            if (sibling.params.some(p => p.name === paramName)) {
              siblingParams[paramName] = [pt[0], pt[1], pt[2] ?? 0, pt[3] ?? 1];
            }
          }
        });
        handleCommandEvaluate(sibling.name, siblingParams, { silent: true });
      }
    }
  }, [activeTab, updateParam, handleCommandEvaluate]);

  // Render the first sketch entity (multi-sketch support can come later)
  const entity = sketchEntities[0];

  // Extract current point2d param values to pass as control point overrides
  const liveControlPoints = (() => {
    if (!activeTab) return undefined;
    const cmd = activeTab.selectedCommand;
    if (!cmd) return undefined;
    const cmdDef = activeTab.commands.find(c => c.name === cmd);
    if (!cmdDef) return undefined;
    const point2dParams = cmdDef.params.filter(
      p => p.type === 'point2d' || (p.ui_hint as Record<string, unknown>)?.widget === 'point2d'
    );
    if (point2dParams.length === 0) return undefined;
    const paramVals = activeTab.paramValues[cmd] ?? {};
    return point2dParams.map(p => {
      const v = paramVals[p.name];
      if (Array.isArray(v)) return v as number[];
      // Fall back to default from eval primitives
      return null;
    }).filter(Boolean) as number[][];
  })();

  return (
    <div
      ref={containerRef}
      style={{ position: 'absolute', inset: 0, overflow: 'hidden' }}
    >
      {entity && (
        <SketchViewer
          sketch={entity}
          controlPoints={liveControlPoints && liveControlPoints.length > 0 ? liveControlPoints : undefined}
          width={size.width}
          height={size.height}
          interactive={true}
          onControlPointsChanged={handleControlPointsChanged}
        />
      )}
    </div>
  );
}

export function App() {
  const [selectedPackage, setSelectedPackage] = useState<PackageEntry | null>(null);
  // Per-mode render settings — each view mode remembers its own lighting + render style.
  // Button state always reflects the active mode; switching modes restores that mode's settings.
  const [lightingPresets, setLightingPresets] = useState<Record<'single' | '4-up', string>>({
    single: 'studio',
    '4-up': 'studio',
  });
  const [renderModes, setRenderModes] = useState<Record<'single' | '4-up', string>>({
    single: 'solid',
    '4-up': 'solid',
  });
  const [leftColumnWidth, setLeftColumnWidth] = useState(() => {
    const saved = localStorage.getItem('yapcad-left-col-width');
    return saved ? parseInt(saved, 10) : 320;
  });
  const [rightColumnWidth, setRightColumnWidth] = useState(() => {
    const saved = localStorage.getItem('yapcad-right-col-width');
    return saved ? parseInt(saved, 10) : 280;
  });
  const [layers, setLayers] = useState<string[]>([]);
  const [layerVisibility, setLayerVisibility] = useState<Record<string, boolean>>({});
  const [viewMode, setViewMode] = useState<'single' | '4-up'>('single');
  // Derived — always reflect the active mode's settings in the toolbar buttons
  const currentLightingPreset = lightingPresets[viewMode];
  const currentRenderMode = renderModes[viewMode];

  // Measurement state
  const [measureMode, setMeasureMode] = useState(false);
  const [measureUnit, setMeasureUnitState] = useState<MeasureUnit>('mm');
  const [measurements, setMeasurements] = useState<Measurement[]>([]);
  const [bboxInfo, setBboxInfo] = useState<BoundingBoxInfo | null>(null);
  
  // Sketch / 2D view state
  const [sketchEntities, setSketchEntities] = useState<SketchEntity[]>([]);
  const [viewerPane, setViewerPane] = useState<'3d' | '2d'>('3d');

  // DSL Editor state
  const [backendUrl, setBackendUrl] = useState(() => {
    const saved = localStorage.getItem('yapcad-backend-url');
    return saved || 'http://localhost:8000';
  });
  const [isEvaluating, setIsEvaluating] = useState(false);
  const [evaluationError, setEvaluationError] = useState<string | null>(null);
  const [isConnected, setIsConnected] = useState(false);
  const [editorChatSplit, setEditorChatSplit] = useState(() => {
    const saved = localStorage.getItem('yapcad-editor-chat-split');
    return saved ? parseFloat(saved) : 0.6;
  });

  // ── Tab state ──────────────────────────────────────────────────────────────
  const {
    tabs, activeTab, activeTabId,
    addTab, closeTab, switchTab, renameTab,
    updateSource, parseCommands, selectCommand, updateParam,
    loadExternalSource,
  } = useTabState({ backendUrl });

  // Derived dslSource from active tab (for chat context and other consumers)
  const dslSource = activeTab?.source ?? '';

  // ── WebSocket eval hook ────────────────────────────────────────────────────
  const wsEval = useWsEval({ backendUrl });

  // ── Eval result state (for parameter panel + chat feedback) ────────────────
  const [lastEvalResult, setLastEvalResult] = useState<EvalResult | null>(null);
  const lastEvalContextRef = useRef<{ command: string; result: EvalResult } | null>(null);

  // FileBar: external source to push into editor
  const [externalSource, setExternalSource] = useState<string | null>(null);

  // FileBar: unsaved-changes gate (prevents ?file= load stomping dirty editor)
  const [hasUnsavedChanges, setHasUnsavedChanges] = useState(false);

  // WBS file requested via ?file= query param — handed to FileBar for loading
  const [pendingWbsFile, setPendingWbsFile] = useState<string | null>(null);

  // On mount: check for ?file=wbs://... query param
  useEffect(() => {
    const params = new URLSearchParams(window.location.search);
    const file = params.get('file');
    if (file && file.startsWith('wbs://')) {
      if (hasUnsavedChanges) {
        const ok = window.confirm(
          `You have unsaved changes. Discard them and open "${file}"?`
        );
        if (!ok) return;
      }
      setPendingWbsFile(file);
      // Clean the URL so a page refresh doesn't re-trigger the load
      const url = new URL(window.location.href);
      url.searchParams.delete('file');
      window.history.replaceState({}, '', url.toString());
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // BroadcastChannel listener — receives 'open-file' messages from the WBS dashboard
  // when the user clicks '⬡ DSL' on a part row while this tab is already open.
  useEffect(() => {
    let ch: BroadcastChannel | null = null;
    try {
      ch = new BroadcastChannel('yapcad-wbs');
      ch.onmessage = (ev) => {
        if (ev.data?.type === 'open-file' && typeof ev.data.path === 'string') {
          const path = ev.data.path as string;
          if (!path.startsWith('wbs://')) return;
          if (hasUnsavedChanges) {
            const ok = window.confirm(
              `You have unsaved changes. Discard them and open "${path}"?`
            );
            if (!ok) return;
          }
          setPendingWbsFile(path);
        }
      };
    } catch (e) {
      console.warn('BroadcastChannel not available:', e);
    }
    return () => { ch?.close(); };
  }, [hasUnsavedChanges]);
  
  const viewerContainerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<YapCADViewer | null>(null);
  // Persists last-loaded geometry across view mode switches (single ↔ 4-up).
  // Viewer disposal wipes the Three.js scene; this ref lets us restore it.
  const lastEntitiesRef = useRef<GeometryEntity[] | null>(null);
  // Refs for per-mode render settings — read by initializeMultiView's RAF closure
  // (which has [] deps and can't capture updating state directly).
  const lightingPresetsRef = useRef<Record<'single' | '4-up', string>>({ single: 'studio', '4-up': 'studio' });
  const renderModesRef = useRef<Record<'single' | '4-up', string>>({ single: 'solid', '4-up': 'solid' });
  const multiViewRef = useRef<{
    perspective: YapCADViewer;
    front: YapCADViewer;
    top: YapCADViewer;
    right: YapCADViewer;
  } | null>(null);

  // DSL Editor handlers
  const handleBackendUrlChange = useCallback((url: string) => {
    // Auto-prepend http:// if no protocol specified
    const normalizedUrl = url.match(/^https?:\/\//) ? url : `http://${url}`;
    setBackendUrl(normalizedUrl);
    localStorage.setItem('yapcad-backend-url', normalizedUrl);
    setIsConnected(false);
  }, []);

  const handleTestConnection = useCallback(async () => {
    try {
      const response = await fetch(`${backendUrl}/health`, { 
        method: 'GET',
        signal: typeof AbortSignal.timeout === "function" ? AbortSignal.timeout(5000) : undefined
      });
      setIsConnected(response.ok);
    } catch (error) {
      console.warn('Connection test failed:', error);
      setIsConnected(false);
    }
  }, [backendUrl]);

  // Load entities into whichever viewer(s) are currently active + persist for mode switches.
  // Uses refs directly so it's stable across renders (safe to call from RAF callbacks).
  const applyToActiveViewers = useCallback((entities: GeometryEntity[]) => {
    lastEntitiesRef.current = entities;

    const updateLayers = (viewer: YapCADViewer) => {
      const newLayers = viewer.getLayers();
      setLayers(newLayers);
      const vis: Record<string, boolean> = {};
      for (const l of newLayers) vis[l] = true;
      setLayerVisibility(vis);
    };

    if (viewerRef.current) {
      viewerRef.current.loadGeometry(entities);
      viewerRef.current.fitToGeometry();
      updateLayers(viewerRef.current);
      setBboxInfo(viewerRef.current.getBoundingBoxInfo());
    } else if (multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(v => {
        v.loadGeometry(entities);
        v.fitToGeometry();
      });
      updateLayers(multiViewRef.current.perspective);
    }
  // setLayers and setLayerVisibility are stable; refs are always current → no deps needed
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Extract first command name for full-editor evaluate
  const extractCommandName = useCallback((source: string): string => {
    const match = source.match(/command\s+(\w+)/);
    return match ? match[1] : 'main';
  }, []);

  const handleDSLEvaluate = useCallback(async (source: string) => {
    return handleCommandEvaluate(extractCommandName(source), {});
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [backendUrl, viewMode]);

  const handleCommandEvaluate = useCallback(async (command: string, parameters: Record<string, unknown>, opts?: { silent?: boolean }) => {
    const silent = opts?.silent ?? false;
    if (!silent) setIsEvaluating(true);
    if (!silent) setEvaluationError(null);
    
    try {
      // Use WebSocket eval (with REST fallback)
      const evalResult = await wsEval.eval(dslSource, command, parameters);

      // Record latency for tier detection
      wsEval.recordLatency(command, evalResult.elapsed_ms);

      // Store eval result for parameter panel + chat feedback
      setLastEvalResult(evalResult);
      lastEvalContextRef.current = { command, result: evalResult };

      if (!evalResult.success) {
        throw new Error(evalResult.error_message || 'DSL evaluation failed');
      }

      if (!evalResult.geometry) {
        throw new Error('No geometry returned');
      }
      
      // Cast to YapCADDocument — the backend returns this structure
      const result = loadYapCADGeometry(evalResult.geometry as unknown as Parameters<typeof loadYapCADGeometry>[0]);
      const entities: GeometryEntity[] = [];
      
      for (const solid of result.solids) {
        for (const surface of solid.surfaces) {
          const positionAttr = surface.geometry.getAttribute('position');
          const normalAttr = surface.geometry.getAttribute('normal');
          const indexAttr = surface.geometry.getIndex();
          
          entities.push({
            id: `${solid.id}-${surface.id}`,
            name: `${solid.id}/${surface.id}`,
            vertices: new Float32Array(positionAttr.array),
            normals: normalAttr ? new Float32Array(normalAttr.array) : undefined,
            indices: indexAttr ? new Uint32Array(indexAttr.array) : undefined,
            material: surface.materialId || 'default',
            layer: surface.layer,
          });
        }
      }
      
      for (const surface of result.surfaces) {
        const positionAttr = surface.geometry.getAttribute('position');
        const normalAttr = surface.geometry.getAttribute('normal');
        const indexAttr = surface.geometry.getIndex();
        
        entities.push({
          id: surface.id,
          name: surface.id,
          vertices: new Float32Array(positionAttr.array),
          normals: normalAttr ? new Float32Array(normalAttr.array) : undefined,
          indices: indexAttr ? new Uint32Array(indexAttr.array) : undefined,
          material: surface.materialId || 'default',
          layer: surface.layer,
        });
      }

      applyToActiveViewers(entities);

      // Collect sketch entities for the 2D viewer
      const newSketches: SketchEntity[] = result.sketches.map((s: LoadedSketch) => ({
        id: s.id,
        type: 'sketch' as const,
        name: s.name,
        boundingBox: s.boundingBox,
        polylines: s.polylines,
        primitives: s.primitives,
        metadata: s.metadata,
      }));
      if (!silent) {
        setSketchEntities(newSketches);
        // Auto-switch pane: if only sketches (no 3D geometry), show 2D; if only 3D, show 3D
        if (newSketches.length > 0 && entities.length === 0) {
          setViewerPane('2d');
        } else if (entities.length > 0 && newSketches.length === 0) {
          setViewerPane('3d');
        }
        setSelectedPackage(null);
        setIsConnected(true);
      }
      
    } catch (error) {
      if (error instanceof Error && error.message === 'Cancelled') return;
      if (!silent) {
        console.error('DSL evaluation failed:', error);
        const errorMessage = error instanceof Error ? error.message : String(error);
        setEvaluationError(errorMessage);
        setIsConnected(false);
      }
    } finally {
      if (!silent) setIsEvaluating(false);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [backendUrl, dslSource, applyToActiveViewers, wsEval]);

  // ── Tab source change → re-parse commands ──────────────────────────────────
  const handleTabSourceChange = useCallback((source: string) => {
    if (!activeTab) return;
    updateSource(activeTab.id, source);
    // Debounced re-parse is handled by a useEffect below
  }, [activeTab, updateSource]);

  // Debounced command parsing when active tab source changes
  const parseTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  useEffect(() => {
    if (!activeTab) return;
    if (parseTimerRef.current) clearTimeout(parseTimerRef.current);
    parseTimerRef.current = setTimeout(() => {
      parseCommands(activeTab.id, activeTab.source);
    }, 800);
    return () => { if (parseTimerRef.current) clearTimeout(parseTimerRef.current); };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [activeTab?.id, activeTab?.source, backendUrl]);

  // ── Parameter panel eval trigger ───────────────────────────────────────────
  const handleParamEval = useCallback((command: string, parameters: Record<string, unknown>) => {
    handleCommandEvaluate(command, parameters);
  }, [handleCommandEvaluate]);

  const handleEditorChatSplitChange = useCallback((split: number) => {
    setEditorChatSplit(split);
    localStorage.setItem('yapcad-editor-chat-split', String(split));
  }, []);

  // FileBar handlers
  const handleFileLoad = useCallback((source: string, _filename: string) => {
    setExternalSource(source);
    // Also load into active tab
    loadExternalSource(source);
  }, [loadExternalSource]);

  const handleExternalSourceConsumed = useCallback(() => {
    setExternalSource(null);
  }, []);

  // Initialize 4-up view
  const initializeMultiView = useCallback(() => {
    if (!viewerContainerRef.current) return;
    
    viewerContainerRef.current.innerHTML = '';
    
    const quadrants = [
      { name: 'perspective', label: 'Perspective' },
      { name: 'front', label: 'Front' },
      { name: 'top', label: 'Top' },
      { name: 'right', label: 'Right' },
    ] as const;
    
    const viewers: any = {};
    
    quadrants.forEach((quad) => {
      const quadContainer = document.createElement('div');
      quadContainer.style.cssText = `
        position: relative;
        border: 2px solid #666;
        background: #0f0f1a;
        overflow: hidden;
        min-height: 0;
        min-width: 0;
        width: 100%;
        height: 100%;
        display: flex;
        flex-direction: column;
      `;
      
      const label = document.createElement('div');
      label.textContent = quad.label;
      label.style.cssText = `
        position: absolute;
        top: 4px;
        left: 8px;
        font-size: 12px;
        color: #888;
        z-index: 10;
        pointer-events: none;
      `;
      quadContainer.appendChild(label);
      
      const viewerDiv = document.createElement('div');
      viewerDiv.style.cssText = 'width: 100%; height: 100%; flex: 1; position: relative;';
      quadContainer.appendChild(viewerDiv);
      
      viewerContainerRef.current!.appendChild(quadContainer);
      
      const viewer = new YapCADViewer({
        container: viewerDiv,
        antialias: true,
        orthographic: quad.name !== 'perspective',
      });
      
      if (quad.name !== 'perspective') {
        setupViewCamera(viewer, quad.name);
        const controls = (viewer as any).controls;
        if (controls && controls.enableRotate !== undefined) {
          controls.enableRotate = false;
        }
      }
      
      viewer.start();
      // Wire up measurement callback — merges measurements from all 4 viewers
      viewer.onMeasureChange(() => {
        if (!multiViewRef.current) return;
        const all: Measurement[] = Object.values(multiViewRef.current).flatMap(v => v.getMeasurements());
        setMeasurements([...all]);
      });
      viewers[quad.name] = viewer;
    });
    
    multiViewRef.current = viewers;
    
    // Double-RAF ensures browser layout is fully resolved before reading dimensions.
    requestAnimationFrame(() => {
      requestAnimationFrame(() => {
        if (!multiViewRef.current) return;
        Object.values(multiViewRef.current).forEach(v => v.handleResize());
        window.dispatchEvent(new Event('resize'));
        // Restore this mode's render settings (read from refs — closure has [] deps)
        Object.values(multiViewRef.current).forEach(v => {
          v.setLightingPreset(lightingPresetsRef.current['4-up']);
          v.setRenderMode(renderModesRef.current['4-up']);
        });
        // Restore last geometry
        if (lastEntitiesRef.current) {
          Object.values(multiViewRef.current).forEach(v => {
            v.loadGeometry(lastEntitiesRef.current!);
            v.fitToGeometry();
          });
        }
      });
    });

    // Belt-and-suspenders: retry resize after 200ms for slow layout environments
    setTimeout(() => {
      if (multiViewRef.current) {
        Object.values(multiViewRef.current).forEach(v => v.handleResize());
      }
    }, 200);
  }, []);
  
  // Initialize viewer when container is ready
  useEffect(() => {
    if (!viewerContainerRef.current) return;

    if (viewMode === 'single') {
      if (multiViewRef.current) {
        Object.values(multiViewRef.current).forEach(viewer => viewer.dispose());
        multiViewRef.current = null;
      }
      
      if (!viewerRef.current) {
        // Clear any leftover imperative DOM nodes from 4-up mode (quadrant divs)
        // before appending the single viewer's canvas.
        viewerContainerRef.current.innerHTML = '';
        viewerRef.current = new YapCADViewer({
          container: viewerContainerRef.current,
          antialias: true,
        });
        viewerRef.current.start();
        // Restore this mode's render settings
        viewerRef.current.setLightingPreset(lightingPresetsRef.current['single']);
        viewerRef.current.setRenderMode(renderModesRef.current['single']);
        // Wire up measurement callback
        viewerRef.current.onMeasureChange((ms) => {
          setMeasurements([...ms]);
        });
        // Restore last geometry
        if (lastEntitiesRef.current) {
          viewerRef.current.loadGeometry(lastEntitiesRef.current);
          viewerRef.current.fitToGeometry();
        }
      }
    } else {
      if (viewerRef.current) {
        viewerRef.current.dispose();
        viewerRef.current = null;
      }
      
      if (!multiViewRef.current) {
        initializeMultiView();
      }
    }
    
    return () => {
      if (viewerRef.current) {
        viewerRef.current.dispose();
        viewerRef.current = null;
      }
      if (multiViewRef.current) {
        Object.values(multiViewRef.current).forEach(viewer => viewer.dispose());
        multiViewRef.current = null;
      }
    };
  }, [viewMode]);
  
  // Load geometry when package changes
  useEffect(() => {
    if (!selectedPackage?.geometry) return;
    
    const entities: GeometryEntity[] = selectedPackage.geometry.entities.map(e => ({
      id: e.id,
      name: e.name,
      vertices: e.vertices,
      normals: e.normals,
      indices: e.indices,
      material: e.material,
      layer: e.layer,
    }));
    
    applyToActiveViewers(entities);
  // viewMode intentionally excluded: viewer restoration on mode switch is handled
  // inside initializeMultiView (RAF) and the viewer-init useEffect below.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedPackage, applyToActiveViewers]);
  
  const handleSelect = useCallback((pkg: PackageEntry | null) => {
    setSelectedPackage(pkg);
  }, []);

  const handleLightingChange = useCallback((preset: string) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setLightingPreset(preset);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(v => v.setLightingPreset(preset));
    }
    setLightingPresets(prev => {
      const next = { ...prev, [viewMode]: preset };
      lightingPresetsRef.current = next;
      return next;
    });
  }, [viewMode]);

  const handleRenderModeChange = useCallback((mode: string) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setRenderMode(mode);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(v => v.setRenderMode(mode));
    }
    setRenderModes(prev => {
      const next = { ...prev, [viewMode]: mode };
      renderModesRef.current = next;
      return next;
    });
  }, [viewMode]);

  const handleClipPlaneChange = useCallback((
    axis: 'x' | 'y' | 'z',
    enabled: boolean,
    position: number,
    invert: boolean
  ) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setClipPlane(axis, enabled, position, invert);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => viewer.setClipPlane(axis, enabled, position, invert));
    }
  }, [viewMode]);

  const handleResetClipPlanes = useCallback(() => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.resetClipPlanes();
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => viewer.resetClipPlanes());
    }
  }, [viewMode]);

  const getClippingState = useCallback(() => {
    if (viewMode === 'single' && viewerRef.current) {
      return viewerRef.current.getClippingState();
    } else if (viewMode === '4-up' && multiViewRef.current) {
      return multiViewRef.current.perspective.getClippingState();
    }
    return [];
  }, [viewMode]);

  const handleLeftResize = useCallback((delta: number) => {
    setLeftColumnWidth(w => {
      const newWidth = Math.max(250, Math.min(500, w + delta));
      localStorage.setItem('yapcad-left-col-width', String(newWidth));
      return newWidth;
    });
  }, []);

  const handleRightResize = useCallback((delta: number) => {
    setRightColumnWidth(w => {
      const newWidth = Math.max(200, Math.min(450, w - delta));
      localStorage.setItem('yapcad-right-col-width', String(newWidth));
      return newWidth;
    });
  }, []);

  const handleLayerToggle = useCallback((layer: string, visible: boolean) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setLayerVisibility(layer, visible);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => viewer.setLayerVisibility(layer, visible));
    }
    setLayerVisibility(prev => ({ ...prev, [layer]: visible }));
  }, [viewMode]);

  const handleShowAllLayers = useCallback(() => {
    const vis: Record<string, boolean> = {};
    for (const l of layers) vis[l] = true;
    if (viewMode === 'single' && viewerRef.current) {
      for (const l of layers) viewerRef.current.setLayerVisibility(l, true);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => {
        for (const l of layers) viewer.setLayerVisibility(l, true);
      });
    }
    setLayerVisibility(vis);
  }, [layers, viewMode]);

  const handleHideAllLayers = useCallback(() => {
    const vis: Record<string, boolean> = {};
    for (const l of layers) vis[l] = false;
    if (viewMode === 'single' && viewerRef.current) {
      for (const l of layers) viewerRef.current.setLayerVisibility(l, false);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => {
        for (const l of layers) viewer.setLayerVisibility(l, false);
      });
    }
    setLayerVisibility(vis);
  }, [layers, viewMode]);

  // Measurement handlers — single view: one viewer; 4-up: all four viewers
  const _allViewers = useCallback((): YapCADViewer[] => {
    if (viewerRef.current) return [viewerRef.current];
    if (multiViewRef.current) return Object.values(multiViewRef.current);
    return [];
  }, []);

  const handleToggleMeasure = useCallback(() => {
    const viewers = _allViewers();
    if (viewers.length === 0) return;
    const next = !measureMode;
    viewers.forEach(v => v.setMeasureMode(next));
    setMeasureMode(next);
  }, [measureMode, _allViewers]);

  const handleDeleteMeasurement = useCallback((id: string) => {
    _allViewers().forEach(v => v.deleteMeasurement(id));
  }, [_allViewers]);

  const handleClearMeasurements = useCallback(() => {
    _allViewers().forEach(v => v.clearMeasurements());
    setMeasurements([]);
  }, [_allViewers]);

  const handleUnitToggle = useCallback(() => {
    const next: MeasureUnit = measureUnit === 'mm' ? 'in' : 'mm';
    _allViewers().forEach(v => v.setMeasureUnit(next));
    setMeasureUnitState(next);
    // Refresh measurements list (labels rebuilt in viewer, update React state too)
    const all = _allViewers().flatMap(v => v.getMeasurements());
    setMeasurements([...all]);
  }, [measureUnit, _allViewers]);

  // Format a raw mm value for display in the measurement panel
  const fmtDist = useCallback((mm: number): string => {
    if (measureUnit === 'in') return `${(mm / 25.4).toFixed(3)}"`;
    return `${mm.toFixed(2)} mm`;
  }, [measureUnit]);

  const fmtBbox = useCallback((mm: number): string => {
    if (measureUnit === 'in') return `${(mm / 25.4).toFixed(3)}"`;
    return `${mm.toFixed(1)} mm`;
  }, [measureUnit]);
  
  return (
    <div style={{
      display: 'flex',
      height: '100vh',
      width: '100vw',
      overflow: 'hidden',
      fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
    }}>
      {/* LEFT COLUMN: FileBar + SessionBar + Editor/Chat */}
      <div style={{
        width: `${leftColumnWidth}px`,
        flexShrink: 0,
        display: 'flex',
        flexDirection: 'column',
        minWidth: '250px',
        maxWidth: '500px',
        borderRight: '1px solid #333',
        backgroundColor: '#1a1a2e',
      }}>
        {/* File Management Bar */}
        <FileBar
          currentSource={dslSource}
          onFileLoad={handleFileLoad}
          onUnsavedChange={setHasUnsavedChanges}
          pendingWbsFile={pendingWbsFile}
          onPendingWbsFileConsumed={() => setPendingWbsFile(null)}
        />
        
        {/* Session Bar */}
        <SessionBar
          backendUrl={backendUrl}
          onBackendUrlChange={handleBackendUrlChange}
          isConnected={isConnected}
          onTestConnection={handleTestConnection}
        />
        
        {/* Resizable Editor/Chat Split */}
        <ResizableSplit
          direction="vertical"
          initialSplit={editorChatSplit}
          minSize={200}
          onSplitChange={handleEditorChatSplitChange}
          style={{ flex: 1 }}
        >
          {/* Tab Bar + DSL Editor + Parameter Panel */}
          <div style={{ display: 'flex', flexDirection: 'column', height: '100%', overflow: 'hidden' }}>
            {/* Tab Bar */}
            <TabBar
              tabs={tabs}
              activeTabId={activeTabId}
              onSwitch={switchTab}
              onClose={closeTab}
              onAdd={addTab}
              onRename={renameTab}
            />
            
            {/* DSL Editor */}
            <div style={{ flex: 1, overflow: 'hidden' }}>
              <DSLEditor
                onEvaluate={handleDSLEvaluate}
                onSourceChange={handleTabSourceChange}
                externalSource={externalSource}
                onExternalSourceConsumed={handleExternalSourceConsumed}
                isEvaluating={isEvaluating}
                evaluationError={evaluationError}
                key={activeTabId}
                initialSource={activeTab?.source}
              />
            </div>
            
            {/* Parameter Panel */}
            <ParameterPanel
              commands={activeTab?.commands ?? []}
              selectedCommand={activeTab?.selectedCommand ?? ''}
              paramValues={activeTab?.paramValues ?? {}}
              onSelectCommand={(cmd) => activeTab && selectCommand(activeTab.id, cmd)}
              onParamChange={(cmd, param, val) => activeTab && updateParam(activeTab.id, cmd, param, val)}
              onEval={handleParamEval}
              isEvaluating={isEvaluating}
              evalResult={lastEvalResult}
              latencyTier={wsEval.getLatencyTier(activeTab?.selectedCommand ?? '')}
              parseError={activeTab?.parseError ?? ''}
              isParsing={activeTab?.isParsing ?? false}
            />
          </div>
          
          {/* Chat Panel — with eval context ref for feedback */}
          <ChatPanel
            dslSource={dslSource}
            onDSLUpdate={(source) => {
              setExternalSource(source);
              loadExternalSource(source);
            }}
            evalContextRef={lastEvalContextRef}
          />
        </ResizableSplit>
      </div>
      
      <ResizeHandle direction="vertical" onResize={handleLeftResize} />
      
      {/* CENTER COLUMN: Toolbar + 3D Viewer */}
      <div style={{
        flex: 1,
        minWidth: 0,
        display: 'flex',
        flexDirection: 'column',
        backgroundColor: '#0f0f1a',
        overflow: 'hidden',
      }}>
        {/* Toolbar */}
        <div style={{
          minHeight: '48px',
          display: 'flex',
          alignItems: 'center',
          flexWrap: 'wrap',
          padding: '8px 16px',
          borderBottom: '1px solid #333',
          backgroundColor: '#1a1a2e',
          gap: '8px 16px',
        }}>
          <span style={{ flex: 1, fontSize: '14px', color: '#eee', fontWeight: 500 }}>
            {selectedPackage ? selectedPackage.name : 'yapCAD Viewer'}
          </span>
          
          <LightingSelector
            currentPreset={currentLightingPreset}
            onPresetChange={handleLightingChange}
          />
          
          <RenderModeSelector
            currentMode={currentRenderMode}
            onModeChange={handleRenderModeChange}
          />

          <button
            style={{
              padding: '6px 12px',
              fontSize: '13px',
              backgroundColor: viewMode === '4-up' ? '#3b82f6' : '#333',
              color: '#eee',
              border: 'none',
              borderRadius: '4px',
              cursor: 'pointer',
            }}
            onClick={() => setViewMode(viewMode === 'single' ? '4-up' : 'single')}
          >
            {viewMode === 'single' ? '4-Up' : 'Single'} View
          </button>

          {/* 2D / 3D toggle — only show when sketch entities are available */}
          {sketchEntities.length > 0 && (
            <button
              style={{
                padding: '6px 12px',
                fontSize: '13px',
                backgroundColor: viewerPane === '2d' ? '#8b5cf6' : '#333',
                color: '#eee',
                border: viewerPane === '2d' ? '1px solid #a78bfa' : 'none',
                borderRadius: '4px',
                cursor: 'pointer',
                fontWeight: viewerPane === '2d' ? 600 : 400,
              }}
              onClick={() => setViewerPane(viewerPane === '2d' ? '3d' : '2d')}
              title={viewerPane === '2d' ? 'Switch to 3D viewer' : 'Switch to 2D sketch viewer'}
            >
              {viewerPane === '2d' ? '✏️ 2D' : '📐 2D'}
            </button>
          )}

          <button
            title={measureMode ? 'Exit Measure mode (click 2 points to measure distance)' : 'Enter Measure mode (click 2 points to measure distance)'}
            style={{
              padding: '6px 12px',
              fontSize: '13px',
              backgroundColor: measureMode ? '#f59e0b' : '#333',
              color: measureMode ? '#1a1a00' : '#eee',
              border: measureMode ? '1px solid #fbbf24' : 'none',
              borderRadius: '4px',
              cursor: 'pointer',
              fontWeight: measureMode ? 600 : 400,
            }}
            onClick={handleToggleMeasure}
          >
            📏 {measureMode ? 'Measuring…' : 'Measure'}
          </button>
        </div>
        
        {/* Clipping Controls */}
        <div style={{ padding: '8px 16px' }}>
          <ClippingControls
            onClipPlaneChange={handleClipPlaneChange}
            onResetClipPlanes={handleResetClipPlanes}
            getClippingState={getClippingState}
          />
        </div>
        
        {/* Viewer area — 3D or 2D depending on pane selection */}
        <div style={{ flex: 1, minWidth: 0, minHeight: 0, position: 'relative', overflow: 'hidden' }}>

          {/* 3D Imperative viewer container — hidden when 2D pane is active */}
          <div 
            ref={viewerContainerRef}
            style={{
              position: 'absolute',
              inset: 0,
              overflow: 'hidden',
              display: viewerPane === '2d' ? 'none' : undefined,
              ...(viewMode === '4-up' && viewerPane !== '2d' ? {
                display: 'grid',
                gridTemplateColumns: '1fr 1fr',
                gridTemplateRows: '1fr 1fr',
                gap: '2px',
              } : {}),
            }}
          />

          {/* 2D Sketch Viewer — shown when 2D pane is active and sketches exist */}
          {viewerPane === '2d' && sketchEntities.length > 0 && (
            <SketchViewerContainer
              sketchEntities={sketchEntities}
              activeTab={activeTab}
              updateParam={updateParam}
              handleCommandEvaluate={handleCommandEvaluate}
            />
          )}

          {/* React-managed overlays — absolutely positioned, never interfere with viewer grid */}
          {viewerPane === '3d' && !selectedPackage && !isEvaluating && sketchEntities.length === 0 && (
            <div style={{
              position: 'absolute',
              inset: 0,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              color: '#666',
              fontSize: '16px',
              pointerEvents: 'none',
              flexDirection: 'column',
              gap: '12px',
              zIndex: 10,
            }}>
              <div>Select a package or write DSL code to view</div>
            </div>
          )}
          
          {isEvaluating && (
            <div style={{
              position: 'absolute',
              inset: 0,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              color: '#4c9f38',
              fontSize: '16px',
              pointerEvents: 'none',
              zIndex: 10,
            }}>
              ⏳ Evaluating DSL...
            </div>
          )}

          {/* Measurement Panel — bottom-right overlay */}
          {(measureMode || measurements.length > 0 || bboxInfo) && (
            <div style={{
              position: 'absolute',
              bottom: '16px',
              right: '16px',
              zIndex: 20,
              backgroundColor: 'rgba(15,15,30,0.92)',
              border: '1px solid #334',
              borderRadius: '8px',
              padding: '12px 14px',
              minWidth: '240px',
              maxWidth: '320px',
              color: '#ccd',
              fontSize: '13px',
              pointerEvents: 'all',
            }}>
              {/* Header row: title + unit toggle + clear */}
              <div style={{ display: 'flex', alignItems: 'center', gap: '6px', marginBottom: '8px' }}>
                <span style={{ fontWeight: 600, color: '#f59e0b', fontSize: '12px', letterSpacing: '0.05em', textTransform: 'uppercase', flex: 1 }}>
                  📏 Measurements
                </span>
                {/* mm / in toggle */}
                <div style={{ display: 'flex', borderRadius: '4px', overflow: 'hidden', border: '1px solid #445', fontSize: '11px' }}>
                  {(['mm', 'in'] as const).map(u => (
                    <button
                      key={u}
                      onClick={handleUnitToggle}
                      style={{
                        padding: '2px 7px',
                        background: measureUnit === u ? '#334' : 'transparent',
                        color: measureUnit === u ? '#aaccff' : '#556',
                        border: 'none',
                        cursor: 'pointer',
                        fontWeight: measureUnit === u ? 600 : 400,
                      }}
                    >
                      {u}
                    </button>
                  ))}
                </div>
                {measurements.length > 0 && (
                  <button
                    onClick={handleClearMeasurements}
                    style={{ background: 'none', border: 'none', color: '#556', cursor: 'pointer', fontSize: '11px', padding: '0' }}
                    title="Clear all measurements"
                  >
                    Clear
                  </button>
                )}
              </div>

              {/* Bounding box info */}
              {bboxInfo && (
                <div style={{ marginBottom: measurements.length > 0 ? '10px' : 0, paddingBottom: measurements.length > 0 ? '10px' : 0, borderBottom: measurements.length > 0 ? '1px solid #223' : 'none' }}>
                  <div style={{ color: '#556', fontSize: '10px', marginBottom: '4px', textTransform: 'uppercase', letterSpacing: '0.06em' }}>Bounding Box</div>
                  <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '4px' }}>
                    {(['width', 'height', 'depth'] as const).map(dim => (
                      <div key={dim} style={{ textAlign: 'center', backgroundColor: '#111828', borderRadius: '4px', padding: '4px 6px' }}>
                        <div style={{ fontSize: '10px', color: '#445', textTransform: 'uppercase', marginBottom: '2px' }}>{dim === 'width' ? 'W' : dim === 'height' ? 'H' : 'D'}</div>
                        <div style={{ fontWeight: 600, color: '#99bbdd', fontFamily: 'monospace', fontSize: '12px' }}>{fmtBbox(bboxInfo[dim])}</div>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {/* Hint when in measure mode with no measurements yet */}
              {measureMode && measurements.length === 0 && (
                <div style={{ color: '#f59e0b', fontSize: '12px', textAlign: 'center', padding: '6px 0', opacity: 0.85 }}>
                  Click a surface point to start
                </div>
              )}

              {/* Measurements list */}
              {measurements.map((m, i) => (
                <div
                  key={m.id}
                  style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    padding: '5px 6px',
                    marginBottom: '3px',
                    backgroundColor: '#111828',
                    borderRadius: '4px',
                    borderLeft: '2px solid #44aaff',
                  }}
                >
                  <span style={{ fontSize: '12px', color: '#aac' }}>
                    <span style={{ color: '#445', marginRight: '6px' }}>#{i + 1}</span>
                    <span style={{ fontWeight: 600, color: '#ddf', fontFamily: 'monospace' }}>{fmtDist(m.distance)}</span>
                  </span>
                  <button
                    onClick={() => handleDeleteMeasurement(m.id)}
                    style={{ background: 'none', border: 'none', color: '#445', cursor: 'pointer', fontSize: '14px', lineHeight: 1, padding: '0 2px' }}
                    title="Delete"
                  >
                    ×
                  </button>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
      
      <ResizeHandle direction="vertical" onResize={handleRightResize} />
      
      {/* RIGHT COLUMN: Package Gallery + Layers */}
      <div style={{
        width: `${rightColumnWidth}px`,
        flexShrink: 0,
        display: 'flex',
        flexDirection: 'column',
        minWidth: '200px',
        maxWidth: '450px',
        borderLeft: '1px solid #333',
        backgroundColor: '#1a1a2e',
        overflow: 'hidden',
      }}>
        <div style={{ flex: 1, minHeight: 0, overflow: 'hidden' }}>
          <PackageSelector onSelect={handleSelect} />
        </div>
        {layers.length > 0 && (
          <div style={{ borderTop: '1px solid #333', maxHeight: '250px', overflow: 'hidden' }}>
            <LayerPanel
              layers={layers}
              visibility={layerVisibility}
              onToggle={handleLayerToggle}
              onShowAll={handleShowAllLayers}
              onHideAll={handleHideAllLayers}
            />
          </div>
        )}
      </div>
    </div>
  );
}

export default App;
