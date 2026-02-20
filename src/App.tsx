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
import { CommandSelector } from './components/CommandSelector';
import { ChatPanel } from './components/ChatPanel';
import { SessionBar } from './components/SessionBar';
import { FileBar } from './components/FileBar';

/** Extract the first command name from DSL source */
function extractCommandName(source: string): string {
  const match = source.match(/command\s+(\w+)/);
  return match ? match[1] : 'main';
}
import { ResizableSplit } from './components/ResizableSplit';
import { YapCADViewer, type GeometryEntity } from './viewer';
import { loadYapCADGeometry } from './yapcad-loader';
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

export function App() {
  const [selectedPackage, setSelectedPackage] = useState<PackageEntry | null>(null);
  const [currentLightingPreset, setCurrentLightingPreset] = useState('studio');
  const [currentRenderMode, setCurrentRenderMode] = useState('solid');
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
  
  // DSL Editor state
  const [backendUrl, setBackendUrl] = useState(() => {
    const saved = localStorage.getItem('yapcad-backend-url');
    return saved || 'http://localhost:8000';
  });
  const [isEvaluating, setIsEvaluating] = useState(false);
  const [evaluationError, setEvaluationError] = useState<string | null>(null);
  const [dslSource, setDslSource] = useState('');
  const [isConnected, setIsConnected] = useState(false);
  const [editorChatSplit, setEditorChatSplit] = useState(() => {
    const saved = localStorage.getItem('yapcad-editor-chat-split');
    return saved ? parseFloat(saved) : 0.6;
  });
  
  // FileBar: external source to push into editor
  const [externalSource, setExternalSource] = useState<string | null>(null);
  
  const viewerContainerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<YapCADViewer | null>(null);
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

  const handleDSLEvaluate = useCallback(async (source: string) => {
    return handleCommandEvaluate(extractCommandName(source), {});
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [backendUrl, viewMode]);

  const handleCommandEvaluate = useCallback(async (command: string, parameters: Record<string, unknown>) => {
    setIsEvaluating(true);
    setEvaluationError(null);
    
    try {
      const response = await fetch(`${backendUrl}/dsl/eval`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          source: dslSource,
          command,
          parameters,
          format: 'json'
        }),
        signal: typeof AbortSignal.timeout === "function" ? AbortSignal.timeout(300000) : undefined
      });

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Server error (${response.status}): ${errorText}`);
      }

      const responseJson = await response.json();
      
      if (!responseJson.success) {
        throw new Error(responseJson.error_message || 'DSL evaluation failed');
      }
      
      const result = loadYapCADGeometry(responseJson.geometry);
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

      if (viewMode === 'single' && viewerRef.current) {
        viewerRef.current.loadGeometry(entities);
        viewerRef.current.fitToGeometry();
        const newLayers = viewerRef.current.getLayers();
        setLayers(newLayers);
        const vis: Record<string, boolean> = {};
        for (const l of newLayers) vis[l] = true;
        setLayerVisibility(vis);
      } else if (viewMode === '4-up' && multiViewRef.current) {
        Object.values(multiViewRef.current).forEach(viewer => {
          viewer.loadGeometry(entities);
          viewer.fitToGeometry();
        });
        const newLayers = multiViewRef.current.perspective.getLayers();
        setLayers(newLayers);
        const vis: Record<string, boolean> = {};
        for (const l of newLayers) vis[l] = true;
        setLayerVisibility(vis);
      }

      setSelectedPackage(null);
      setIsConnected(true);
      
    } catch (error) {
      console.error('DSL evaluation failed:', error);
      const errorMessage = error instanceof Error ? error.message : String(error);
      setEvaluationError(errorMessage);
      setIsConnected(false);
    } finally {
      setIsEvaluating(false);
    }
  }, [backendUrl, dslSource, viewMode]);

  const handleEditorChatSplitChange = useCallback((split: number) => {
    setEditorChatSplit(split);
    localStorage.setItem('yapcad-editor-chat-split', String(split));
  }, []);

  // FileBar handlers
  const handleFileLoad = useCallback((source: string, _filename: string) => {
    setExternalSource(source);
  }, []);

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
      viewers[quad.name] = viewer;
    });
    
    multiViewRef.current = viewers;
    
    requestAnimationFrame(() => {
      if (multiViewRef.current) {
        Object.values(multiViewRef.current).forEach(viewer => {
          viewer.handleResize();
        });
        window.dispatchEvent(new Event('resize'));
      }
    });
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
        viewerRef.current = new YapCADViewer({
          container: viewerContainerRef.current,
          antialias: true,
        });
        viewerRef.current.start();
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
    
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.loadGeometry(entities);
      viewerRef.current.fitToGeometry();
      const newLayers = viewerRef.current.getLayers();
      setLayers(newLayers);
      const vis: Record<string, boolean> = {};
      for (const l of newLayers) vis[l] = true;
      setLayerVisibility(vis);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => {
        viewer.loadGeometry(entities);
        viewer.fitToGeometry();
      });
      const newLayers = multiViewRef.current.perspective.getLayers();
      setLayers(newLayers);
      const vis: Record<string, boolean> = {};
      for (const l of newLayers) vis[l] = true;
      setLayerVisibility(vis);
    }
  }, [selectedPackage, viewMode]);
  
  const handleSelect = useCallback((pkg: PackageEntry | null) => {
    setSelectedPackage(pkg);
  }, []);

  const handleLightingChange = useCallback((preset: string) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setLightingPreset(preset);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => viewer.setLightingPreset(preset));
    }
    setCurrentLightingPreset(preset);
  }, [viewMode]);

  const handleRenderModeChange = useCallback((mode: string) => {
    if (viewMode === 'single' && viewerRef.current) {
      viewerRef.current.setRenderMode(mode);
    } else if (viewMode === '4-up' && multiViewRef.current) {
      Object.values(multiViewRef.current).forEach(viewer => viewer.setRenderMode(mode));
    }
    setCurrentRenderMode(mode);
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
          {/* DSL Editor + Command Selector */}
          <div style={{ display: 'flex', flexDirection: 'column', height: '100%', overflow: 'hidden' }}>
            <div style={{ flex: 1, overflow: 'hidden' }}>
              <DSLEditor
                onEvaluate={handleDSLEvaluate}
                onSourceChange={setDslSource}
                externalSource={externalSource}
                onExternalSourceConsumed={handleExternalSourceConsumed}
                isEvaluating={isEvaluating}
                evaluationError={evaluationError}
              />
            </div>
            <CommandSelector
              source={dslSource}
              backendUrl={backendUrl}
              onEvaluate={handleCommandEvaluate}
              isEvaluating={isEvaluating}
            />
          </div>
          
          {/* Chat Panel */}
          <ChatPanel
            dslSource={dslSource}
            onDSLUpdate={(source) => setExternalSource(source)}
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
        </div>
        
        {/* Clipping Controls */}
        <div style={{ padding: '8px 16px' }}>
          <ClippingControls
            onClipPlaneChange={handleClipPlaneChange}
            onResetClipPlanes={handleResetClipPlanes}
            getClippingState={getClippingState}
          />
        </div>
        
        {/* 3D Viewer */}
        <div 
          ref={viewerContainerRef}
          style={{
            flex: 1,
            minWidth: 0,
            minHeight: 0,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            color: '#666',
            fontSize: '16px',
            overflow: 'hidden',
            position: 'relative',
            ...(viewMode === '4-up' ? {
              display: 'grid',
              gridTemplateColumns: '1fr 1fr',
              gridTemplateRows: '1fr 1fr',
              gap: '2px',
              alignItems: 'stretch',
              justifyContent: 'stretch',
              alignContent: 'stretch',
              justifyItems: 'stretch',
            } : {}),
          }}
        >
          {!selectedPackage && !isEvaluating && (
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
            }}>
              ⏳ Evaluating DSL...
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
