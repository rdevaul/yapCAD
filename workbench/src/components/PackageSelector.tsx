/**
 * PackageSelector - Gallery view for browsing and selecting yapCAD packages.
 */

import React, { useCallback, useReducer, useRef, useState } from 'react';
import type { PackageEntry, PackageSelectorState, PackageSelectorAction } from '../types/package';
import { parsePackageFile, generateThumbnail, loadPackageFromUrl, loadPackageFromPath } from '../utils/packageParser';

// Styles (inline for now, can extract to CSS module later)
const styles = {
  container: {
    display: 'flex',
    flexDirection: 'column' as const,
    height: '100%',
    backgroundColor: '#1a1a2e',
    color: '#eee',
  },
  header: {
    padding: '16px',
    borderBottom: '1px solid #333',
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
  },
  title: {
    margin: 0,
    fontSize: '18px',
    fontWeight: 600,
  },
  gallery: {
    flex: 1,
    padding: '16px',
    display: 'grid',
    gridTemplateColumns: 'repeat(auto-fill, minmax(150px, 1fr))',
    gap: '16px',
    overflowY: 'auto' as const,
    alignContent: 'start',
  },
  card: {
    display: 'flex',
    flexDirection: 'column' as const,
    backgroundColor: '#252540',
    borderRadius: '8px',
    overflow: 'hidden',
    cursor: 'pointer',
    transition: 'transform 0.2s, box-shadow 0.2s',
    border: '2px solid transparent',
  },
  cardSelected: {
    border: '2px solid #3b82f6',
    boxShadow: '0 0 12px rgba(59, 130, 246, 0.4)',
  },
  cardHover: {
    transform: 'translateY(-2px)',
    boxShadow: '0 4px 12px rgba(0, 0, 0, 0.3)',
  },
  thumbnail: {
    width: '100%',
    aspectRatio: '1',
    objectFit: 'cover' as const,
    backgroundColor: '#1a1a2e',
  },
  cardInfo: {
    padding: '12px',
  },
  cardName: {
    margin: 0,
    fontSize: '14px',
    fontWeight: 500,
    whiteSpace: 'nowrap' as const,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
  },
  cardVersion: {
    margin: '4px 0 0',
    fontSize: '12px',
    color: '#888',
  },
  addCard: {
    display: 'flex',
    flexDirection: 'column' as const,
    alignItems: 'center',
    justifyContent: 'center',
    backgroundColor: '#252540',
    borderRadius: '8px',
    border: '2px dashed #444',
    cursor: 'pointer',
    minHeight: '150px',
    transition: 'border-color 0.2s, background-color 0.2s',
  },
  addCardHover: {
    borderColor: '#3b82f6',
    backgroundColor: '#2a2a4a',
  },
  addIcon: {
    fontSize: '32px',
    marginBottom: '8px',
    color: '#666',
  },
  addText: {
    fontSize: '14px',
    color: '#888',
  },
  dropOverlay: {
    position: 'absolute' as const,
    inset: 0,
    backgroundColor: 'rgba(59, 130, 246, 0.2)',
    border: '3px dashed #3b82f6',
    borderRadius: '8px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    fontSize: '18px',
    fontWeight: 500,
    color: '#3b82f6',
    zIndex: 10,
  },
  error: {
    padding: '8px 16px',
    backgroundColor: '#7f1d1d',
    color: '#fca5a5',
    fontSize: '14px',
  },
};

// Reducer for package state management
function packageReducer(
  state: PackageSelectorState,
  action: PackageSelectorAction
): PackageSelectorState {
  switch (action.type) {
    case 'ADD_PACKAGE':
      // Check for duplicate UUID
      if (state.packages.some(p => p.uuid === action.payload.uuid)) {
        // Update existing instead
        return {
          ...state,
          packages: state.packages.map(p =>
            p.uuid === action.payload.uuid ? action.payload : p
          ),
        };
      }
      return {
        ...state,
        packages: [...state.packages, action.payload],
      };
    
    case 'REMOVE_PACKAGE':
      return {
        ...state,
        packages: state.packages.filter(p => p.uuid !== action.payload),
        selectedId: state.selectedId === action.payload ? null : state.selectedId,
      };
    
    case 'SELECT_PACKAGE':
      return { ...state, selectedId: action.payload };
    
    case 'UPDATE_PACKAGE':
      return {
        ...state,
        packages: state.packages.map(p =>
          p.uuid === action.payload.uuid ? { ...p, ...action.payload } : p
        ),
      };
    
    case 'SET_LOADING':
      return { ...state, loading: action.payload };
    
    case 'SET_ERROR':
      return { ...state, error: action.payload };
    
    default:
      return state;
  }
}

interface PackageSelectorProps {
  onSelect?: (pkg: PackageEntry | null) => void;
  initialPackages?: PackageEntry[];
}

export function PackageSelector({ onSelect, initialPackages = [] }: PackageSelectorProps) {
  const [state, dispatch] = useReducer(packageReducer, {
    packages: initialPackages,
    selectedId: null,
    loading: false,
  });
  
  const [isDragging, setIsDragging] = React.useState(false);
  const [hoveredCard, setHoveredCard] = React.useState<string | null>(null);
  const [showUrlInput, setShowUrlInput] = useState(false);
  const [urlInput, setUrlInput] = useState('');
  const fileInputRef = useRef<HTMLInputElement>(null);
  
  // Handle file selection
  const handleFiles = useCallback(async (files: FileList | File[]) => {
    dispatch({ type: 'SET_LOADING', payload: true });
    dispatch({ type: 'SET_ERROR', payload: undefined });
    
    for (const file of Array.from(files)) {
      if (!file.name.endsWith('.ycpkg')) {
        dispatch({ type: 'SET_ERROR', payload: `Invalid file type: ${file.name}` });
        continue;
      }
      
      try {
        const pkg = await parsePackageFile(file);
        
        // Generate thumbnail if missing
        if (!pkg.thumbnail && pkg.geometry) {
          pkg.thumbnail = await generateThumbnail(pkg.geometry);
        }
        
        dispatch({ type: 'ADD_PACKAGE', payload: pkg });
      } catch (err) {
        console.error('Failed to parse package:', err);
        dispatch({ 
          type: 'SET_ERROR', 
          payload: `Failed to parse ${file.name}: ${err instanceof Error ? err.message : 'Unknown error'}` 
        });
      }
    }
    
    dispatch({ type: 'SET_LOADING', payload: false });
  }, []);
  
  // Drag and drop handlers
  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(true);
  }, []);
  
  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);
  
  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
    
    if (e.dataTransfer.files.length > 0) {
      handleFiles(e.dataTransfer.files);
    }
  }, [handleFiles]);
  
  // File input handler
  const handleFileInput = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      handleFiles(e.target.files);
      e.target.value = ''; // Reset for re-selection
    }
  }, [handleFiles]);
  
  // Package selection
  const handleSelect = useCallback((pkg: PackageEntry) => {
    const newSelectedId = state.selectedId === pkg.uuid ? null : pkg.uuid;
    dispatch({ type: 'SELECT_PACKAGE', payload: newSelectedId });
    onSelect?.(newSelectedId ? pkg : null);
  }, [state.selectedId, onSelect]);
  
  // Open file picker
  const openFilePicker = useCallback(() => {
    fileInputRef.current?.click();
  }, []);
  
  // Load from URL or local path
  const handleUrlLoad = useCallback(async () => {
    if (!urlInput.trim()) return;
    
    dispatch({ type: 'SET_LOADING', payload: true });
    dispatch({ type: 'SET_ERROR', payload: undefined });
    
    try {
      let pkg: PackageEntry;
      
      if (urlInput.startsWith('http://') || urlInput.startsWith('https://')) {
        // Remote URL
        pkg = await loadPackageFromUrl(urlInput);
      } else {
        // Local path - use backend proxy
        pkg = await loadPackageFromPath(urlInput);
      }
      
      // Generate thumbnail if missing
      if (!pkg.thumbnail && pkg.geometry) {
        pkg.thumbnail = await generateThumbnail(pkg.geometry);
      }
      
      dispatch({ type: 'ADD_PACKAGE', payload: pkg });
      setUrlInput('');
      setShowUrlInput(false);
    } catch (err) {
      console.error('Failed to load from URL/path:', err);
      dispatch({
        type: 'SET_ERROR',
        payload: `Failed to load: ${err instanceof Error ? err.message : 'Unknown error'}`
      });
    }
    
    dispatch({ type: 'SET_LOADING', payload: false });
  }, [urlInput]);
  
  return (
    <div 
      style={styles.container}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      {/* Header */}
      <div style={styles.header}>
        <h2 style={styles.title}>Packages</h2>
        <span style={{ fontSize: '14px', color: '#888' }}>
          {state.packages.length} package{state.packages.length !== 1 ? 's' : ''}
        </span>
      </div>
      
      {/* Error message */}
      {state.error && (
        <div style={styles.error}>{state.error}</div>
      )}
      
      {/* Gallery grid */}
      <div style={{ ...styles.gallery, position: 'relative' }}>
        {/* Drop overlay */}
        {isDragging && (
          <div style={styles.dropOverlay}>
            Drop .ycpkg files here
          </div>
        )}
        
        {/* Package cards */}
        {state.packages.map(pkg => (
          <div
            key={pkg.uuid}
            style={{
              ...styles.card,
              ...(state.selectedId === pkg.uuid ? styles.cardSelected : {}),
              ...(hoveredCard === pkg.uuid ? styles.cardHover : {}),
            }}
            onClick={() => handleSelect(pkg)}
            onMouseEnter={() => setHoveredCard(pkg.uuid)}
            onMouseLeave={() => setHoveredCard(null)}
          >
            {pkg.thumbnail ? (
              <img 
                src={pkg.thumbnail} 
                alt={pkg.name}
                style={styles.thumbnail}
              />
            ) : (
              <div style={{ ...styles.thumbnail, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                <span style={{ fontSize: '32px' }}>📦</span>
              </div>
            )}
            <div style={styles.cardInfo}>
              <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                <p style={{ ...styles.cardName, flex: 1 }} title={pkg.name}>{pkg.name}</p>
                <button
                  onClick={(e) => {
                    e.stopPropagation();
                    dispatch({ type: 'REMOVE_PACKAGE', payload: pkg.uuid });
                    if (state.selectedId === pkg.uuid) {
                      onSelect?.(null);
                    }
                  }}
                  title="Remove package"
                  style={{
                    backgroundColor: '#333',
                    border: 'none',
                    color: '#aaa',
                    cursor: 'pointer',
                    fontSize: '16px',
                    padding: '2px 6px',
                    borderRadius: '50%',
                    lineHeight: 1,
                    flexShrink: 0,
                    transition: 'background-color 0.2s, color 0.2s',
                  }}
                  onMouseEnter={(e) => {
                    e.currentTarget.style.backgroundColor = '#ff4444';
                    e.currentTarget.style.color = '#fff';
                  }}
                  onMouseLeave={(e) => {
                    e.currentTarget.style.backgroundColor = '#333';
                    e.currentTarget.style.color = '#aaa';
                  }}
                >
                  ✕
                </button>
              </div>
              <p style={styles.cardVersion}>v{pkg.version}</p>
            </div>
          </div>
        ))}
        
        {/* Add package card */}
        <div
          style={{
            ...styles.addCard,
            ...(hoveredCard === 'add' ? styles.addCardHover : {}),
          }}
          onClick={openFilePicker}
          onMouseEnter={() => setHoveredCard('add')}
          onMouseLeave={() => setHoveredCard(null)}
        >
          <span style={styles.addIcon}>+</span>
          <span style={styles.addText}>Add File</span>
        </div>
        
        {/* Load from URL/path card */}
        <div
          style={{
            ...styles.addCard,
            ...(hoveredCard === 'url' ? styles.addCardHover : {}),
          }}
          onClick={() => setShowUrlInput(!showUrlInput)}
          onMouseEnter={() => setHoveredCard('url')}
          onMouseLeave={() => setHoveredCard(null)}
        >
          <span style={styles.addIcon}>🔗</span>
          <span style={styles.addText}>URL/Path</span>
        </div>
        
        {/* Hidden file input */}
        <input
          ref={fileInputRef}
          type="file"
          accept=".ycpkg"
          multiple
          style={{ display: 'none' }}
          onChange={handleFileInput}
        />
      </div>
      
      {/* URL/Path input panel */}
      {showUrlInput && (
        <div style={{
          padding: '12px 16px',
          borderTop: '1px solid #333',
          backgroundColor: '#252540',
        }}>
          <div style={{ fontSize: '12px', color: '#888', marginBottom: '8px' }}>
            Enter URL or local path (e.g., ~/Projects/yapCAD/my.ycpkg)
          </div>
          <div style={{ display: 'flex', gap: '8px' }}>
            <input
              type="text"
              value={urlInput}
              onChange={(e) => setUrlInput(e.target.value)}
              onKeyDown={(e) => e.key === 'Enter' && handleUrlLoad()}
              placeholder="https://... or /path/to/package.ycpkg"
              style={{
                flex: 1,
                padding: '8px 12px',
                fontSize: '14px',
                backgroundColor: '#1a1a2e',
                border: '1px solid #444',
                borderRadius: '4px',
                color: '#eee',
                outline: 'none',
              }}
              autoFocus
            />
            <button
              onClick={handleUrlLoad}
              disabled={!urlInput.trim() || state.loading}
              style={{
                padding: '8px 16px',
                fontSize: '14px',
                backgroundColor: urlInput.trim() ? '#3b82f6' : '#333',
                color: '#eee',
                border: 'none',
                borderRadius: '4px',
                cursor: urlInput.trim() ? 'pointer' : 'default',
                opacity: state.loading ? 0.5 : 1,
              }}
            >
              {state.loading ? 'Loading...' : 'Load'}
            </button>
          </div>
        </div>
      )}
    </div>
  );
}

export default PackageSelector;
