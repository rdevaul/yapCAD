/**
 * FileBar Component — file management bar with localStorage persistence
 */

import { useState, useCallback, useRef, useEffect } from 'react';

interface DslFileEntry {
  source: string;
  lastModified: string;
}

interface FileBarProps {
  currentSource: string;
  onFileLoad: (source: string, filename: string) => void;
}

const STORAGE_KEY_FILES = 'yapcad-dsl-files';
const STORAGE_KEY_CURRENT = 'yapcad-dsl-current';

function loadFileMap(): Record<string, DslFileEntry> {
  try {
    const raw = localStorage.getItem(STORAGE_KEY_FILES);
    return raw ? JSON.parse(raw) : {};
  } catch {
    return {};
  }
}

function saveFileMap(map: Record<string, DslFileEntry>) {
  localStorage.setItem(STORAGE_KEY_FILES, JSON.stringify(map));
}

export function FileBar({ currentSource, onFileLoad }: FileBarProps) {
  const [filename, setFilename] = useState(() => {
    return localStorage.getItem(STORAGE_KEY_CURRENT) || 'untitled.dsl';
  });
  const [fileList, setFileList] = useState<string[]>(() => Object.keys(loadFileMap()));
  const [dropdownOpen, setDropdownOpen] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);

  // Close dropdown on outside click
  useEffect(() => {
    const handleClick = (e: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(e.target as Node)) {
        setDropdownOpen(false);
      }
    };
    document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, []);

  const handleSave = useCallback(() => {
    const map = loadFileMap();
    map[filename] = { source: currentSource, lastModified: new Date().toISOString() };
    saveFileMap(map);
    localStorage.setItem(STORAGE_KEY_CURRENT, filename);
    setFileList(Object.keys(map));
  }, [filename, currentSource]);

  const handleSelectFile = useCallback((name: string) => {
    const map = loadFileMap();
    const entry = map[name];
    if (entry) {
      setFilename(name);
      localStorage.setItem(STORAGE_KEY_CURRENT, name);
      onFileLoad(entry.source, name);
    }
    setDropdownOpen(false);
  }, [onFileLoad]);

  const handleLoadFromDisk = useCallback(() => {
    fileInputRef.current?.click();
  }, []);

  const handleFileInputChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const source = reader.result as string;
      const name = file.name;
      // Store in localStorage
      const map = loadFileMap();
      map[name] = { source, lastModified: new Date().toISOString() };
      saveFileMap(map);
      setFilename(name);
      localStorage.setItem(STORAGE_KEY_CURRENT, name);
      setFileList(Object.keys(map));
      onFileLoad(source, name);
    };
    reader.readAsText(file);
    // Reset so same file can be loaded again
    e.target.value = '';
  }, [onFileLoad]);

  const handleDeleteFile = useCallback((name: string, e: React.MouseEvent) => {
    e.stopPropagation();
    const map = loadFileMap();
    delete map[name];
    saveFileMap(map);
    setFileList(Object.keys(map));
  }, []);

  return (
    <div style={{
      display: 'flex',
      alignItems: 'center',
      gap: '4px',
      padding: '4px 8px',
      backgroundColor: '#1e1e36',
      borderBottom: '1px solid #333',
      fontSize: '12px',
      flexShrink: 0,
    }}>
      {/* Filename combo box */}
      <div ref={dropdownRef} style={{ position: 'relative', flex: 1, minWidth: 0 }}>
        <input
          type="text"
          value={filename}
          onChange={e => setFilename(e.target.value)}
          onFocus={() => setDropdownOpen(true)}
          style={{
            width: '100%',
            padding: '4px 8px',
            fontSize: '12px',
            backgroundColor: '#1a1a2e',
            color: '#eee',
            border: '1px solid #444',
            borderRadius: '3px',
            outline: 'none',
            boxSizing: 'border-box',
          }}
        />
        {dropdownOpen && fileList.length > 0 && (
          <div style={{
            position: 'absolute',
            top: '100%',
            left: 0,
            right: 0,
            backgroundColor: '#252540',
            border: '1px solid #444',
            borderRadius: '0 0 3px 3px',
            maxHeight: '200px',
            overflowY: 'auto',
            zIndex: 100,
          }}>
            {fileList.map(name => (
              <div
                key={name}
                onClick={() => handleSelectFile(name)}
                style={{
                  padding: '4px 8px',
                  cursor: 'pointer',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  backgroundColor: name === filename ? '#333355' : 'transparent',
                  color: '#eee',
                }}
                onMouseEnter={e => (e.currentTarget.style.backgroundColor = '#333355')}
                onMouseLeave={e => (e.currentTarget.style.backgroundColor = name === filename ? '#333355' : 'transparent')}
              >
                <span style={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>{name}</span>
                <span
                  onClick={(e) => handleDeleteFile(name, e)}
                  style={{ color: '#888', cursor: 'pointer', marginLeft: '8px', fontSize: '10px' }}
                  title="Delete"
                >✕</span>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Buttons */}
      <button onClick={handleLoadFromDisk} style={btnStyle} title="Load from disk">📂</button>
      <button onClick={handleSave} style={btnStyle} title="Save">💾</button>

      <input
        ref={fileInputRef}
        type="file"
        accept=".dsl,.txt"
        style={{ display: 'none' }}
        onChange={handleFileInputChange}
      />
    </div>
  );
}

const btnStyle: React.CSSProperties = {
  padding: '3px 8px',
  fontSize: '13px',
  backgroundColor: '#333',
  color: '#eee',
  border: '1px solid #444',
  borderRadius: '3px',
  cursor: 'pointer',
  flexShrink: 0,
};
