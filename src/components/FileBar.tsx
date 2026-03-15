/**
 * FileBar Component — file management bar with localStorage + WBS remote file support.
 *
 * File namespaces:
 *   Local (default)  — stored in browser localStorage, loaded/saved by the user's browser.
 *   WBS remote       — paths prefixed with "wbs://", fetched/saved through the WBS server
 *                      at /wbs/api/dsl/<rel>.  Optimistic locking via ETag / If-Match.
 *
 * Conflict handling:
 *   On PUT 409, the user is shown a dialog: overwrite, discard, or cancel.
 */

import { useState, useCallback, useRef, useEffect } from 'react';

// ── WBS helpers ───────────────────────────────────────────────────────────────

const WBS_SCHEME   = 'wbs://';
const WBS_PROJECT  = 'agentic-1';          // currently the only project
const WBS_API_BASE = '/wbs/api/dsl';       // proxied through Vite → localhost:8765

/** Return true if this filename lives in the WBS namespace. */
function isWbsPath(name: string): boolean {
  return name.startsWith(WBS_SCHEME);
}

/**
 * Convert a wbs:// URI to the API sub-path used in fetch calls.
 * e.g. "wbs://agentic-1/designs/tanks/dsl/tank.dsl" → "designs/tanks/dsl/tank.dsl"
 */
function wbsToRel(wbsPath: string): string {
  // strip "wbs://<project>/"
  const prefix = `${WBS_SCHEME}${WBS_PROJECT}/`;
  return wbsPath.startsWith(prefix) ? wbsPath.slice(prefix.length) : wbsPath.slice(WBS_SCHEME.length);
}

interface WbsResult {
  content: string;
  etag: string;
  modified: string;
}

interface WbsConflict {
  currentContent: string;
  currentEtag: string;
  modified: string;
}

async function wbsFetch(wbsPath: string, wbsToken: string): Promise<WbsResult> {
  const rel = wbsToRel(wbsPath);
  const headers: Record<string, string> = {};
  if (wbsToken) headers['X-Dashboard-Token'] = wbsToken;
  const resp = await fetch(`${WBS_API_BASE}/${rel}`, { headers });
  if (!resp.ok) throw new Error(`WBS fetch failed: ${resp.status} ${await resp.text()}`);
  const content  = await resp.text();
  const etag     = resp.headers.get('ETag')?.replace(/^"|"$/g, '') ?? '';
  const modified = resp.headers.get('Last-Modified') ?? '';
  return { content, etag, modified };
}

/** Returns new ETag on success, throws on error, returns WbsConflict on 409. */
async function wbsPut(
  wbsPath: string,
  content: string,
  etag: string,
  wbsToken: string,
): Promise<{ ok: true; etag: string } | { ok: false; conflict: WbsConflict }> {
  const rel = wbsToRel(wbsPath);
  const headers: Record<string, string> = {
    'Content-Type': 'text/plain; charset=utf-8',
    'If-Match': `"${etag}"`,
  };
  if (wbsToken) headers['X-Dashboard-Token'] = wbsToken;
  const resp = await fetch(`${WBS_API_BASE}/${rel}`, {
    method: 'PUT',
    headers,
    body: content,
  });
  if (resp.ok) {
    const json = await resp.json();
    return { ok: true, etag: json.etag };
  }
  if (resp.status === 409) {
    const json = await resp.json();
    return {
      ok: false,
      conflict: {
        currentContent: json.current_content,
        currentEtag:    json.current_etag,
        modified:       json.modified,
      },
    };
  }
  throw new Error(`WBS PUT failed: ${resp.status} ${await resp.text()}`);
}

// ── Local storage helpers ─────────────────────────────────────────────────────

interface DslFileEntry {
  source: string;
  lastModified: string;
}

const STORAGE_KEY_FILES   = 'yapcad-dsl-files';
const STORAGE_KEY_CURRENT = 'yapcad-dsl-current';
const STORAGE_KEY_WBSTOKEN = 'yapcad-wbs-token';

function loadFileMap(): Record<string, DslFileEntry> {
  try {
    const raw = localStorage.getItem(STORAGE_KEY_FILES);
    return raw ? JSON.parse(raw) : {};
  } catch { return {}; }
}

function saveFileMap(map: Record<string, DslFileEntry>) {
  localStorage.setItem(STORAGE_KEY_FILES, JSON.stringify(map));
}

// ── Conflict dialog ───────────────────────────────────────────────────────────

interface ConflictDialogProps {
  filename: string;
  conflict: WbsConflict;
  onOverwrite: () => void;
  onDiscard:   () => void;
  onCancel:    () => void;
}

function ConflictDialog({ filename, conflict, onOverwrite, onDiscard, onCancel }: ConflictDialogProps) {
  return (
    <div style={dlgStyles.overlay}>
      <div style={dlgStyles.box}>
        <div style={dlgStyles.title}>⚠️ Edit Conflict</div>
        <div style={dlgStyles.body}>
          <p style={dlgStyles.p}>
            <strong>{filename}</strong> was modified on the server since you opened it.
          </p>
          <p style={dlgStyles.p}>
            Last server save: <code style={dlgStyles.code}>{conflict.modified}</code>
          </p>
          <p style={dlgStyles.p}>Choose how to resolve:</p>
        </div>
        <div style={dlgStyles.actions}>
          <button style={dlgStyles.cancelBtn}    onClick={onCancel}>Cancel</button>
          <button style={dlgStyles.discardBtn}   onClick={onDiscard}>Discard my changes</button>
          <button style={dlgStyles.overwriteBtn} onClick={onOverwrite}>Overwrite server version</button>
        </div>
      </div>
    </div>
  );
}

const dlgStyles = {
  overlay:      { position: 'fixed' as const, inset: 0, backgroundColor: 'rgba(0,0,0,0.7)',
                  display: 'flex', alignItems: 'center', justifyContent: 'center', zIndex: 9999 },
  box:          { backgroundColor: '#1e1e36', border: '1px solid #555', borderRadius: '8px',
                  padding: '20px 24px', maxWidth: '440px', width: '90%', color: '#ddd' },
  title:        { fontSize: '15px', fontWeight: 700, marginBottom: '12px', color: '#f59e0b' },
  body:         { fontSize: '13px', lineHeight: 1.6 },
  p:            { margin: '0 0 8px' },
  code:         { fontSize: '11px', backgroundColor: '#111', padding: '2px 5px',
                  borderRadius: '3px', color: '#adf' },
  actions:      { display: 'flex', gap: '8px', justifyContent: 'flex-end', marginTop: '18px',
                  flexWrap: 'wrap' as const },
  cancelBtn:    { padding: '6px 14px', fontSize: '12px', backgroundColor: '#333', color: '#aaa',
                  border: '1px solid #555', borderRadius: '4px', cursor: 'pointer' } as React.CSSProperties,
  discardBtn:   { padding: '6px 14px', fontSize: '12px', backgroundColor: '#7f1d1d', color: '#fca5a5',
                  border: '1px solid #991b1b', borderRadius: '4px', cursor: 'pointer' } as React.CSSProperties,
  overwriteBtn: { padding: '6px 14px', fontSize: '12px', backgroundColor: '#1d4ed8', color: '#bfdbfe',
                  border: '1px solid #2563eb', borderRadius: '4px', cursor: 'pointer' } as React.CSSProperties,
};

// ── Main component ────────────────────────────────────────────────────────────

export interface FileBarProps {
  currentSource: string;
  onFileLoad: (source: string, filename: string) => void;
  /** Called whenever unsaved-changes state flips, so App can gate ?file= loads. */
  onUnsavedChange?: (hasChanges: boolean) => void;
  /** A wbs:// path passed from App (e.g. from ?file= query param). FileBar loads it then calls onPendingWbsFileConsumed. */
  pendingWbsFile?: string | null;
  onPendingWbsFileConsumed?: () => void;
}

export function FileBar({ currentSource, onFileLoad, onUnsavedChange, pendingWbsFile, onPendingWbsFileConsumed }: FileBarProps) {
  const [filename, setFilename] = useState(() =>
    localStorage.getItem(STORAGE_KEY_CURRENT) || 'untitled.dsl');
  const [fileList, setFileList]     = useState<string[]>(() => Object.keys(loadFileMap()));
  const [dropdownOpen, setDropdownOpen] = useState(false);
  const [wbsToken, setWbsToken]     = useState(() =>
    localStorage.getItem(STORAGE_KEY_WBSTOKEN) || '');
  const [showTokenInput, setShowTokenInput] = useState(false);
  const [status, setStatus]         = useState<string>('');
  const [isLoading, setIsLoading]   = useState(false);

  // Unsaved-changes tracking
  const [savedSource, setSavedSource] = useState(currentSource);
  const hasUnsaved = currentSource !== savedSource;

  // ETag for current WBS file (updated on every successful GET or PUT)
  const etagRef = useRef<string>('');

  // Conflict dialog state
  const [conflict, setConflict] = useState<WbsConflict | null>(null);
  const pendingSaveRef = useRef<string>('');   // content we tried to save when conflict hit

  const fileInputRef  = useRef<HTMLInputElement>(null);
  const dropdownRef   = useRef<HTMLDivElement>(null);

  // Bubble unsaved state to parent
  useEffect(() => {
    onUnsavedChange?.(hasUnsaved);
  }, [hasUnsaved, onUnsavedChange]);

  // Load a wbs:// file when App requests it (e.g. from ?file= query param)
  useEffect(() => {
    if (pendingWbsFile) {
      handleWbsLoad(pendingWbsFile).finally(() => onPendingWbsFileConsumed?.());
    }
  // handleWbsLoad is stable (useCallback); pendingWbsFile drives this
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pendingWbsFile]);

  // Close dropdown on outside click
  useEffect(() => {
    const handleClick = (e: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(e.target as Node))
        setDropdownOpen(false);
    };
    document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, []);

  // ── Save ─────────────────────────────────────────────────────────────────────

  const handleSave = useCallback(async () => {
    if (isWbsPath(filename)) {
      // ── WBS save ────────────────────────────────────────────────────────────
      setIsLoading(true);
      setStatus('Saving…');
      try {
        const result = await wbsPut(filename, currentSource, etagRef.current, wbsToken);
        if (result.ok) {
          etagRef.current = result.etag;
          setSavedSource(currentSource);
          setStatus('✓ Saved to WBS');
          setTimeout(() => setStatus(''), 2000);
        } else {
          // 409 conflict — show dialog
          pendingSaveRef.current = currentSource;
          setConflict(result.conflict);
          setStatus('');
        }
      } catch (e) {
        setStatus(`✗ ${e instanceof Error ? e.message : e}`);
      } finally {
        setIsLoading(false);
      }
    } else {
      // ── Local save ──────────────────────────────────────────────────────────
      const map = loadFileMap();
      map[filename] = { source: currentSource, lastModified: new Date().toISOString() };
      saveFileMap(map);
      localStorage.setItem(STORAGE_KEY_CURRENT, filename);
      setFileList(Object.keys(map));
      setSavedSource(currentSource);
      setStatus('✓ Saved');
      setTimeout(() => setStatus(''), 1500);
    }
  }, [filename, currentSource, wbsToken]);

  // ── Load (local dropdown) ────────────────────────────────────────────────────

  const handleSelectFile = useCallback((name: string) => {
    const map = loadFileMap();
    const entry = map[name];
    if (entry) {
      setFilename(name);
      localStorage.setItem(STORAGE_KEY_CURRENT, name);
      setSavedSource(entry.source);
      etagRef.current = '';
      onFileLoad(entry.source, name);
    }
    setDropdownOpen(false);
  }, [onFileLoad]);

  // ── Load from disk ───────────────────────────────────────────────────────────

  const handleLoadFromDisk = useCallback(() => { fileInputRef.current?.click(); }, []);

  const handleFileInputChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const source = reader.result as string;
      const name   = file.name;
      const map    = loadFileMap();
      map[name] = { source, lastModified: new Date().toISOString() };
      saveFileMap(map);
      setFilename(name);
      localStorage.setItem(STORAGE_KEY_CURRENT, name);
      setFileList(Object.keys(map));
      setSavedSource(source);
      etagRef.current = '';
      onFileLoad(source, name);
    };
    reader.readAsText(file);
    e.target.value = '';
  }, [onFileLoad]);

  // ── Load from WBS ─────────────────────────────────────────────────────────────

  /**
   * Load a wbs:// path into the editor.  Called when:
   *   (a) user types a wbs:// path in the filename box and presses Enter
   *   (b) App passes a ?file= query param on startup
   */
  const handleWbsLoad = useCallback(async (wbsPath: string) => {
    setIsLoading(true);
    setStatus('Loading from WBS…');
    try {
      const result = await wbsFetch(wbsPath, wbsToken);
      etagRef.current = result.etag;
      setFilename(wbsPath);
      localStorage.setItem(STORAGE_KEY_CURRENT, wbsPath);
      setSavedSource(result.content);
      onFileLoad(result.content, wbsPath);
      setStatus(`✓ Loaded (etag: ${result.etag.slice(0, 8)})`);
      setTimeout(() => setStatus(''), 2500);
    } catch (e) {
      setStatus(`✗ ${e instanceof Error ? e.message : e}`);
    } finally {
      setIsLoading(false);
    }
  }, [wbsToken, onFileLoad]);

  // Expose wbsLoad so App.tsx can call it via ref
  // (Parent passes a callback via onFileLoad; we use a different channel for external triggers)

  // ── Filename input key handler ───────────────────────────────────────────────

  const handleFilenameKeyDown = useCallback((e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter' && isWbsPath(filename)) {
      handleWbsLoad(filename);
    }
  }, [filename, handleWbsLoad]);

  // ── Delete (local only) ──────────────────────────────────────────────────────

  const handleDeleteFile = useCallback((name: string, e: React.MouseEvent) => {
    e.stopPropagation();
    const map = loadFileMap();
    delete map[name];
    saveFileMap(map);
    setFileList(Object.keys(map));
  }, []);

  // ── Conflict resolution ──────────────────────────────────────────────────────

  const handleConflictOverwrite = useCallback(async () => {
    // Force-write our content using the server's current ETag
    setConflict(null);
    setIsLoading(true);
    setStatus('Overwriting…');
    try {
      const result = await wbsPut(filename, pendingSaveRef.current, conflict!.currentEtag, wbsToken);
      if (result.ok) {
        etagRef.current = result.etag;
        setSavedSource(pendingSaveRef.current);
        setStatus('✓ Overwritten');
        setTimeout(() => setStatus(''), 2000);
      } else {
        // Another conflict — rarer but possible
        setConflict(result.conflict);
        setStatus('');
      }
    } catch (e) {
      setStatus(`✗ ${e instanceof Error ? e.message : e}`);
    } finally {
      setIsLoading(false);
    }
  }, [filename, conflict, wbsToken]);

  const handleConflictDiscard = useCallback(() => {
    if (!conflict) return;
    setConflict(null);
    etagRef.current = conflict.currentEtag;
    setSavedSource(conflict.currentContent);
    onFileLoad(conflict.currentContent, filename);
    setStatus('Loaded server version');
    setTimeout(() => setStatus(''), 2000);
  }, [conflict, filename, onFileLoad]);

  const handleConflictCancel = useCallback(() => {
    setConflict(null);
    setStatus('Save cancelled');
    setTimeout(() => setStatus(''), 1500);
  }, []);

  // ── Render ────────────────────────────────────────────────────────────────────

  const wbs = isWbsPath(filename);

  return (
    <>
      {conflict && (
        <ConflictDialog
          filename={filename}
          conflict={conflict}
          onOverwrite={handleConflictOverwrite}
          onDiscard={handleConflictDiscard}
          onCancel={handleConflictCancel}
        />
      )}

      <div style={{
        display: 'flex',
        alignItems: 'center',
        gap: '4px',
        padding: '4px 8px',
        backgroundColor: wbs ? '#0f1a2e' : '#1e1e36',
        borderBottom: `1px solid ${wbs ? '#1e4a7a' : '#333'}`,
        fontSize: '12px',
        flexShrink: 0,
      }}>
        {/* WBS indicator badge */}
        {wbs && (
          <span title="WBS-managed file — saved to server" style={{
            fontSize: '10px', fontWeight: 700, color: '#60a5fa',
            backgroundColor: '#1e3a5f', padding: '1px 5px', borderRadius: '3px',
            letterSpacing: '0.04em', flexShrink: 0,
          }}>WBS</span>
        )}

        {/* Unsaved indicator */}
        {hasUnsaved && (
          <span title="Unsaved changes" style={{
            color: '#f59e0b', fontSize: '14px', lineHeight: 1, flexShrink: 0,
          }}>●</span>
        )}

        {/* Filename combo box */}
        <div ref={dropdownRef} style={{ position: 'relative', flex: 1, minWidth: 0 }}>
          <input
            type="text"
            value={filename}
            onChange={e => { setFilename(e.target.value); etagRef.current = ''; }}
            onFocus={() => !wbs && setDropdownOpen(true)}
            onKeyDown={handleFilenameKeyDown}
            title={wbs ? 'Press Enter to load this WBS file' : 'Filename — click to select a saved file'}
            style={{
              width: '100%',
              padding: '4px 8px',
              fontSize: '12px',
              backgroundColor: wbs ? '#0c1424' : '#1a1a2e',
              color: wbs ? '#93c5fd' : '#eee',
              border: `1px solid ${wbs ? '#1e4a7a' : '#444'}`,
              borderRadius: '3px',
              outline: 'none',
              boxSizing: 'border-box' as const,
            }}
          />
          {dropdownOpen && fileList.length > 0 && (
            <div style={{
              position: 'absolute', top: '100%', left: 0, right: 0,
              backgroundColor: '#252540', border: '1px solid #444',
              borderRadius: '0 0 3px 3px', maxHeight: '200px',
              overflowY: 'auto', zIndex: 100,
            }}>
              {fileList.map(name => (
                <div
                  key={name}
                  onClick={() => handleSelectFile(name)}
                  style={{
                    padding: '4px 8px', cursor: 'pointer',
                    display: 'flex', alignItems: 'center',
                    justifyContent: 'space-between',
                    backgroundColor: name === filename ? '#333355' : 'transparent',
                    color: '#eee',
                  }}
                  onMouseEnter={e => (e.currentTarget.style.backgroundColor = '#333355')}
                  onMouseLeave={e => (e.currentTarget.style.backgroundColor =
                    name === filename ? '#333355' : 'transparent')}
                >
                  <span style={{ overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
                    {name}
                  </span>
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

        {/* Action buttons */}
        {!wbs && (
          <button onClick={handleLoadFromDisk} style={btnStyle} title="Load from disk" disabled={isLoading}>
            📂
          </button>
        )}
        {wbs && (
          <button
            onClick={() => handleWbsLoad(filename)}
            style={{ ...btnStyle, color: '#60a5fa' }}
            title="Reload from WBS server"
            disabled={isLoading}
          >↓</button>
        )}
        <button
          onClick={handleSave}
          style={{ ...btnStyle, ...(hasUnsaved ? { color: '#f59e0b', borderColor: '#78350f' } : {}) }}
          title={wbs ? 'Save to WBS server (Ctrl+S)' : 'Save to browser storage (Ctrl+S)'}
          disabled={isLoading}
        >💾</button>

        {/* WBS token config */}
        <button
          onClick={() => setShowTokenInput(t => !t)}
          style={{ ...btnStyle, color: wbsToken ? '#4ade80' : '#888' }}
          title={wbsToken ? 'WBS token configured' : 'Set WBS auth token'}
        >🔑</button>

        {/* Status */}
        {status && (
          <span style={{ fontSize: '11px', color: '#888', flexShrink: 0, maxWidth: '140px',
                         overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
            {status}
          </span>
        )}

        <input
          ref={fileInputRef}
          type="file"
          accept=".dsl,.txt"
          style={{ display: 'none' }}
          onChange={handleFileInputChange}
        />
      </div>

      {/* Token input panel */}
      {showTokenInput && (
        <div style={{
          padding: '4px 8px', backgroundColor: '#111828',
          borderBottom: '1px solid #1e4a7a', display: 'flex', gap: '6px', alignItems: 'center',
        }}>
          <span style={{ fontSize: '11px', color: '#60a5fa', flexShrink: 0 }}>WBS token:</span>
          <input
            type="password"
            value={wbsToken}
            onChange={e => {
              setWbsToken(e.target.value);
              localStorage.setItem(STORAGE_KEY_WBSTOKEN, e.target.value);
            }}
            placeholder="X-Dashboard-Token value (leave blank if none)"
            style={{
              flex: 1, padding: '3px 7px', fontSize: '11px',
              backgroundColor: '#1a1a2e', color: '#ddd',
              border: '1px solid #334', borderRadius: '3px', outline: 'none',
            }}
          />
          <button
            onClick={() => setShowTokenInput(false)}
            style={{ ...btnStyle, fontSize: '11px', padding: '2px 8px' }}
          >✓</button>
        </div>
      )}
    </>
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
