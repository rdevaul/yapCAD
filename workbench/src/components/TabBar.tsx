/**
 * TabBar — Multi-file tab strip for the DSL editor.
 *
 * Shows one tab per open DSL file with:
 *   - Editable filename (double-click to rename)
 *   - Close (×) button (disabled when only one tab)
 *   - "+" button to add a new tab
 */

import { useState, useRef, useCallback } from 'react';

interface Tab {
  id: string;
  filename: string;
}

interface TabBarProps {
  tabs: Tab[];
  activeTabId: string;
  onSwitch: (tabId: string) => void;
  onClose: (tabId: string) => void;
  onAdd: () => void;
  onRename: (tabId: string, newName: string) => void;
}

export function TabBar({ tabs, activeTabId, onSwitch, onClose, onAdd, onRename }: TabBarProps) {
  const [editingTabId, setEditingTabId] = useState<string | null>(null);
  const [editValue, setEditValue] = useState('');
  const inputRef = useRef<HTMLInputElement>(null);

  const handleDoubleClick = useCallback((tab: Tab) => {
    setEditingTabId(tab.id);
    setEditValue(tab.filename);
    setTimeout(() => inputRef.current?.select(), 0);
  }, []);

  const handleRenameSubmit = useCallback(() => {
    if (editingTabId && editValue.trim()) {
      onRename(editingTabId, editValue.trim());
    }
    setEditingTabId(null);
  }, [editingTabId, editValue, onRename]);

  return (
    <div style={styles.container}>
      <div style={styles.tabStrip}>
        {tabs.map(tab => {
          const isActive = tab.id === activeTabId;
          const isEditing = tab.id === editingTabId;

          return (
            <div
              key={tab.id}
              style={{
                ...styles.tab,
                ...(isActive ? styles.activeTab : {}),
              }}
              onClick={() => onSwitch(tab.id)}
              title={tab.filename}
            >
              {isEditing ? (
                <input
                  ref={inputRef}
                  style={styles.renameInput}
                  value={editValue}
                  onChange={e => setEditValue(e.target.value)}
                  onBlur={handleRenameSubmit}
                  onKeyDown={e => {
                    if (e.key === 'Enter') handleRenameSubmit();
                    if (e.key === 'Escape') setEditingTabId(null);
                  }}
                  onClick={e => e.stopPropagation()}
                  autoFocus
                />
              ) : (
                <span
                  style={styles.tabName}
                  onDoubleClick={(e) => {
                    e.stopPropagation();
                    handleDoubleClick(tab);
                  }}
                >
                  {tab.filename}
                </span>
              )}

              {tabs.length > 1 && (
                <button
                  style={styles.closeBtn}
                  onClick={(e) => {
                    e.stopPropagation();
                    onClose(tab.id);
                  }}
                  title="Close tab"
                >
                  ×
                </button>
              )}
            </div>
          );
        })}
      </div>

      <button
        style={styles.addBtn}
        onClick={onAdd}
        title="New tab"
      >
        +
      </button>
    </div>
  );
}

// ── Styles ────────────────────────────────────────────────────────────────────

const styles = {
  container: {
    display: 'flex',
    alignItems: 'stretch',
    backgroundColor: '#15152a',
    borderBottom: '1px solid #2a2a4a',
    flexShrink: 0,
    minHeight: '30px',
    overflow: 'hidden',
  } as React.CSSProperties,

  tabStrip: {
    display: 'flex',
    flex: 1,
    overflowX: 'auto',
    overflowY: 'hidden',
    scrollbarWidth: 'none',
    msOverflowStyle: 'none',
  } as React.CSSProperties,

  tab: {
    display: 'flex',
    alignItems: 'center',
    gap: '4px',
    padding: '4px 10px',
    fontSize: '11px',
    color: '#888',
    cursor: 'pointer',
    borderRight: '1px solid #1a1a30',
    whiteSpace: 'nowrap',
    userSelect: 'none',
    transition: 'background-color 0.1s',
    minWidth: 0,
    flexShrink: 0,
  } as React.CSSProperties,

  activeTab: {
    backgroundColor: '#1a1a2e',
    color: '#ddd',
    borderBottom: '2px solid #3b82f6',
  } as React.CSSProperties,

  tabName: {
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    maxWidth: '120px',
    fontFamily: 'monospace',
    fontSize: '11px',
  } as React.CSSProperties,

  renameInput: {
    width: '80px',
    padding: '1px 4px',
    fontSize: '11px',
    fontFamily: 'monospace',
    backgroundColor: '#2a2a4a',
    color: '#ddd',
    border: '1px solid #3b82f6',
    borderRadius: '2px',
    outline: 'none',
  } as React.CSSProperties,

  closeBtn: {
    background: 'none',
    border: 'none',
    color: '#555',
    cursor: 'pointer',
    fontSize: '13px',
    lineHeight: 1,
    padding: '0 2px',
    borderRadius: '2px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
  } as React.CSSProperties,

  addBtn: {
    padding: '0 10px',
    background: 'none',
    border: 'none',
    borderLeft: '1px solid #2a2a4a',
    color: '#666',
    cursor: 'pointer',
    fontSize: '16px',
    lineHeight: 1,
    flexShrink: 0,
    display: 'flex',
    alignItems: 'center',
  } as React.CSSProperties,
};
