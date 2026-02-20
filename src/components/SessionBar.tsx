/**
 * Session Bar Component - Shows connection status and backend configuration
 */

import { useState, useCallback, useEffect } from 'react';

interface SessionBarProps {
  backendUrl: string;
  onBackendUrlChange: (url: string) => void;
  isConnected: boolean;
  onTestConnection: () => Promise<void>;
}

const styles = {
  container: {
    display: 'flex',
    alignItems: 'center',
    padding: '8px 12px',
    backgroundColor: '#252540',
    borderBottom: '1px solid #333',
    gap: '12px',
    flexWrap: 'wrap' as const,
  },
  statusIndicator: {
    width: '8px',
    height: '8px',
    borderRadius: '50%',
    flexShrink: 0,
  },
  statusConnected: {
    backgroundColor: '#4caf50',
    boxShadow: '0 0 6px #4caf50',
  },
  statusDisconnected: {
    backgroundColor: '#f44336',
    boxShadow: '0 0 6px #f44336',
  },
  statusText: {
    fontSize: '12px',
    color: '#888',
    flexShrink: 0,
  },
  urlContainer: {
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    flex: 1,
    minWidth: '200px',
  },
  urlLabel: {
    fontSize: '12px',
    color: '#aaa',
    flexShrink: 0,
  },
  urlInput: {
    flex: 1,
    padding: '4px 8px',
    fontSize: '12px',
    backgroundColor: '#333',
    border: '1px solid #444',
    borderRadius: '3px',
    color: '#eee',
    outline: 'none',
    minWidth: '120px',
  },
  testButton: {
    padding: '4px 8px',
    fontSize: '11px',
    backgroundColor: '#4c9f38',
    color: '#fff',
    border: 'none',
    borderRadius: '3px',
    cursor: 'pointer',
    flexShrink: 0,
  },
  authStatus: {
    fontSize: '12px',
    color: '#888',
    flexShrink: 0,
  },
};

export function SessionBar({ backendUrl, onBackendUrlChange, isConnected, onTestConnection }: SessionBarProps) {
  const [localUrl, setLocalUrl] = useState(backendUrl);
  const [isEditing, setIsEditing] = useState(false);
  const [isTesting, setIsTesting] = useState(false);

  // Update local state when prop changes
  useEffect(() => {
    setLocalUrl(backendUrl);
  }, [backendUrl]);

  const handleUrlSubmit = useCallback(() => {
    if (localUrl !== backendUrl) {
      onBackendUrlChange(localUrl);
    }
    setIsEditing(false);
  }, [localUrl, backendUrl, onBackendUrlChange]);

  const handleKeyPress = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleUrlSubmit();
    } else if (e.key === 'Escape') {
      setLocalUrl(backendUrl);
      setIsEditing(false);
    }
  }, [handleUrlSubmit, backendUrl]);

  const handleTestConnection = useCallback(async () => {
    setIsTesting(true);
    try {
      await onTestConnection();
    } finally {
      setIsTesting(false);
    }
  }, [onTestConnection]);

  return (
    <div style={styles.container}>
      {/* Connection Status */}
      <div
        style={{
          ...styles.statusIndicator,
          ...(isConnected ? styles.statusConnected : styles.statusDisconnected),
        }}
        title={isConnected ? 'Connected' : 'Disconnected'}
      />
      
      <span style={styles.statusText}>
        {isConnected ? 'Connected' : 'Disconnected'}
      </span>
      
      {/* Backend URL Configuration */}
      <div style={styles.urlContainer}>
        <span style={styles.urlLabel}>Backend:</span>
        <input
          style={{
            ...styles.urlInput,
            ...(isEditing ? { borderColor: '#4c9f38' } : {}),
          }}
          type="text"
          value={localUrl}
          onChange={(e) => setLocalUrl(e.target.value)}
          onFocus={() => setIsEditing(true)}
          onBlur={handleUrlSubmit}
          onKeyDown={handleKeyPress}
          placeholder="http://localhost:8000"
        />
        <button
          style={styles.testButton}
          onClick={handleTestConnection}
          disabled={isTesting}
          title="Test connection"
        >
          {isTesting ? '⏳' : '🔄'}
        </button>
      </div>
      
      {/* Auth Status */}
      <span style={styles.authStatus}>
        Not authenticated
      </span>
    </div>
  );
}