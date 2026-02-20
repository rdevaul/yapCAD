/**
 * Chat Panel — Geometry-aware agent chat for yapCAD workbench.
 * Streams responses from the OpenClaw gateway.
 * Detects DSL code blocks and surfaces "Apply" buttons.
 */

import { useState, useRef, useEffect, useCallback } from 'react';
import { useChat } from '../hooks/useChat';

interface ChatPanelProps {
  dslSource?: string;
  onDSLUpdate?: (source: string) => void;
}

const GATEWAY_URL_KEY = 'yapcad-chat-gateway-url';
const GATEWAY_TOKEN_KEY = 'yapcad-chat-token';
const DEFAULT_GATEWAY = '/openclaw';

// Extract DSL code blocks from assistant message
function extractDSLBlocks(content: string): { before: string; code: string; after: string }[] {
  const blocks: { before: string; code: string; after: string }[] = [];
  const regex = /```dsl\n([\s\S]*?)```/g;
  let lastIndex = 0;
  let match;
  while ((match = regex.exec(content)) !== null) {
    blocks.push({
      before: content.slice(lastIndex, match.index),
      code: match[1].trim(),
      after: '',
    });
    lastIndex = match.index + match[0].length;
  }
  if (blocks.length === 0) return [];
  blocks[blocks.length - 1].after = content.slice(lastIndex);
  return blocks;
}

// Render text with DSL blocks surfaced as apply buttons
function MessageContent({
  content,
  onApply,
}: {
  content: string;
  onApply?: (code: string) => void;
}) {
  const blocks = extractDSLBlocks(content);
  if (blocks.length === 0) {
    return <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>{content}</span>;
  }

  return (
    <>
      {blocks.map((block, i) => (
        <span key={i}>
          <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>{block.before}</span>
          <div style={styles.dslBlock}>
            <div style={styles.dslHeader}>
              <span style={styles.dslLabel}>DSL</span>
              {onApply && (
                <button
                  style={styles.applyBtn}
                  onClick={() => onApply(block.code)}
                  title="Apply this DSL to the editor"
                >
                  ▶ Apply
                </button>
              )}
            </div>
            <pre style={styles.dslCode}>{block.code}</pre>
          </div>
          {i === blocks.length - 1 && (
            <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>{block.after}</span>
          )}
        </span>
      ))}
    </>
  );
}

export function ChatPanel({ dslSource = '', onDSLUpdate }: ChatPanelProps) {
  const [gatewayUrl, setGatewayUrl] = useState(() =>
    localStorage.getItem(GATEWAY_URL_KEY) || DEFAULT_GATEWAY
  );
  const [token, setToken] = useState(() =>
    localStorage.getItem(GATEWAY_TOKEN_KEY) || ''
  );
  const [showSettings, setShowSettings] = useState(false);
  const [inputValue, setInputValue] = useState('');
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const { messages, isStreaming, error, send, clear } = useChat({
    gatewayUrl,
    token,
    dslSource,
  });

  // Scroll to bottom on new content
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  const handleSaveSettings = useCallback((url: string, tok: string) => {
    const cleanUrl = url.replace(/\/$/, '');
    setGatewayUrl(cleanUrl);
    setToken(tok);
    localStorage.setItem(GATEWAY_URL_KEY, cleanUrl);
    localStorage.setItem(GATEWAY_TOKEN_KEY, tok);
    setShowSettings(false);
  }, []);

  const handleSend = useCallback(async () => {
    const text = inputValue.trim();
    if (!text || isStreaming) return;
    setInputValue('');
    await send(text);
  }, [inputValue, isStreaming, send]);

  const handleKeyDown = useCallback((e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  }, [handleSend]);

  const handleApplyDSL = useCallback((code: string) => {
    if (onDSLUpdate) {
      onDSLUpdate(code);
    }
  }, [onDSLUpdate]);

  if (showSettings) {
    return (
      <SettingsPanel
        initialUrl={gatewayUrl}
        initialToken={token}
        onSave={handleSaveSettings}
        onCancel={() => setShowSettings(false)}
      />
    );
  }

  const isConfigured = !!token;

  return (
    <div style={styles.container}>
      {/* Header */}
      <div style={styles.header}>
        <span style={styles.statusDot(isConfigured)} title={isConfigured ? 'Connected' : 'Not configured'} />
        <span style={styles.title}>Agent Chat</span>
        <div style={styles.headerActions}>
          {messages.length > 0 && (
            <button style={styles.iconBtn} onClick={clear} title="Clear conversation">
              ✕
            </button>
          )}
          <button style={styles.iconBtn} onClick={() => setShowSettings(true)} title="Settings">
            ⚙
          </button>
        </div>
      </div>

      {/* Messages */}
      <div style={styles.messages}>
        {messages.length === 0 && (
          <div style={styles.emptyState}>
            {isConfigured ? (
              <>
                <div style={styles.emptyIcon}>💬</div>
                <div style={styles.emptyText}>Ask me anything about your geometry.</div>
                <div style={styles.emptyHint}>I can see your DSL and suggest improvements.</div>
              </>
            ) : (
              <>
                <div style={styles.emptyIcon}>🔑</div>
                <div style={styles.emptyText}>Configure your OpenClaw token to start chatting.</div>
                <button style={styles.setupBtn} onClick={() => setShowSettings(true)}>
                  Open Settings
                </button>
              </>
            )}
          </div>
        )}

        {messages.map(msg => (
          <div key={msg.id} style={styles.message(msg.role)}>
            <div style={styles.messageRole(msg.role)}>
              {msg.role === 'user' ? 'You' : 'Jarvis'}
            </div>
            <div style={styles.messageContent}>
              {msg.role === 'assistant' ? (
                <>
                  <MessageContent
                    content={msg.content}
                    onApply={onDSLUpdate ? handleApplyDSL : undefined}
                  />
                  {msg.pending && msg.content === '' && (
                    <span style={styles.cursor}>▊</span>
                  )}
                  {msg.pending && msg.content !== '' && (
                    <span style={styles.cursor}>▊</span>
                  )}
                </>
              ) : (
                <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
                  {msg.content}
                </span>
              )}
            </div>
          </div>
        ))}

        {error && (
          <div style={styles.errorBanner}>{error}</div>
        )}

        <div ref={messagesEndRef} />
      </div>

      {/* Input */}
      <div style={styles.inputContainer}>
        <input
          ref={inputRef}
          style={styles.input(isConfigured)}
          type="text"
          value={inputValue}
          onChange={e => setInputValue(e.target.value)}
          onKeyDown={handleKeyDown}
          placeholder={isConfigured ? 'Ask about this geometry...' : 'Configure token first'}
          disabled={!isConfigured || isStreaming}
          autoComplete="off"
        />
        <button
          style={styles.sendBtn(isConfigured && !isStreaming && !!inputValue.trim())}
          onClick={handleSend}
          disabled={!isConfigured || isStreaming || !inputValue.trim()}
          title="Send (Enter)"
        >
          {isStreaming ? '⏸' : '▶'}
        </button>
      </div>
    </div>
  );
}

// ─── Settings panel ───────────────────────────────────────────────────────────

function SettingsPanel({
  initialUrl,
  initialToken,
  onSave,
  onCancel,
}: {
  initialUrl: string;
  initialToken: string;
  onSave: (url: string, token: string) => void;
  onCancel: () => void;
}) {
  const [url, setUrl] = useState(initialUrl);
  const [tok, setTok] = useState(initialToken);

  return (
    <div style={styles.container}>
      <div style={styles.header}>
        <span style={styles.title}>Chat Settings</span>
      </div>
      <div style={styles.settingsBody}>
        <label style={styles.settingsLabel}>Gateway URL</label>
        <input
          style={styles.settingsInput}
          value={url}
          onChange={e => setUrl(e.target.value)}
          placeholder="http://localhost:18789"
        />
        <label style={styles.settingsLabel}>Auth Token</label>
        <input
          style={styles.settingsInput}
          type="password"
          value={tok}
          onChange={e => setTok(e.target.value)}
          placeholder="OpenClaw gateway token"
        />
        <div style={styles.settingsHint}>
          Find your token: <code>openclaw status</code> → gateway auth token
        </div>
        <div style={styles.settingsBtns}>
          <button style={styles.cancelBtn} onClick={onCancel}>Cancel</button>
          <button style={styles.saveBtn} onClick={() => onSave(url, tok)}>Save</button>
        </div>
      </div>
    </div>
  );
}

// ─── Styles ───────────────────────────────────────────────────────────────────

const styles = {
  container: {
    display: 'flex',
    flexDirection: 'column' as const,
    height: '100%',
    backgroundColor: '#13131f',
    overflow: 'hidden',
  },
  header: {
    padding: '8px 10px',
    backgroundColor: '#1e1e32',
    borderBottom: '1px solid #2a2a42',
    display: 'flex',
    alignItems: 'center',
    gap: '7px',
    flexShrink: 0,
  },
  title: {
    fontSize: '13px',
    color: '#ccc',
    fontWeight: 500,
    flex: 1,
  },
  statusDot: (active: boolean) => ({
    width: '7px',
    height: '7px',
    borderRadius: '50%',
    backgroundColor: active ? '#4ade80' : '#555',
    flexShrink: 0,
  }),
  headerActions: {
    display: 'flex',
    gap: '4px',
  },
  iconBtn: {
    background: 'none',
    border: 'none',
    color: '#777',
    cursor: 'pointer',
    fontSize: '13px',
    padding: '2px 5px',
    borderRadius: '3px',
    lineHeight: 1,
  } as React.CSSProperties,
  messages: {
    flex: 1,
    overflowY: 'auto' as const,
    padding: '8px',
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '8px',
  },
  emptyState: {
    flex: 1,
    display: 'flex',
    flexDirection: 'column' as const,
    alignItems: 'center',
    justifyContent: 'center',
    padding: '24px 12px',
    textAlign: 'center' as const,
    color: '#555',
  },
  emptyIcon: {
    fontSize: '28px',
    marginBottom: '10px',
  },
  emptyText: {
    fontSize: '13px',
    color: '#666',
    marginBottom: '6px',
  },
  emptyHint: {
    fontSize: '11px',
    color: '#555',
  },
  setupBtn: {
    marginTop: '14px',
    padding: '7px 16px',
    fontSize: '12px',
    backgroundColor: '#3b82f6',
    color: '#fff',
    border: 'none',
    borderRadius: '5px',
    cursor: 'pointer',
  } as React.CSSProperties,
  message: (role: 'user' | 'assistant') => ({
    display: 'flex',
    flexDirection: 'column' as const,
    alignItems: role === 'user' ? 'flex-end' : 'flex-start',
  }),
  messageRole: (role: 'user' | 'assistant') => ({
    fontSize: '10px',
    color: role === 'user' ? '#3b82f6' : '#7c7cf8',
    marginBottom: '3px',
    fontWeight: 600,
    letterSpacing: '0.04em',
    textTransform: 'uppercase' as const,
  }),
  messageContent: {
    maxWidth: '92%',
    fontSize: '12px',
    color: '#ddd',
    lineHeight: 1.55,
  },
  cursor: {
    display: 'inline-block',
    animation: 'blink 1s step-end infinite',
    color: '#7c7cf8',
    marginLeft: '1px',
  },
  dslBlock: {
    backgroundColor: '#0e0e1c',
    border: '1px solid #2a2a48',
    borderRadius: '5px',
    marginTop: '6px',
    marginBottom: '6px',
    overflow: 'hidden',
  },
  dslHeader: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    padding: '4px 8px',
    backgroundColor: '#1a1a30',
    borderBottom: '1px solid #2a2a48',
  },
  dslLabel: {
    fontSize: '10px',
    color: '#7c7cf8',
    fontWeight: 600,
    letterSpacing: '0.05em',
    textTransform: 'uppercase' as const,
  },
  applyBtn: {
    padding: '2px 8px',
    fontSize: '11px',
    backgroundColor: '#4ade80',
    color: '#000',
    border: 'none',
    borderRadius: '3px',
    cursor: 'pointer',
    fontWeight: 600,
  } as React.CSSProperties,
  dslCode: {
    margin: 0,
    padding: '8px 10px',
    fontSize: '11px',
    color: '#a8d8a8',
    fontFamily: '"JetBrains Mono", "Fira Code", monospace',
    overflowX: 'auto' as const,
    whiteSpace: 'pre' as const,
  },
  errorBanner: {
    backgroundColor: '#3f1818',
    color: '#f87171',
    fontSize: '11px',
    padding: '7px 10px',
    borderRadius: '4px',
    border: '1px solid #5a2020',
  },
  inputContainer: {
    display: 'flex',
    gap: '5px',
    padding: '8px',
    borderTop: '1px solid #2a2a42',
    backgroundColor: '#1a1a2e',
    flexShrink: 0,
  },
  input: (enabled: boolean) => ({
    flex: 1,
    padding: '7px 10px',
    fontSize: '12px',
    backgroundColor: enabled ? '#23233a' : '#1a1a2e',
    border: '1px solid ' + (enabled ? '#3b3b5a' : '#2a2a42'),
    borderRadius: '5px',
    color: enabled ? '#ddd' : '#555',
    cursor: enabled ? 'text' : 'not-allowed',
    outline: 'none',
  } as React.CSSProperties),
  sendBtn: (active: boolean) => ({
    padding: '7px 10px',
    fontSize: '13px',
    backgroundColor: active ? '#3b82f6' : '#2a2a42',
    color: active ? '#fff' : '#555',
    border: 'none',
    borderRadius: '5px',
    cursor: active ? 'pointer' : 'not-allowed',
    flexShrink: 0,
    transition: 'background-color 0.15s',
  } as React.CSSProperties),
  // Settings
  settingsBody: {
    padding: '14px',
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '6px',
  },
  settingsLabel: {
    fontSize: '11px',
    color: '#888',
    fontWeight: 500,
    marginTop: '6px',
  },
  settingsInput: {
    padding: '7px 10px',
    fontSize: '12px',
    backgroundColor: '#1e1e32',
    border: '1px solid #3b3b5a',
    borderRadius: '5px',
    color: '#ddd',
    outline: 'none',
  } as React.CSSProperties,
  settingsHint: {
    fontSize: '11px',
    color: '#555',
    marginTop: '4px',
  },
  settingsBtns: {
    display: 'flex',
    justifyContent: 'flex-end',
    gap: '8px',
    marginTop: '14px',
  },
  cancelBtn: {
    padding: '7px 14px',
    fontSize: '12px',
    backgroundColor: '#2a2a42',
    color: '#aaa',
    border: 'none',
    borderRadius: '5px',
    cursor: 'pointer',
  } as React.CSSProperties,
  saveBtn: {
    padding: '7px 14px',
    fontSize: '12px',
    backgroundColor: '#3b82f6',
    color: '#fff',
    border: 'none',
    borderRadius: '5px',
    cursor: 'pointer',
  } as React.CSSProperties,
};
