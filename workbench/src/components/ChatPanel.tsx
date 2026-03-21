/**
 * Chat Panel — Geometry-aware agent chat for yapCAD workbench.
 * Streams responses from the OpenClaw gateway.
 *
 * Two modes:
 *   • Agent session (default): routes through jarvis-rich named session —
 *     full memory, workspace tools, and file access.
 *   • Stateless (legacy): raw completions, local history only.
 */

import { useState, useRef, useEffect, useCallback } from 'react';
import {
  useChat,
  AVAILABLE_MODELS,
  DEFAULT_MODEL,
  DEFAULT_USER,
  KNOWN_USERS,
  resolveUserSession,
} from '../hooks/useChat';
import CommandPalette from './CommandPalette';
import { SkillEditor } from './SkillEditor';

interface SkillParam { name: string; default: string; unit: string; }
interface Skill { name: string; description: string; category: string; promptTemplate: string; params: SkillParam[]; }

interface EvalContext {
  command: string;
  result: {
    success: boolean;
    error_message: string | null;
    volume: number | null;
    elapsed_ms: number;
  };
}

interface ChatPanelProps {
  dslSource?: string;
  onDSLUpdate?: (source: string) => void;
  evalContextRef?: React.MutableRefObject<EvalContext | null>;
}

const GATEWAY_URL_KEY   = 'yapcad-chat-gateway-url';
const GATEWAY_TOKEN_KEY = 'yapcad-chat-token';
const CHAT_MODEL_KEY    = 'yapcad-chat-model';
const CHAT_USER_KEY     = 'yapcad-chat-user';
const DEFAULT_GATEWAY   = '/openclaw';

// ── DSL block extraction ───────────────────────────────────────────────────────

function extractDSLBlocks(content: string): { before: string; code: string; after: string }[] {
  const blocks: { before: string; code: string; after: string }[] = [];
  const regex = /```dsl\n([\s\S]*?)```/g;
  let lastIndex = 0;
  let match;
  while ((match = regex.exec(content)) !== null) {
    blocks.push({ before: content.slice(lastIndex, match.index), code: match[1].trim(), after: '' });
    lastIndex = match.index + match[0].length;
  }
  if (blocks.length === 0) return [];
  blocks[blocks.length - 1].after = content.slice(lastIndex);
  return blocks;
}

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
                <button style={styles.applyBtn} onClick={() => onApply(block.code)} title="Apply to editor">
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

// ── Main panel ────────────────────────────────────────────────────────────────

export function ChatPanel({ dslSource = '', onDSLUpdate, evalContextRef }: ChatPanelProps) {
  const [gatewayUrl, setGatewayUrl] = useState(() =>
    localStorage.getItem(GATEWAY_URL_KEY) || DEFAULT_GATEWAY);
  const [token, setToken] = useState(() =>
    localStorage.getItem(GATEWAY_TOKEN_KEY) || '');
  const [model, setModel] = useState(() =>
    localStorage.getItem(CHAT_MODEL_KEY) || DEFAULT_MODEL);
  const [userName, setUserName] = useState(() =>
    localStorage.getItem(CHAT_USER_KEY) || DEFAULT_USER);

  // Derived from userName — these are never stored separately
  const { agentId, sessionKey } = resolveUserSession(userName);

  const [showSettings, setShowSettings] = useState(false);
  const [inputValue, setInputValue] = useState('');
  const messagesEndRef = useRef<HTMLDivElement>(null);

  // Image attachment
  const [imageAttachment, setImageAttachment] = useState<{
    preview: string; base64: string; mimeType: string;
  } | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Slash command palette
  const [showPalette, setShowPalette] = useState(false);
  const paletteFilter = inputValue.startsWith('/') ? inputValue.slice(1) : '';

  // Skill editor
  const [editingSkill, setEditingSkill] = useState<Skill | null>(null);

  const agentMode = !!(agentId && sessionKey);

  const { messages, isStreaming, error, send, cancel, clear } = useChat({
    gatewayUrl,
    token,
    dslSource,
    model,
    agentId,
    sessionKey,
  });

  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  const handleSaveSettings = useCallback((
    url: string, tok: string, mdl: string, user: string
  ) => {
    const cleanUrl = url.replace(/\/$/, '');
    setGatewayUrl(cleanUrl);  localStorage.setItem(GATEWAY_URL_KEY,   cleanUrl);
    setToken(tok);            localStorage.setItem(GATEWAY_TOKEN_KEY, tok);
    setModel(mdl);            localStorage.setItem(CHAT_MODEL_KEY,    mdl);
    setUserName(user);        localStorage.setItem(CHAT_USER_KEY,     user);
    setShowSettings(false);
  }, []);

  const handleFileSelect = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      const dataUrl = ev.target?.result as string;
      const base64 = dataUrl.split(',')[1];
      setImageAttachment({ preview: dataUrl, base64, mimeType: file.type });
    };
    reader.readAsDataURL(file);
    e.target.value = '';
  }, []);

  const handleSend = useCallback(async () => {
    const text = inputValue.trim();
    if (!text || isStreaming) return;
    setInputValue('');
    setShowPalette(false);
    const img = imageAttachment;
    setImageAttachment(null);

    // Prepend eval context if available (Task 4: eval result feedback to chat)
    let messageText = text;
    if (evalContextRef?.current) {
      const { command, result } = evalContextRef.current;
      if (result.success) {
        const volStr = result.volume != null ? ` | volume=${result.volume.toFixed(1)}mm³` : '';
        messageText = `[Eval: SUCCESS | cmd=${command}${volStr} | time=${result.elapsed_ms.toFixed(0)}ms]\n${text}`;
      } else {
        const errStr = result.error_message || 'unknown error';
        messageText = `[Eval: ERROR | cmd=${command} | ${errStr}]\n${text}`;
      }
      // Clear after consuming
      evalContextRef.current = null;
    }

    await send(messageText, img ?? undefined);
  }, [inputValue, isStreaming, send, imageAttachment, evalContextRef]);

  const handleKeyDown = useCallback((e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Escape') { setShowPalette(false); return; }
    if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); handleSend(); }
  }, [handleSend]);

  const handleCommandSelect = useCallback((cmd: string) => {
    setInputValue('');
    setShowPalette(false);
    if (cmd === '/new') { clear(); }
    else if (cmd === '/status') { send('/status'); }
    else if (cmd === '/compact') { send('/compact'); }
    else if (cmd === '/help') { send('/help'); }
  }, [clear, send]);

  const handleSkillSelect = useCallback((skill: Skill) => {
    setShowPalette(false);
    setInputValue('');
    setEditingSkill(skill);
  }, []);

  const handleSkillSubmit = useCallback((prompt: string) => {
    setEditingSkill(null);
    send(prompt);
  }, [send]);

  const handleApplyDSL = useCallback((code: string) => {
    onDSLUpdate?.(code);
  }, [onDSLUpdate]);

  if (showSettings) {
    return (
      <SettingsPanel
        initialUrl={gatewayUrl}
        initialToken={token}
        initialModel={model}
        initialUserName={userName}
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
        <span style={styles.title}>
          {agentMode ? 'Jarvis' : 'Agent Chat'}
        </span>
        {agentMode ? (
          <span style={styles.sessionBadge} title={`Agent: ${agentId} | Session: ${sessionKey}`}>
            🧠 {userName}
          </span>
        ) : (
          <span style={styles.modelBadge}>
            {AVAILABLE_MODELS.find(m => m.id === model)?.label.split(' ')[0] ?? 'Opus'}
          </span>
        )}
        <div style={styles.headerActions}>
          {messages.length > 0 && (
            <button style={styles.iconBtn} onClick={clear} title="Clear local display">✕</button>
          )}
          <button style={styles.iconBtn} onClick={() => setShowSettings(true)} title="Settings">⚙</button>
        </div>
      </div>

      {/* Session mode info bar */}
      {agentMode && isConfigured && (
        <div style={styles.sessionBar}>
          <span style={styles.sessionBarText}>
            ⚡ Agent session · {agentId} · full tools &amp; memory
          </span>
        </div>
      )}

      {/* Messages */}
      <div style={styles.messages}>
        {messages.length === 0 && (
          <div style={styles.emptyState}>
            {isConfigured ? (
              <>
                <div style={styles.emptyIcon}>{agentMode ? '🧠' : '💬'}</div>
                <div style={styles.emptyText}>
                  {agentMode
                    ? 'Connected to Jarvis agent session.'
                    : 'Ask me anything about your geometry.'}
                </div>
                <div style={styles.emptyHint}>
                  {agentMode
                    ? 'I have access to workspace files, memory, and tools. DSL context is sent with every message.'
                    : 'I can see your DSL and suggest improvements.'}
                </div>
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
              {msg.role === 'user' ? 'You' : (agentMode ? 'Jarvis' : 'Assistant')}
            </div>
            <div style={styles.messageContent}>
              {msg.role === 'assistant' ? (
                <>
                  <MessageContent
                    content={msg.content}
                    onApply={onDSLUpdate ? handleApplyDSL : undefined}
                  />
                  {msg.pending && <span style={styles.cursor}>▊</span>}
                </>
              ) : (
                <>
                  {msg.imagePreview && (
                    <img src={msg.imagePreview} style={styles.messageImage} alt="attachment" />
                  )}
                  <span style={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
                    {msg.content}
                  </span>
                </>
              )}
            </div>
          </div>
        ))}

        {error && <div style={styles.errorBanner}>{error}</div>}
        <div ref={messagesEndRef} />
      </div>

      {/* Hidden file input */}
      <input
        ref={fileInputRef}
        type="file"
        accept="image/png,image/jpeg,image/gif,image/webp"
        style={{ display: 'none' }}
        onChange={handleFileSelect}
      />

      {/* Skill editor modal */}
      {editingSkill && (
        <SkillEditor
          skill={editingSkill}
          onSubmit={handleSkillSubmit}
          onClose={() => setEditingSkill(null)}
        />
      )}

      {/* Input area */}
      <div style={{ position: 'relative', flexShrink: 0 }}>
        {/* Command palette */}
        {showPalette && (
          <CommandPalette
            filter={paletteFilter}
            onSelectCommand={handleCommandSelect}
            onSelectSkill={handleSkillSelect}
            onDismiss={() => setShowPalette(false)}
          />
        )}

        {/* Image preview strip */}
        {imageAttachment && (
          <div style={styles.imagePreview}>
            <img src={imageAttachment.preview} style={styles.imageThumb} alt="attachment" />
            <button style={styles.removeImageBtn} onClick={() => setImageAttachment(null)} title="Remove">✕</button>
          </div>
        )}

        <div style={styles.inputContainer}>
          <button
            style={styles.attachBtn}
            onClick={() => fileInputRef.current?.click()}
            title="Attach image"
            disabled={!isConfigured}
          >📎</button>
          <input
            style={styles.input(isConfigured)}
            type="text"
            value={inputValue}
            onChange={e => {
              setInputValue(e.target.value);
              setShowPalette(e.target.value.startsWith('/'));
            }}
            onKeyDown={handleKeyDown}
            placeholder={isConfigured
              ? (agentMode ? 'Message Jarvis... (/ for commands)' : 'Ask about this geometry...')
              : 'Configure token first'}
            disabled={!isConfigured}
            autoComplete="off"
          />
          {isStreaming ? (
            <button style={styles.stopBtn} onClick={cancel} title="Stop generation">◼ Stop</button>
          ) : (
            <button
              style={styles.sendBtn(isConfigured && !isStreaming && !!inputValue.trim())}
              onClick={handleSend}
              disabled={!isConfigured || isStreaming || !inputValue.trim()}
              title="Send (Enter)"
            >▶</button>
          )}
        </div>
      </div>
    </div>
  );
}

// ── Settings panel ────────────────────────────────────────────────────────────

function SettingsPanel({
  initialUrl, initialToken, initialModel, initialUserName,
  onSave, onCancel,
}: {
  initialUrl: string;
  initialToken: string;
  initialModel: string;
  initialUserName: string;
  onSave: (url: string, token: string, model: string, userName: string) => void;
  onCancel: () => void;
}) {
  const [url, setUrl]   = useState(initialUrl);
  const [tok, setTok]   = useState(initialToken);
  const [mdl, setMdl]   = useState(initialModel);
  const [user, setUser] = useState(initialUserName);

  const knownUsers = Object.keys(KNOWN_USERS);
  const isKnown    = knownUsers.includes(user.trim().toLowerCase());
  const { agentId, sessionKey } = resolveUserSession(user);

  return (
    <div style={styles.container}>
      <div style={styles.header}>
        <span style={styles.title}>Chat Settings</span>
      </div>
      <div style={styles.settingsBody}>

        <div style={styles.settingsSection}>IDENTITY</div>
        <label style={styles.settingsLabel}>Your Name</label>
        <div style={{ display: 'flex', gap: '6px' }}>
          <select
            style={{ ...styles.settingsInput, flex: 1 }}
            value={isKnown ? user.trim().toLowerCase() : '__custom'}
            onChange={e => { if (e.target.value !== '__custom') setUser(e.target.value); }}
          >
            {knownUsers.map(u => (
              <option key={u} value={u}>{u.charAt(0).toUpperCase() + u.slice(1)}</option>
            ))}
            {!isKnown && <option value="__custom">{user} (custom)</option>}
          </select>
          <input
            style={{ ...styles.settingsInput, flex: 1 }}
            value={user}
            onChange={e => setUser(e.target.value)}
            placeholder="your name"
          />
        </div>
        <div style={styles.settingsHint}>
          {isKnown
            ? `✅ Agent: ${agentId} → Session: ${sessionKey}`
            : `⚠️ Unknown user — will create session: ${sessionKey}`}
        </div>

        <div style={styles.settingsSection}>CONNECTION</div>
        <label style={styles.settingsLabel}>Model</label>
        <select style={styles.settingsInput} value={mdl} onChange={e => setMdl(e.target.value)}>
          {AVAILABLE_MODELS.map(m => <option key={m.id} value={m.id}>{m.label}</option>)}
        </select>
        <label style={styles.settingsLabel}>Gateway URL</label>
        <input style={styles.settingsInput} value={url} onChange={e => setUrl(e.target.value)}
          placeholder="/openclaw" />
        <label style={styles.settingsLabel}>Auth Token</label>
        <input style={styles.settingsInput} type="password" value={tok}
          onChange={e => setTok(e.target.value)} placeholder="OpenClaw gateway token" />
        <div style={styles.settingsHint}>
          Token: <code>openclaw status</code> → gateway auth token
        </div>

        <div style={styles.settingsBtns}>
          <button style={styles.cancelBtn} onClick={onCancel}>Cancel</button>
          <button style={styles.saveBtn} onClick={() => onSave(url, tok, mdl, user)}>Save</button>
        </div>
      </div>
    </div>
  );
}

// ── Styles ────────────────────────────────────────────────────────────────────

const styles = {
  container: { display: 'flex', flexDirection: 'column' as const, height: '100%',
    backgroundColor: '#13131f', overflow: 'hidden' },
  header: { padding: '8px 10px', backgroundColor: '#1e1e32',
    borderBottom: '1px solid #2a2a42', display: 'flex', alignItems: 'center',
    gap: '7px', flexShrink: 0 },
  title: { fontSize: '13px', color: '#ccc', fontWeight: 500, flex: 1 },
  modelBadge: { fontSize: '10px', color: '#7c7cf8', backgroundColor: '#1e1e38',
    padding: '2px 6px', borderRadius: '3px', fontWeight: 600, letterSpacing: '0.03em' },
  sessionBadge: { fontSize: '10px', color: '#4ade80', backgroundColor: '#0f2318',
    padding: '2px 7px', borderRadius: '3px', fontWeight: 600, letterSpacing: '0.02em',
    border: '1px solid #1a4a2a', cursor: 'default' },
  sessionBar: { padding: '4px 10px', backgroundColor: '#0f1a14',
    borderBottom: '1px solid #1a3a20', flexShrink: 0 },
  sessionBarText: { fontSize: '10px', color: '#4ade80', letterSpacing: '0.02em' },
  statusDot: (active: boolean) => ({ width: '7px', height: '7px', borderRadius: '50%',
    backgroundColor: active ? '#4ade80' : '#555', flexShrink: 0 }),
  headerActions: { display: 'flex', gap: '4px' },
  iconBtn: { background: 'none', border: 'none', color: '#777', cursor: 'pointer',
    fontSize: '13px', padding: '2px 5px', borderRadius: '3px', lineHeight: 1 } as React.CSSProperties,
  messages: { flex: 1, overflowY: 'auto' as const, padding: '8px',
    display: 'flex', flexDirection: 'column' as const, gap: '8px' },
  emptyState: { flex: 1, display: 'flex', flexDirection: 'column' as const,
    alignItems: 'center', justifyContent: 'center', padding: '24px 12px',
    textAlign: 'center' as const, color: '#555' },
  emptyIcon: { fontSize: '28px', marginBottom: '10px' },
  emptyText: { fontSize: '13px', color: '#666', marginBottom: '6px' },
  emptyHint: { fontSize: '11px', color: '#555', maxWidth: '220px', lineHeight: 1.5 },
  setupBtn: { marginTop: '14px', padding: '7px 16px', fontSize: '12px',
    backgroundColor: '#3b82f6', color: '#fff', border: 'none', borderRadius: '5px',
    cursor: 'pointer' } as React.CSSProperties,
  message: (role: 'user' | 'assistant') => ({ display: 'flex',
    flexDirection: 'column' as const, alignItems: role === 'user' ? 'flex-end' : 'flex-start' }),
  messageRole: (role: 'user' | 'assistant') => ({ fontSize: '10px',
    color: role === 'user' ? '#3b82f6' : '#7c7cf8', marginBottom: '3px', fontWeight: 600,
    letterSpacing: '0.04em', textTransform: 'uppercase' as const }),
  messageContent: { maxWidth: '92%', fontSize: '12px', color: '#ddd', lineHeight: 1.55 },
  cursor: { display: 'inline-block', animation: 'blink 1s step-end infinite',
    color: '#7c7cf8', marginLeft: '1px' },
  dslBlock: { backgroundColor: '#0e0e1c', border: '1px solid #2a2a48', borderRadius: '5px',
    marginTop: '6px', marginBottom: '6px', overflow: 'hidden' },
  dslHeader: { display: 'flex', alignItems: 'center', justifyContent: 'space-between',
    padding: '4px 8px', backgroundColor: '#1a1a30', borderBottom: '1px solid #2a2a48' },
  dslLabel: { fontSize: '10px', color: '#7c7cf8', fontWeight: 600, letterSpacing: '0.05em',
    textTransform: 'uppercase' as const },
  applyBtn: { padding: '2px 8px', fontSize: '11px', backgroundColor: '#4ade80', color: '#000',
    border: 'none', borderRadius: '3px', cursor: 'pointer', fontWeight: 600 } as React.CSSProperties,
  dslCode: { margin: 0, padding: '8px 10px', fontSize: '11px', color: '#a8d8a8',
    fontFamily: '"JetBrains Mono", "Fira Code", monospace', overflowX: 'auto' as const,
    whiteSpace: 'pre' as const },
  errorBanner: { backgroundColor: '#3f1818', color: '#f87171', fontSize: '11px',
    padding: '7px 10px', borderRadius: '4px', border: '1px solid #5a2020' },
  inputContainer: { display: 'flex', gap: '5px', padding: '8px',
    borderTop: '1px solid #2a2a42', backgroundColor: '#1a1a2e', flexShrink: 0 },
  input: (enabled: boolean) => ({ flex: 1, padding: '7px 10px', fontSize: '12px',
    backgroundColor: enabled ? '#23233a' : '#1a1a2e',
    border: '1px solid ' + (enabled ? '#3b3b5a' : '#2a2a42'), borderRadius: '5px',
    color: enabled ? '#ddd' : '#555', cursor: enabled ? 'text' : 'not-allowed',
    outline: 'none' } as React.CSSProperties),
  sendBtn: (active: boolean) => ({ padding: '7px 10px', fontSize: '13px',
    backgroundColor: active ? '#3b82f6' : '#2a2a42', color: active ? '#fff' : '#555',
    border: 'none', borderRadius: '5px', cursor: active ? 'pointer' : 'not-allowed',
    flexShrink: 0, transition: 'background-color 0.15s' } as React.CSSProperties),
  // Settings
  settingsBody: { padding: '14px', display: 'flex', flexDirection: 'column' as const, gap: '5px',
    overflowY: 'auto' as const, flex: 1 },
  settingsSection: { fontSize: '10px', color: '#555', fontWeight: 700, letterSpacing: '0.08em',
    textTransform: 'uppercase' as const, marginTop: '12px', marginBottom: '2px',
    paddingBottom: '4px', borderBottom: '1px solid #2a2a42' },
  settingsLabel: { fontSize: '11px', color: '#888', fontWeight: 500, marginTop: '4px' },
  settingsInput: { padding: '7px 10px', fontSize: '12px', backgroundColor: '#1e1e32',
    border: '1px solid #3b3b5a', borderRadius: '5px', color: '#ddd',
    outline: 'none' } as React.CSSProperties,
  settingsHint: { fontSize: '11px', color: '#555', lineHeight: 1.45 },
  settingsBtns: { display: 'flex', justifyContent: 'flex-end', gap: '8px', marginTop: '16px' },
  cancelBtn: { padding: '7px 14px', fontSize: '12px', backgroundColor: '#2a2a42', color: '#aaa',
    border: 'none', borderRadius: '5px', cursor: 'pointer' } as React.CSSProperties,
  saveBtn: { padding: '7px 14px', fontSize: '12px', backgroundColor: '#3b82f6', color: '#fff',
    border: 'none', borderRadius: '5px', cursor: 'pointer' } as React.CSSProperties,
  // New styles for image attach, stop button, command palette
  imagePreview: { display: 'flex', alignItems: 'center', gap: '6px', padding: '4px 8px',
    backgroundColor: '#1a1a2e', borderTop: '1px solid #2a2a42' } as React.CSSProperties,
  imageThumb: { width: '48px', height: '48px', objectFit: 'cover' as const,
    borderRadius: '4px', border: '1px solid #3b3b5a' } as React.CSSProperties,
  removeImageBtn: { background: 'none', border: 'none', color: '#888',
    cursor: 'pointer', fontSize: '14px', padding: '2px' } as React.CSSProperties,
  attachBtn: { padding: '7px 9px', fontSize: '13px', backgroundColor: '#23233a',
    border: '1px solid #3b3b5a', borderRadius: '5px', cursor: 'pointer',
    color: '#888', flexShrink: 0 } as React.CSSProperties,
  stopBtn: { padding: '7px 12px', fontSize: '12px', backgroundColor: '#7f1d1d',
    color: '#fca5a5', border: '1px solid #991b1b', borderRadius: '5px',
    cursor: 'pointer', flexShrink: 0, fontWeight: 600 } as React.CSSProperties,
  messageImage: { maxWidth: '200px', maxHeight: '150px', borderRadius: '4px',
    border: '1px solid #3b3b5a', marginBottom: '4px', display: 'block' } as React.CSSProperties,
};
