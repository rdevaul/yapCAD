import { useState, useEffect, useCallback } from 'react';

interface CommandParam {
  name: string;
  type: string;
  default: number | string | boolean | null;
}

interface CommandInfo {
  name: string;
  params: CommandParam[];
}

interface CommandSelectorProps {
  source: string;
  backendUrl: string;
  onEvaluate: (command: string, parameters: Record<string, unknown>) => Promise<void>;
  isEvaluating: boolean;
}

export function CommandSelector({ source, backendUrl, onEvaluate, isEvaluating }: CommandSelectorProps) {
  const [commands, setCommands] = useState<CommandInfo[]>([]);
  const [selectedCommand, setSelectedCommand] = useState<string>('');
  const [paramValues, setParamValues] = useState<Record<string, string>>({});
  const [parseError, setParseError] = useState<string>('');
  const [isParsing, setIsParsing] = useState(false);

  // Parse commands from source (debounced)
  const parseCommands = useCallback(async () => {
    if (!source.trim() || !backendUrl) return;
    setIsParsing(true);
    setParseError('');
    try {
      const response = await fetch(`${backendUrl}/dsl/commands`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ source }),
        signal: typeof AbortSignal.timeout === 'function' ? AbortSignal.timeout(5000) : undefined,
      });
      if (!response.ok) {
        const err = await response.text();
        setParseError(err);
        return;
      }
      const data = await response.json();
      const cmds: CommandInfo[] = data.commands || [];
      setCommands(cmds);

      // Auto-select first command if none selected or selection no longer valid
      if (cmds.length > 0 && (!selectedCommand || !cmds.find(c => c.name === selectedCommand))) {
        setSelectedCommand(cmds[0].name);
        // Initialize param values with defaults
        const defaults: Record<string, string> = {};
        for (const p of cmds[0].params) {
          defaults[p.name] = p.default !== null && p.default !== undefined ? String(p.default) : '';
        }
        setParamValues(defaults);
      }
    } catch (err) {
      console.error('CommandSelector parse error:', err);
      setParseError(err instanceof Error ? err.message : 'Failed to connect to backend');
    } finally {
      setIsParsing(false);
    }
  }, [source, backendUrl, selectedCommand]);

  // Debounce parsing: 800ms after source changes
  useEffect(() => {
    const timer = setTimeout(parseCommands, 800);
    return () => clearTimeout(timer);
  }, [parseCommands]);

  // When selected command changes, update param defaults
  useEffect(() => {
    const cmd = commands.find(c => c.name === selectedCommand);
    if (cmd) {
      const defaults: Record<string, string> = {};
      for (const p of cmd.params) {
        // Keep existing value if user already typed something, else use default
        defaults[p.name] = paramValues[p.name] ?? (p.default !== null && p.default !== undefined ? String(p.default) : '');
      }
      setParamValues(defaults);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedCommand, commands]);

  const handleEvaluate = () => {
    if (!selectedCommand || isEvaluating) return;
    // Convert param strings to typed values
    const cmd = commands.find(c => c.name === selectedCommand);
    const params: Record<string, unknown> = {};
    if (cmd) {
      for (const p of cmd.params) {
        const raw = paramValues[p.name];
        if (raw === undefined || raw === '') continue;
        if (p.type === 'float') params[p.name] = parseFloat(raw);
        else if (p.type === 'int') params[p.name] = parseInt(raw, 10);
        else if (p.type === 'bool') params[p.name] = raw === 'true';
        else params[p.name] = raw;
      }
    }
    onEvaluate(selectedCommand, params);
  };

  const selectedCmd = commands.find(c => c.name === selectedCommand);

  const styles = {
    container: {
      display: 'flex',
      flexDirection: 'column' as const,
      gap: '6px',
      padding: '8px 12px',
      background: '#1e1e3a',
      borderTop: '1px solid #333',
    },
    topRow: {
      display: 'flex',
      alignItems: 'center',
      gap: '8px',
    },
    select: {
      flex: 1,
      padding: '6px 8px',
      background: '#2a2a4a',
      color: '#eee',
      border: '1px solid #444',
      borderRadius: '4px',
      fontSize: '13px',
      fontFamily: 'monospace',
      cursor: 'pointer',
    },
    evalButton: {
      padding: '6px 16px',
      background: '#2d7d46',
      color: '#fff',
      border: 'none',
      borderRadius: '4px',
      fontSize: '13px',
      fontWeight: 'bold' as const,
      cursor: 'pointer',
      whiteSpace: 'nowrap' as const,
      opacity: 1,
    },
    evalButtonDisabled: {
      opacity: 0.5,
      cursor: 'not-allowed',
    },
    paramsRow: {
      display: 'flex',
      flexWrap: 'wrap' as const,
      gap: '6px',
      alignItems: 'center',
    },
    paramGroup: {
      display: 'flex',
      alignItems: 'center',
      gap: '3px',
    },
    paramLabel: {
      fontSize: '11px',
      color: '#aaa',
      fontFamily: 'monospace',
    },
    paramInput: {
      width: '60px',
      padding: '3px 6px',
      background: '#2a2a4a',
      color: '#eee',
      border: '1px solid #444',
      borderRadius: '3px',
      fontSize: '12px',
      fontFamily: 'monospace',
    },
    paramType: {
      fontSize: '10px',
      color: '#666',
      fontFamily: 'monospace',
    },
    noCommands: {
      fontSize: '12px',
      color: '#666',
      fontStyle: 'italic' as const,
    },
    error: {
      fontSize: '11px',
      color: '#ff6b6b',
    },
    parsing: {
      fontSize: '11px',
      color: '#666',
    },
  };

  return (
    <div style={styles.container}>
      <div style={styles.topRow}>
        {commands.length > 0 ? (
          <select
            style={styles.select}
            value={selectedCommand}
            onChange={e => setSelectedCommand(e.target.value)}
          >
            {commands.map(cmd => (
              <option key={cmd.name} value={cmd.name}>
                {cmd.name}({cmd.params.map(p => p.name).join(', ')})
              </option>
            ))}
          </select>
        ) : (
          <span style={styles.noCommands}>
            {isParsing ? '⏳ Parsing...' : source.trim() ? 'No commands found' : 'Write DSL above'}
          </span>
        )}
        <button
          style={{
            ...styles.evalButton,
            ...(isEvaluating || !selectedCommand ? styles.evalButtonDisabled : {}),
          }}
          onClick={handleEvaluate}
          disabled={isEvaluating || !selectedCommand}
          title="Ctrl+Enter"
        >
          {isEvaluating ? '⏳' : '▶'} Evaluate
        </button>
      </div>

      {selectedCmd && selectedCmd.params.length > 0 && (
        <div style={styles.paramsRow}>
          {selectedCmd.params.map(p => (
            <div key={p.name} style={styles.paramGroup}>
              <span style={styles.paramLabel}>{p.name}</span>
              <span style={styles.paramType}>:{p.type}</span>
              <input
                style={styles.paramInput}
                value={paramValues[p.name] ?? ''}
                onChange={e => setParamValues(prev => ({ ...prev, [p.name]: e.target.value }))}
                placeholder={p.default !== null ? String(p.default) : ''}
                onKeyDown={e => { if (e.key === 'Enter') handleEvaluate(); }}
              />
            </div>
          ))}
        </div>
      )}

      {parseError && <div style={styles.error}>{parseError}</div>}
    </div>
  );
}
