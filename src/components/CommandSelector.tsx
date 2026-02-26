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

      if (cmds.length > 0 && (!selectedCommand || !cmds.find(c => c.name === selectedCommand))) {
        setSelectedCommand(cmds[0].name);
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

  useEffect(() => {
    const timer = setTimeout(parseCommands, 800);
    return () => clearTimeout(timer);
  }, [parseCommands]);

  useEffect(() => {
    const cmd = commands.find(c => c.name === selectedCommand);
    if (cmd) {
      const defaults: Record<string, string> = {};
      for (const p of cmd.params) {
        defaults[p.name] = paramValues[p.name] ?? (p.default !== null && p.default !== undefined ? String(p.default) : '');
      }
      setParamValues(defaults);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedCommand, commands]);

  const handleEvaluate = () => {
    if (!selectedCommand || isEvaluating) return;
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
  const hasParams = selectedCmd && selectedCmd.params.length > 0;

  return (
    <div style={{
      display: 'flex',
      flexDirection: 'column',
      gap: '4px',
      padding: '6px 10px',
      background: '#1e1e3a',
      borderTop: '1px solid #2a2a50',
      flexShrink: 0,
    }}>

      {/* Row 1: Command selector + status */}
      <div style={{ display: 'flex', alignItems: 'center', gap: '6px', minHeight: '28px' }}>
        {commands.length > 0 ? (
          <select
            style={{
              flex: 1,
              padding: '4px 6px',
              background: '#2a2a4a',
              color: '#eee',
              border: '1px solid #444',
              borderRadius: '4px',
              fontSize: '12px',
              fontFamily: 'monospace',
              cursor: 'pointer',
              minWidth: 0,
            }}
            value={selectedCommand}
            onChange={e => setSelectedCommand(e.target.value)}
          >
            {commands.map(cmd => (
              /* Show only command name — params shown separately below */
              <option key={cmd.name} value={cmd.name}>
                {cmd.name}
              </option>
            ))}
          </select>
        ) : (
          <span style={{ flex: 1, fontSize: '11px', color: '#555', fontStyle: 'italic' }}>
            {isParsing ? '⏳ Parsing…' : source.trim() ? 'No commands found' : 'Load DSL to see commands'}
          </span>
        )}

        {/* Param count badge */}
        {selectedCmd && selectedCmd.params.length > 0 && (
          <span style={{ fontSize: '10px', color: '#666', whiteSpace: 'nowrap' }}>
            {selectedCmd.params.length} params
          </span>
        )}
      </div>

      {/* Row 2: Parameters — compact 2-column grid, scrollable if tall */}
      {hasParams && (
        <div style={{
          display: 'grid',
          gridTemplateColumns: '1fr 1fr',
          gap: '3px 8px',
          maxHeight: '110px',
          overflowY: 'auto',
          padding: '2px 0',
        }}>
          {selectedCmd.params.map(p => (
            <div key={p.name} style={{ display: 'flex', alignItems: 'center', gap: '3px', minWidth: 0 }}>
              <label
                title={`${p.name}: ${p.type}`}
                style={{
                  fontSize: '10px',
                  color: '#999',
                  fontFamily: 'monospace',
                  whiteSpace: 'nowrap',
                  overflow: 'hidden',
                  textOverflow: 'ellipsis',
                  flexShrink: 1,
                  minWidth: '30px',
                  cursor: 'default',
                }}
              >
                {p.name}
              </label>
              <input
                style={{
                  flex: 1,
                  minWidth: '40px',
                  maxWidth: '70px',
                  padding: '2px 4px',
                  background: '#2a2a4a',
                  color: '#eee',
                  border: '1px solid #444',
                  borderRadius: '3px',
                  fontSize: '11px',
                  fontFamily: 'monospace',
                }}
                value={paramValues[p.name] ?? ''}
                onChange={e => setParamValues(prev => ({ ...prev, [p.name]: e.target.value }))}
                placeholder={p.default !== null ? String(p.default) : ''}
                onKeyDown={e => { if (e.key === 'Enter') handleEvaluate(); }}
                title={`${p.name}: ${p.type}${p.default !== null ? ` = ${p.default}` : ''}`}
              />
            </div>
          ))}
        </div>
      )}

      {/* Row 3: Error message */}
      {parseError && (
        <div style={{ fontSize: '10px', color: '#ff6b6b', lineHeight: 1.3 }}>
          {parseError}
        </div>
      )}

      {/* Row 4: Evaluate button — full width, always at bottom */}
      <button
        style={{
          width: '100%',
          padding: '6px 0',
          background: isEvaluating || !selectedCommand ? '#1e3d2a' : '#2d7d46',
          color: isEvaluating || !selectedCommand ? '#5a8a6a' : '#fff',
          border: 'none',
          borderRadius: '4px',
          fontSize: '13px',
          fontWeight: 'bold',
          cursor: isEvaluating || !selectedCommand ? 'not-allowed' : 'pointer',
          letterSpacing: '0.3px',
          transition: 'background 0.15s',
          flexShrink: 0,
        }}
        onClick={handleEvaluate}
        disabled={isEvaluating || !selectedCommand}
        title="Ctrl+Enter"
      >
        {isEvaluating ? '⏳  Evaluating…' : '▶  Evaluate'}
      </button>

    </div>
  );
}
