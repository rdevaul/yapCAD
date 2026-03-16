/**
 * ParameterPanel — Live parameter editing for the selected command.
 *
 * Level 1 type-driven widgets:
 *   float → number input with optional slider (if min/max from ui_hint)
 *   int   → number input (integer step)
 *   bool  → checkbox
 *   string → text input
 *   list[*] → textarea (JSON array)
 *
 * Triggers debounced eval on parameter change, shows eval status.
 */

import { useCallback, useRef, useEffect } from 'react';
import type { CommandInfo, CommandParam } from '../hooks/useTabState';
import type { EvalResult, LatencyTier } from '../hooks/useWsEval';

interface ParameterPanelProps {
  commands: CommandInfo[];
  selectedCommand: string;
  paramValues: Record<string, Record<string, unknown>>;
  onSelectCommand: (command: string) => void;
  onParamChange: (command: string, param: string, value: unknown) => void;
  onEval: (command: string, parameters: Record<string, unknown>) => void;
  onSetDefaults?: (command: string, params: Record<string, unknown>) => void;
  isEvaluating: boolean;
  evalResult: EvalResult | null;
  latencyTier: LatencyTier;
  parseError: string;
  isParsing: boolean;
}

// ── Latency tier display ──────────────────────────────────────────────────────

function LatencyBadge({ tier }: { tier: LatencyTier }) {
  const config = {
    fast: { icon: '⚡', label: 'live', color: '#4ade80', bg: '#0f2318' },
    medium: { icon: '⏱', label: 'fast', color: '#fbbf24', bg: '#1a1700' },
    slow: { icon: '🔄', label: 'slow', color: '#f87171', bg: '#1a0f0f' },
  }[tier];

  return (
    <span style={{
      fontSize: '9px',
      padding: '1px 5px',
      borderRadius: '3px',
      backgroundColor: config.bg,
      color: config.color,
      fontWeight: 600,
      letterSpacing: '0.03em',
    }}>
      {config.icon} {config.label}
    </span>
  );
}

// ── Eval status display ───────────────────────────────────────────────────────

function EvalStatus({ isEvaluating, result }: { isEvaluating: boolean; result: EvalResult | null }) {
  if (isEvaluating) {
    return (
      <div style={styles.evalStatus}>
        <span style={{ color: '#fbbf24', fontSize: '11px' }}>⏳ Evaluating...</span>
      </div>
    );
  }

  if (!result) return null;

  if (result.success) {
    return (
      <div style={{ ...styles.evalStatus, borderColor: '#1a4a2a' }}>
        <span style={{ color: '#4ade80', fontSize: '11px' }}>
          ✓ {result.volume != null ? `volume: ${result.volume.toFixed(1)}mm³` : 'OK'}
        </span>
        <span style={{ color: '#555', fontSize: '10px', marginLeft: '6px' }}>
          {result.elapsed_ms.toFixed(0)}ms
        </span>
      </div>
    );
  }

  return (
    <div style={{ ...styles.evalStatus, borderColor: '#4a1a1a' }}>
      <span style={{ color: '#f87171', fontSize: '11px' }}>
        ✗ {result.error_message || 'Eval failed'}
      </span>
    </div>
  );
}

// ── Parameter widget for a single param ───────────────────────────────────────

function ParamWidget({
  param,
  value,
  onChange,
}: {
  param: CommandParam;
  value: unknown;
  onChange: (value: unknown) => void;
}) {
  const type = param.type.toLowerCase();

  // point2d → (x, y) coordinate pair inputs
  if (type === 'point2d' || param.ui_hint?.widget === 'point2d') {
    // Value is [x, y, z, w] homogeneous or [x, y] — normalise to [x, y]
    const arr = Array.isArray(value) ? value as number[] :
                Array.isArray(param.default) ? param.default as number[] : [0, 0];
    const px = Number(arr[0] ?? 0);
    const py = Number(arr[1] ?? 0);
    const label = (param.ui_hint as Record<string, unknown>)?.label as string ?? param.name;

    return (
      <div style={{ ...styles.paramRow, alignItems: 'flex-start', flexDirection: 'column', gap: '4px' }}>
        <label style={{ ...styles.paramLabel, color: '#ff8c00', fontWeight: 600 }}>
          {label}
        </label>
        <div style={{ display: 'flex', gap: '6px', paddingLeft: '8px' }}>
          <label style={{ ...styles.paramLabel, color: '#888', fontSize: '10px', width: 'auto' }}>x</label>
          <input
            type="number"
            value={px}
            step={0.5}
            onChange={e => onChange([parseFloat(e.target.value) || 0, py, 0, 1])}
            style={{ ...styles.numberInput, width: '60px' }}
          />
          <label style={{ ...styles.paramLabel, color: '#888', fontSize: '10px', width: 'auto' }}>y</label>
          <input
            type="number"
            value={py}
            step={0.5}
            onChange={e => onChange([px, parseFloat(e.target.value) || 0, 0, 1])}
            style={{ ...styles.numberInput, width: '60px' }}
          />
        </div>
      </div>
    );
  }

  // Bool → checkbox
  if (type === 'bool') {
    return (
      <div style={styles.paramRow}>
        <label style={styles.paramLabel} title={`${param.name}: ${param.type}`}>
          {param.name}
        </label>
        <input
          type="checkbox"
          checked={Boolean(value ?? param.default ?? false)}
          onChange={e => onChange(e.target.checked)}
          style={{ cursor: 'pointer' }}
        />
      </div>
    );
  }

  // Float → number input + optional slider
  if (type === 'float') {
    const numValue = Number(value ?? param.default ?? 0);
    const hasRange = param.ui_hint?.min != null && param.ui_hint?.max != null;

    return (
      <div style={styles.paramRow}>
        <label style={styles.paramLabel} title={`${param.name}: ${param.type}`}>
          {param.name}
        </label>
        <div style={{ display: 'flex', alignItems: 'center', gap: '4px', flex: 1, minWidth: 0 }}>
          {hasRange && (
            <input
              type="range"
              min={param.ui_hint!.min}
              max={param.ui_hint!.max}
              step={param.ui_hint?.step ?? 0.1}
              value={numValue}
              onChange={e => onChange(parseFloat(e.target.value))}
              style={{ flex: 1, minWidth: '40px', cursor: 'pointer' }}
            />
          )}
          <input
            type="number"
            value={numValue}
            step={param.ui_hint?.step ?? 0.1}
            onChange={e => onChange(parseFloat(e.target.value) || 0)}
            style={{
              ...styles.numberInput,
              width: hasRange ? '55px' : '70px',
            }}
          />
        </div>
      </div>
    );
  }

  // Int → number input (integer step)
  if (type === 'int') {
    const intValue = Number(value ?? param.default ?? 0);

    return (
      <div style={styles.paramRow}>
        <label style={styles.paramLabel} title={`${param.name}: ${param.type}`}>
          {param.name}
        </label>
        <input
          type="number"
          value={intValue}
          step={1}
          onChange={e => onChange(parseInt(e.target.value, 10) || 0)}
          style={styles.numberInput}
        />
      </div>
    );
  }

  // list[*] → textarea (JSON array)
  if (type.startsWith('list')) {
    const strValue = typeof value === 'string'
      ? value
      : value != null ? JSON.stringify(value) : '';

    return (
      <div style={styles.paramRow}>
        <label style={styles.paramLabel} title={`${param.name}: ${param.type}`}>
          {param.name}
        </label>
        <textarea
          value={strValue}
          onChange={e => {
            try {
              const parsed = JSON.parse(e.target.value);
              onChange(parsed);
            } catch {
              onChange(e.target.value); // Store raw string on invalid JSON
            }
          }}
          placeholder="[1, 2, 3]"
          style={styles.textareaInput}
          rows={2}
        />
      </div>
    );
  }

  // Default: string → text input
  return (
    <div style={styles.paramRow}>
      <label style={styles.paramLabel} title={`${param.name}: ${param.type}`}>
        {param.name}
      </label>
      <input
        type="text"
        value={String(value ?? param.default ?? '')}
        onChange={e => onChange(e.target.value)}
        style={styles.textInput}
      />
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────────

export function ParameterPanel({
  commands,
  selectedCommand,
  paramValues,
  onSelectCommand,
  onParamChange,
  onEval,
  onSetDefaults,
  isEvaluating,
  evalResult,
  latencyTier,
  parseError,
  isParsing,
}: ParameterPanelProps) {
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const selectedCmd = commands.find(c => c.name === selectedCommand);
  const hasParams = selectedCmd && selectedCmd.params.length > 0;
  const currentParamValues = paramValues[selectedCommand] ?? {};

  // Debounced eval trigger
  const triggerEval = useCallback((command: string, params: Record<string, unknown>) => {
    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => {
      // Convert param values to typed format
      const typedParams: Record<string, unknown> = {};
      const cmd = commands.find(c => c.name === command);
      if (cmd) {
        for (const p of cmd.params) {
          const val = params[p.name];
          if (val === undefined || val === null) continue;
          typedParams[p.name] = val;
        }
      }
      onEval(command, typedParams);
    }, 300);
  }, [commands, onEval]);

  // Clean up debounce on unmount
  useEffect(() => {
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, []);

  const handleParamChange = useCallback((paramName: string, value: unknown) => {
    onParamChange(selectedCommand, paramName, value);
    // Trigger debounced eval with updated values
    const updatedParams = { ...currentParamValues, [paramName]: value };
    triggerEval(selectedCommand, updatedParams);
  }, [selectedCommand, currentParamValues, onParamChange, triggerEval]);

  const handleManualEval = useCallback(() => {
    if (!selectedCommand || isEvaluating) return;
    const typedParams: Record<string, unknown> = {};
    const cmd = commands.find(c => c.name === selectedCommand);
    if (cmd) {
      for (const p of cmd.params) {
        const val = currentParamValues[p.name];
        if (val === undefined || val === null) continue;
        typedParams[p.name] = val;
      }
    }
    onEval(selectedCommand, typedParams);
  }, [selectedCommand, isEvaluating, commands, currentParamValues, onEval]);

  return (
    <div style={styles.container}>

      {/* Row 1: Command selector + latency badge */}
      <div style={styles.headerRow}>
        {commands.length > 0 ? (
          <select
            style={styles.commandSelect}
            value={selectedCommand}
            onChange={e => onSelectCommand(e.target.value)}
          >
            {commands.map(cmd => (
              <option key={cmd.name} value={cmd.name}>
                {cmd.name}
              </option>
            ))}
          </select>
        ) : (
          <span style={styles.emptyHint}>
            {isParsing ? '⏳ Parsing…' : 'No commands found'}
          </span>
        )}

        {selectedCmd && <LatencyBadge tier={latencyTier} />}

        {selectedCmd && selectedCmd.params.length > 0 && (
          <span style={styles.paramCount}>
            {selectedCmd.params.length}p
          </span>
        )}
      </div>

      {/* Row 2: Parameters */}
      {hasParams && (
        <div style={styles.paramsContainer}>
          {selectedCmd.params.map(p => (
            <ParamWidget
              key={p.name}
              param={p}
              value={currentParamValues[p.name]}
              onChange={(val) => handleParamChange(p.name, val)}
            />
          ))}
        </div>
      )}

      {/* Row 3: Eval status */}
      <EvalStatus isEvaluating={isEvaluating} result={evalResult} />

      {/* Row 4: Parse error */}
      {parseError && (
        <div style={styles.errorLine}>
          {parseError}
        </div>
      )}

      {/* Row 5: Evaluate button + Set as defaults */}
      <div style={{ display: 'flex', gap: '6px' }}>
        <button
          style={{
            ...styles.evalButton,
            flex: 1,
            ...(isEvaluating || !selectedCommand ? styles.evalButtonDisabled : {}),
          }}
          onClick={handleManualEval}
          disabled={isEvaluating || !selectedCommand}
          title="Ctrl+Enter"
        >
          {isEvaluating ? '⏳  Evaluating…' : '▶  Evaluate'}
        </button>
        {onSetDefaults && hasParams && (
          <button
            style={{
              padding: '8px 10px',
              fontSize: '12px',
              backgroundColor: '#1a2a1a',
              color: '#4ade80',
              border: '1px solid #2a4a2a',
              borderRadius: '4px',
              cursor: 'pointer',
              whiteSpace: 'nowrap',
              flexShrink: 0,
            }}
            onClick={() => onSetDefaults(selectedCommand, currentParamValues)}
            title="Bake current parameter values into the DSL source as new defaults"
          >
            ✦ Set defaults
          </button>
        )}
      </div>
    </div>
  );
}

// ── Styles ────────────────────────────────────────────────────────────────────

const styles = {
  container: {
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '4px',
    padding: '6px 10px',
    background: '#1e1e3a',
    borderTop: '1px solid #2a2a50',
    flexShrink: 0,
  },

  headerRow: {
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    minHeight: '28px',
  },

  commandSelect: {
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
  } as React.CSSProperties,

  emptyHint: {
    flex: 1,
    fontSize: '11px',
    color: '#555',
    fontStyle: 'italic' as const,
  },

  paramCount: {
    fontSize: '10px',
    color: '#666',
    whiteSpace: 'nowrap' as const,
  },

  paramsContainer: {
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '4px',
    maxHeight: '160px',
    overflowY: 'auto' as const,
    padding: '2px 0',
  },

  paramRow: {
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    minWidth: 0,
  },

  paramLabel: {
    fontSize: '10px',
    color: '#999',
    fontFamily: 'monospace',
    whiteSpace: 'nowrap' as const,
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    minWidth: '50px',
    maxWidth: '80px',
    cursor: 'default',
  } as React.CSSProperties,

  numberInput: {
    width: '70px',
    padding: '3px 5px',
    background: '#2a2a4a',
    color: '#eee',
    border: '1px solid #444',
    borderRadius: '3px',
    fontSize: '11px',
    fontFamily: 'monospace',
  } as React.CSSProperties,

  textInput: {
    flex: 1,
    minWidth: '40px',
    padding: '3px 5px',
    background: '#2a2a4a',
    color: '#eee',
    border: '1px solid #444',
    borderRadius: '3px',
    fontSize: '11px',
    fontFamily: 'monospace',
  } as React.CSSProperties,

  textareaInput: {
    flex: 1,
    minWidth: '40px',
    padding: '3px 5px',
    background: '#2a2a4a',
    color: '#eee',
    border: '1px solid #444',
    borderRadius: '3px',
    fontSize: '10px',
    fontFamily: 'monospace',
    resize: 'vertical' as const,
  } as React.CSSProperties,

  evalStatus: {
    padding: '3px 8px',
    borderRadius: '3px',
    border: '1px solid #2a2a4a',
    backgroundColor: '#111828',
  },

  errorLine: {
    fontSize: '10px',
    color: '#ff6b6b',
    lineHeight: 1.3,
  },

  evalButton: {
    width: '100%',
    padding: '6px 0',
    background: '#2d7d46',
    color: '#fff',
    border: 'none',
    borderRadius: '4px',
    fontSize: '13px',
    fontWeight: 'bold',
    cursor: 'pointer',
    letterSpacing: '0.3px',
    transition: 'background 0.15s',
    flexShrink: 0,
  } as React.CSSProperties,

  evalButtonDisabled: {
    background: '#1e3d2a',
    color: '#5a8a6a',
    cursor: 'not-allowed',
  } as React.CSSProperties,
};
