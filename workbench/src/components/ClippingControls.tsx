/**
 * ClippingControls - Collapsible panel for managing clipping planes
 */

import { useState, useEffect, useRef } from 'react';

const styles = {
  container: {
    backgroundColor: '#1a1a2e',
    border: '1px solid #333',
    borderRadius: '6px',
    overflow: 'hidden',
  },
  header: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
    padding: '8px 12px',
    backgroundColor: '#252540',
    cursor: 'pointer',
    fontSize: '13px',
    fontWeight: 500,
    color: '#eee',
  },
  headerIcon: {
    fontSize: '10px',
    color: '#888',
    transition: 'transform 0.2s',
  },
  headerIconExpanded: {
    transform: 'rotate(180deg)',
  },
  content: {
    padding: '12px',
    borderTop: '1px solid #333',
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '12px',
  },
  axisControl: {
    display: 'flex',
    flexDirection: 'column' as const,
    gap: '8px',
  },
  axisHeader: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
  },
  axisLabel: {
    fontSize: '13px',
    fontWeight: 500,
    color: '#eee',
    minWidth: '16px',
  },
  checkbox: {
    width: '14px',
    height: '14px',
    accentColor: '#3b82f6',
  },
  sliderContainer: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    marginLeft: '24px',
  },
  slider: {
    flex: 1,
    height: '4px',
    borderRadius: '2px',
    background: '#333',
    appearance: 'none' as const,
    outline: 'none',
  },
  sliderValue: {
    fontSize: '11px',
    color: '#888',
    fontFamily: 'monospace',
    minWidth: '32px',
    textAlign: 'right' as const,
  },
  invertContainer: {
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    marginLeft: '24px',
  },
  invertLabel: {
    fontSize: '11px',
    color: '#aaa',
  },
  separator: {
    height: '1px',
    backgroundColor: '#333',
    margin: '4px 0',
  },
  resetButton: {
    padding: '6px 12px',
    fontSize: '12px',
    backgroundColor: '#333',
    color: '#eee',
    border: 'none',
    borderRadius: '4px',
    cursor: 'pointer',
    transition: 'background-color 0.2s',
    alignSelf: 'flex-start',
  },
  resetButtonHover: {
    backgroundColor: '#444',
  },
};

interface AxisState {
  enabled: boolean;
  position: number;
  invert: boolean;
}

interface ClippingState {
  x: AxisState;
  y: AxisState;
  z: AxisState;
}

interface ClippingControlsProps {
  onClipPlaneChange: (axis: 'x' | 'y' | 'z', enabled: boolean, position: number, invert: boolean) => void;
  onResetClipPlanes: () => void;
  getClippingState: () => { axis: string; enabled: boolean; position: number; invert: boolean }[];
}

export function ClippingControls({
  onClipPlaneChange,
  onResetClipPlanes,
  getClippingState,
}: ClippingControlsProps) {
  const [isExpanded, setIsExpanded] = useState(false);
  const [hoveredReset, setHoveredReset] = useState(false);
  const [state, setState] = useState<ClippingState>({
    x: { enabled: false, position: 0.5, invert: false },
    y: { enabled: false, position: 0.5, invert: false },
    z: { enabled: false, position: 0.5, invert: false },
  });

  // Sync from parent on mount only (avoid re-render loops)
  const didMount = useRef(false);
  useEffect(() => {
    if (didMount.current) return;
    didMount.current = true;
    const clippingState = getClippingState();
    if (clippingState.length === 0) return;
    const newState: ClippingState = { x: { enabled: false, position: 0.5, invert: false }, y: { enabled: false, position: 0.5, invert: false }, z: { enabled: false, position: 0.5, invert: false } };
    for (const axisState of clippingState) {
      const axis = axisState.axis as 'x' | 'y' | 'z';
      newState[axis] = {
        enabled: axisState.enabled,
        position: axisState.position,
        invert: axisState.invert,
      };
    }
    setState(newState);
  }, [getClippingState]);

  const handleToggle = () => {
    setIsExpanded(!isExpanded);
  };

  const handleAxisToggle = (axis: 'x' | 'y' | 'z') => {
    const newEnabled = !state[axis].enabled;
    const newState = {
      ...state,
      [axis]: { ...state[axis], enabled: newEnabled },
    };
    setState(newState);
    onClipPlaneChange(axis, newEnabled, state[axis].position, state[axis].invert);
  };

  const handlePositionChange = (axis: 'x' | 'y' | 'z', position: number) => {
    const newState = {
      ...state,
      [axis]: { ...state[axis], position },
    };
    setState(newState);
    onClipPlaneChange(axis, state[axis].enabled, position, state[axis].invert);
  };

  const handleInvertToggle = (axis: 'x' | 'y' | 'z') => {
    const newInvert = !state[axis].invert;
    const newState = {
      ...state,
      [axis]: { ...state[axis], invert: newInvert },
    };
    setState(newState);
    onClipPlaneChange(axis, state[axis].enabled, state[axis].position, newInvert);
  };

  const handleReset = () => {
    const newState: ClippingState = {
      x: { enabled: false, position: 0.5, invert: false },
      y: { enabled: false, position: 0.5, invert: false },
      z: { enabled: false, position: 0.5, invert: false },
    };
    setState(newState);
    onResetClipPlanes();
  };

  const axes: Array<{ key: 'x' | 'y' | 'z'; label: string; color: string }> = [
    { key: 'x', label: 'X', color: '#ff4444' },
    { key: 'y', label: 'Y', color: '#44ff44' },
    { key: 'z', label: 'Z', color: '#4444ff' },
  ];

  const hasAnyEnabled = axes.some(axis => state[axis.key].enabled);

  return (
    <div style={styles.container}>
      <div style={styles.header} onClick={handleToggle}>
        <span>Clipping Planes</span>
        <span style={{
          ...styles.headerIcon,
          ...(isExpanded ? styles.headerIconExpanded : {}),
        }}>
          ▼
        </span>
      </div>

      {isExpanded && (
        <div style={styles.content}>
          {axes.map((axis, index) => {
            const axisState = state[axis.key];
            
            return (
              <div key={axis.key} style={styles.axisControl}>
                <div style={styles.axisHeader}>
                  <span style={{ ...styles.axisLabel, color: axis.color }}>
                    {axis.label}
                  </span>
                  <input
                    type="checkbox"
                    checked={axisState.enabled}
                    onChange={() => handleAxisToggle(axis.key)}
                    style={styles.checkbox}
                  />
                </div>

                {axisState.enabled && (
                  <>
                    <div style={styles.sliderContainer}>
                      <input
                        type="range"
                        min="0"
                        max="1"
                        step="0.01"
                        value={axisState.position}
                        onChange={(e) => handlePositionChange(axis.key, parseFloat(e.target.value))}
                        style={{
                          ...styles.slider,
                          background: `linear-gradient(to right, ${axis.color}40 0%, ${axis.color}80 ${axisState.position * 100}%, #333 ${axisState.position * 100}%, #333 100%)`,
                        }}
                      />
                      <span style={styles.sliderValue}>
                        {(axisState.position * 100).toFixed(0)}%
                      </span>
                    </div>

                    <div style={styles.invertContainer}>
                      <input
                        type="checkbox"
                        checked={axisState.invert}
                        onChange={() => handleInvertToggle(axis.key)}
                        style={styles.checkbox}
                      />
                      <span style={styles.invertLabel}>Invert</span>
                    </div>
                  </>
                )}

                {index < axes.length - 1 && <div style={styles.separator} />}
              </div>
            );
          })}

          {hasAnyEnabled && (
            <>
              <div style={styles.separator} />
              <button
                style={{
                  ...styles.resetButton,
                  ...(hoveredReset ? styles.resetButtonHover : {}),
                }}
                onMouseEnter={() => setHoveredReset(true)}
                onMouseLeave={() => setHoveredReset(false)}
                onClick={handleReset}
              >
                Reset All
              </button>
            </>
          )}
        </div>
      )}
    </div>
  );
}

export default ClippingControls;