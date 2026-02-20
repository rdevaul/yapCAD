/**
 * RenderModeSelector - Button group for selecting render modes in the toolbar
 */

import { useState } from 'react';

const styles = {
  container: {
    display: 'inline-flex',
    backgroundColor: '#333',
    borderRadius: '4px',
    overflow: 'hidden',
  },
  button: {
    padding: '6px 12px',
    fontSize: '13px',
    backgroundColor: 'transparent',
    color: '#eee',
    border: 'none',
    cursor: 'pointer',
    transition: 'background-color 0.2s',
    whiteSpace: 'nowrap' as const,
  },
  buttonSelected: {
    backgroundColor: '#3b82f6',
    color: 'white',
  },
  buttonHover: {
    backgroundColor: '#444',
  },
  separator: {
    width: '1px',
    backgroundColor: '#555',
  },
};

interface RenderMode {
  id: string;
  name: string;
  description: string;
}

const RENDER_MODES: RenderMode[] = [
  { id: 'solid', name: 'Solid', description: 'Standard solid rendering' },
  { id: 'wireframe', name: 'Wire', description: 'Wireframe only' },
  { id: 'solid+wireframe', name: 'Both', description: 'Solid with wireframe overlay' },
  { id: 'normals', name: 'Normals', description: 'Color-mapped surface normals' },
];

interface RenderModeSelectorProps {
  currentMode: string;
  onModeChange: (mode: string) => void;
}

export function RenderModeSelector({ currentMode, onModeChange }: RenderModeSelectorProps) {
  const [hoveredMode, setHoveredMode] = useState<string | null>(null);

  const handleSelect = (mode: string) => {
    onModeChange(mode);
  };

  return (
    <div style={styles.container}>
      {RENDER_MODES.map((mode, index) => {
        const isSelected = mode.id === currentMode;
        const isHovered = hoveredMode === mode.id;
        const showSeparator = index > 0;

        return (
          <div key={mode.id} style={{ display: 'flex' }}>
            {showSeparator && <div style={styles.separator} />}
            <button
              style={{
                ...styles.button,
                ...(isSelected ? styles.buttonSelected : {}),
                ...(isHovered && !isSelected ? styles.buttonHover : {}),
              }}
              onMouseEnter={() => setHoveredMode(mode.id)}
              onMouseLeave={() => setHoveredMode(null)}
              onClick={() => handleSelect(mode.id)}
              title={mode.description}
            >
              {mode.name}
            </button>
          </div>
        );
      })}
    </div>
  );
}

export default RenderModeSelector;