/**
 * LightingSelector - Dropdown for selecting lighting presets in the toolbar
 */

import { useState, useRef, useEffect } from 'react';
import { getPresetNames, LIGHTING_PRESETS } from '../viewer/lighting';

const styles = {
  container: {
    position: 'relative' as const,
    display: 'inline-block',
  },
  button: {
    padding: '6px 12px',
    fontSize: '13px',
    backgroundColor: '#333',
    color: '#eee',
    border: 'none',
    borderRadius: '4px',
    cursor: 'pointer',
    transition: 'background-color 0.2s',
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    minWidth: '100px',
    justifyContent: 'space-between',
  },
  buttonHover: {
    backgroundColor: '#444',
  },
  dropdown: {
    position: 'absolute' as const,
    top: '100%',
    left: 0,
    right: 0,
    marginTop: '4px',
    backgroundColor: '#2a2a4a',
    border: '1px solid #444',
    borderRadius: '6px',
    boxShadow: '0 4px 12px rgba(0, 0, 0, 0.4)',
    zIndex: 1000,
    overflow: 'hidden',
  },
  option: {
    padding: '8px 12px',
    fontSize: '13px',
    color: '#eee',
    cursor: 'pointer',
    transition: 'background-color 0.15s',
    borderBottom: '1px solid #333',
  },
  optionHover: {
    backgroundColor: '#333',
  },
  optionSelected: {
    backgroundColor: '#3b82f6',
    color: 'white',
  },
  optionLast: {
    borderBottom: 'none',
  },
  arrow: {
    fontSize: '10px',
    color: '#888',
    transition: 'transform 0.2s',
  },
  arrowOpen: {
    transform: 'rotate(180deg)',
  },
};

interface LightingSelectorProps {
  currentPreset: string;
  onPresetChange: (preset: string) => void;
}

export function LightingSelector({ currentPreset, onPresetChange }: LightingSelectorProps) {
  const [isOpen, setIsOpen] = useState(false);
  const [hoveredOption, setHoveredOption] = useState<string | null>(null);

  const presetNames = getPresetNames();
  const currentPresetData = LIGHTING_PRESETS[currentPreset];

  const handleToggle = () => {
    setIsOpen(!isOpen);
  };

  const handleSelect = (preset: string) => {
    onPresetChange(preset);
    setIsOpen(false);
    setHoveredOption(null);
  };

  // Close on click outside
  const containerRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!isOpen) return;
    const handleClickOutside = (e: MouseEvent) => {
      if (containerRef.current && !containerRef.current.contains(e.target as Node)) {
        setIsOpen(false);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, [isOpen]);

  return (
    <div style={styles.container} ref={containerRef}>
      <button
        style={styles.button}
        onClick={handleToggle}
        title={currentPresetData?.description || 'Select lighting preset'}
      >
        <span>{currentPresetData?.name || currentPreset}</span>
        <span style={{
          ...styles.arrow,
          ...(isOpen ? styles.arrowOpen : {}),
        }}>
          ▼
        </span>
      </button>

      {isOpen && (
        <div style={{ ...styles.dropdown, minWidth: '160px' }}>
          {presetNames.map((presetName, index) => {
            const preset = LIGHTING_PRESETS[presetName];
            const isSelected = presetName === currentPreset;
            const isHovered = hoveredOption === presetName;
            const isLast = index === presetNames.length - 1;

            return (
              <div
                key={presetName}
                style={{
                  ...styles.option,
                  ...(isLast ? styles.optionLast : {}),
                  ...(isSelected ? styles.optionSelected : {}),
                  ...(isHovered && !isSelected ? styles.optionHover : {}),
                }}
                onMouseEnter={() => setHoveredOption(presetName)}
                onMouseLeave={() => setHoveredOption(null)}
                onClick={() => handleSelect(presetName)}
                title={preset.description}
              >
                {preset.name}
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}

export default LightingSelector;