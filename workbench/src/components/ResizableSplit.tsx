/**
 * Resizable Split Panel Component
 * Allows resizing between two panels with a draggable handle
 */

import { useState, useCallback, useRef, useEffect } from 'react';

interface ResizableSplitProps {
  direction: 'horizontal' | 'vertical';
  initialSplit?: number; // 0-1, default position
  minSize?: number; // Minimum size in pixels for first panel
  maxSize?: number; // Maximum size in pixels for first panel  
  children: [React.ReactNode, React.ReactNode];
  onSplitChange?: (split: number) => void;
  className?: string;
  style?: React.CSSProperties;
}

const styles = {
  container: {
    display: 'flex',
    width: '100%',
    height: '100%',
    overflow: 'hidden',
  },
  containerHorizontal: {
    flexDirection: 'row' as const,
  },
  containerVertical: {
    flexDirection: 'column' as const,
  },
  panel: {
    overflow: 'hidden',
    display: 'flex',
    flexDirection: 'column' as const,
    minHeight: 0,
    minWidth: 0,
  },
  handle: {
    backgroundColor: '#444',
    flexShrink: 0,
    position: 'relative' as const,
    transition: 'background-color 0.2s',
  },
  handleHorizontal: {
    width: '4px',
    cursor: 'col-resize',
    borderLeft: '1px solid #333',
    borderRight: '1px solid #333',
  },
  handleVertical: {
    height: '4px',
    cursor: 'row-resize',
    borderTop: '1px solid #333',
    borderBottom: '1px solid #333',
  },
  handleActive: {
    backgroundColor: '#4c9f38',
  },
  handleHover: {
    backgroundColor: '#555',
  },
};

export function ResizableSplit({
  direction,
  initialSplit = 0.6,
  minSize = 100,
  maxSize,
  children,
  onSplitChange,
  className,
  style,
}: ResizableSplitProps) {
  const [split, setSplit] = useState(initialSplit);
  const [isDragging, setIsDragging] = useState(false);
  const [isHovered, setIsHovered] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);
  const startPositionRef = useRef<number>(0);
  const startSplitRef = useRef<number>(0);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    setIsDragging(true);
    
    const startPos = direction === 'horizontal' ? e.clientX : e.clientY;
    startPositionRef.current = startPos;
    startSplitRef.current = split;
    
    document.body.style.cursor = direction === 'horizontal' ? 'col-resize' : 'row-resize';
    document.body.style.userSelect = 'none';
  }, [direction, split]);

  const handleMouseMove = useCallback((e: MouseEvent) => {
    if (!isDragging || !containerRef.current) return;

    const container = containerRef.current;
    const containerSize = direction === 'horizontal' 
      ? container.clientWidth 
      : container.clientHeight;
    
    const currentPos = direction === 'horizontal' ? e.clientX : e.clientY;
    const delta = currentPos - startPositionRef.current;
    const deltaRatio = delta / containerSize;
    
    let newSplit = startSplitRef.current + deltaRatio;
    
    // Apply constraints
    const minRatio = minSize / containerSize;
    const maxRatio = maxSize ? maxSize / containerSize : 0.9;
    
    newSplit = Math.max(minRatio, Math.min(maxRatio, newSplit));
    
    setSplit(newSplit);
    onSplitChange?.(newSplit);
  }, [isDragging, direction, minSize, maxSize, onSplitChange]);

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
    document.body.style.cursor = '';
    document.body.style.userSelect = '';
  }, []);

  // Mouse event listeners
  useEffect(() => {
    if (isDragging) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
      
      return () => {
        document.removeEventListener('mousemove', handleMouseMove);
        document.removeEventListener('mouseup', handleMouseUp);
      };
    }
  }, [isDragging, handleMouseMove, handleMouseUp]);

  const containerStyle = {
    ...styles.container,
    ...(direction === 'horizontal' ? styles.containerHorizontal : styles.containerVertical),
    ...style,
  };

  const firstPanelStyle = {
    ...styles.panel,
    ...(direction === 'horizontal' 
      ? { width: `${split * 100}%` }
      : { height: `${split * 100}%` }
    ),
  };

  const secondPanelStyle = {
    ...styles.panel,
    ...(direction === 'horizontal' 
      ? { width: `${(1 - split) * 100}%` }
      : { height: `${(1 - split) * 100}%` }
    ),
  };

  const handleStyle = {
    ...styles.handle,
    ...(direction === 'horizontal' ? styles.handleHorizontal : styles.handleVertical),
    ...(isDragging ? styles.handleActive : {}),
    ...(isHovered && !isDragging ? styles.handleHover : {}),
  };

  return (
    <div ref={containerRef} className={className} style={containerStyle}>
      <div style={firstPanelStyle}>
        {children[0]}
      </div>
      
      <div
        style={handleStyle}
        onMouseDown={handleMouseDown}
        onMouseEnter={() => setIsHovered(true)}
        onMouseLeave={() => setIsHovered(false)}
      />
      
      <div style={secondPanelStyle}>
        {children[1]}
      </div>
    </div>
  );
}