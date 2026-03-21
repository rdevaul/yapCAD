import React, { useState, useEffect, useMemo } from 'react';
import skills from '../data/yapcad-skills.json';

interface Skill {
  name: string;
  description: string;
  category: string;
  promptTemplate: string;
  params: { name: string; default: string; unit: string }[];
}

interface CommandPaletteProps {
  filter: string;
  onSelectCommand: (cmd: string) => void;
  onSelectSkill: (skill: Skill) => void;
  onDismiss: () => void;
}

const BUILTIN_COMMANDS = [
  { name: 'new',     description: 'Start a fresh agent session',   icon: '🔄' },
  { name: 'status',  description: 'Show session status & token usage', icon: '📊' },
  { name: 'compact', description: 'Compact context window',         icon: '🗜' },
  { name: 'help',    description: 'Show available commands',        icon: '❓' },
];

const CommandPalette: React.FC<CommandPaletteProps> = ({ 
  filter, 
  onSelectCommand, 
  onSelectSkill, 
  onDismiss 
}) => {
  const [selectedIndex, setSelectedIndex] = useState(0);
  
  const filteredItems = useMemo(() => {
    const lowerFilter = filter.toLowerCase();
    
    const commands = BUILTIN_COMMANDS.filter(cmd => 
      cmd.name.toLowerCase().includes(lowerFilter) || 
      cmd.description.toLowerCase().includes(lowerFilter)
    ).map(cmd => ({
      type: 'command',
      ...cmd
    }));
    
    const skillItems = skills.filter(skill => 
      skill.name.toLowerCase().includes(lowerFilter) || 
      skill.description.toLowerCase().includes(lowerFilter)
    ).map(skill => ({
      type: 'skill',
      ...skill
    }));
    
    return [...commands, ...skillItems];
  }, [filter]);
  
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if (e.key === 'Escape') { 
        onDismiss(); 
        return; 
      }
      if (e.key === 'ArrowDown') { 
        e.preventDefault(); 
        setSelectedIndex(i => Math.min(i + 1, filteredItems.length - 1)); 
      }
      if (e.key === 'ArrowUp') { 
        e.preventDefault(); 
        setSelectedIndex(i => Math.max(i - 1, 0)); 
      }
      if (e.key === 'Enter') { 
        e.preventDefault(); 
        if (filteredItems.length > 0) {
          const selectedItem = filteredItems[selectedIndex];
          if (selectedItem.type === 'command') {
            onSelectCommand(`/${selectedItem.name}`);
          } else {
            onSelectSkill(selectedItem as Skill);
          }
        }
      }
    };
    
    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, [selectedIndex, filteredItems, onDismiss, onSelectCommand, onSelectSkill]);
  
  const getGroupedItems = () => {
    const groups: any[] = [];
    
    if (filteredItems.some(item => item.type === 'command')) {
      groups.push({
        type: 'header',
        label: 'COMMANDS'
      });
      
      groups.push(...filteredItems.filter(item => item.type === 'command'));
    }
    
    const skillGroups = filteredItems
      .filter(item => item.type === 'skill')
      .reduce((acc, item) => {
        const skill = item as Skill & { type: string };
        if (!acc[skill.category]) {
          acc[skill.category] = [];
        }
        acc[skill.category].push(skill);
        return acc;
      }, {} as Record<string, (Skill & { type: string })[]>);
    
    Object.entries(skillGroups).forEach(([category, categorySkills]) => {
      groups.push({
        type: 'header',
        label: `SKILLS — ${category}`
      });
      groups.push(...categorySkills);
    });
    
    return groups;
  };
  
  const groupedItems = getGroupedItems();
  
  return (
    <div 
      style={{
        position: 'absolute',
        bottom: '100%',
        left: 0,
        right: 0,
        zIndex: 100,
        backgroundColor: '#1e1e32',
        border: '1px solid #3b3b5a',
        borderRadius: '6px',
        boxShadow: '0 -4px 16px rgba(0,0,0,0.4)',
        maxHeight: '240px',
        overflowY: 'auto',
        padding: '4px 0',
        fontFamily: 'monospace',
        fontSize: '12px'
      }}
    >
      {groupedItems.map((item, index) => {
        if (item.type === 'header') {
          return (
            <div 
              key={index}
              style={{
                padding: '4px 8px',
                color: '#555',
                fontSize: '10px',
                textTransform: 'uppercase',
                letterSpacing: '0.08em',
                fontWeight: 600
              }}
            >
              {item.label}
            </div>
          );
        }
        
        const isSelected = index === selectedIndex;
        
        return (
          <div
            key={index}
            style={{
              display: 'flex',
              alignItems: 'center',
              padding: '6px 8px',
              backgroundColor: isSelected ? '#2a2a48' : 'transparent',
              cursor: 'pointer'
            }}
            onClick={() => {
              if (item.type === 'command') {
                onSelectCommand(`/${item.name}`);
              } else {
                onSelectSkill(item);
              }
            }}
          >
            <span style={{ 
              color: item.type === 'command' ? '#7c7cf8' : '#4ade80',
              fontWeight: 600,
              width: '100px',
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap'
            }}>
              {item.type === 'command' ? `/${item.name}` : item.name}
            </span>
            <span style={{ 
              color: '#888',
              fontSize: '11px',
              marginLeft: '8px',
              flex: 1,
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap'
            }}>
              {item.description}
            </span>
          </div>
        );
      })}
    </div>
  );
};

export default CommandPalette;