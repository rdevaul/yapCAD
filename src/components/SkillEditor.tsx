import React, { useState, useEffect } from 'react';

interface SkillParam { name: string; default: string; unit: string; }
interface Skill {
  name: string; description: string; category: string;
  promptTemplate: string; params: SkillParam[];
}
interface SkillEditorProps {
  skill: Skill;
  onSubmit: (prompt: string) => void;
  onClose: () => void;
}

function renderPrompt(template: string, paramValues: Record<string, string>): string {
  let result = template;
  for (const [key, value] of Object.entries(paramValues)) {
    result = result.replace(new RegExp(`\\{\\{${key}\\}\\}`, 'g'), value);
  }
  return result;
}

export function SkillEditor({ skill, onSubmit, onClose }: SkillEditorProps) {
  const [paramValues, setParamValues] = useState<Record<string, string>>(
    () => Object.fromEntries(skill.params.map(p => [p.name, p.default]))
  );
  const [preview, setPreview] = useState(() => renderPrompt(skill.promptTemplate, paramValues));
  const [previewEdited, setPreviewEdited] = useState(false);

  useEffect(() => {
    if (!previewEdited) {
      setPreview(renderPrompt(skill.promptTemplate, paramValues));
    }
  }, [paramValues, previewEdited]);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        onClose();
      }
    };

    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [onClose]);

  const handleParamChange = (name: string, value: string) => {
    setParamValues(prev => ({ ...prev, [name]: value }));
  };

  const handlePreviewChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    setPreview(e.target.value);
    setPreviewEdited(true);
  };

  const handleSubmit = () => {
    onSubmit(preview);
    onClose();
  };

  const handleSave = () => {
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
    const filename = `yapcad-skill-${skill.name}-${timestamp}.json`;
    
    const data = {
      name: `${skill.name}-modified`,
      description: skill.description,
      category: skill.category,
      promptTemplate: preview,
      params: skill.params.map(p => ({
        ...p,
        default: paramValues[p.name] || p.default
      }))
    };

    const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div
      style={{
        position: 'fixed',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        backgroundColor: 'rgba(0,0,0,0.6)',
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        zIndex: 200
      }}
      onClick={onClose}
    >
      <div
        style={{
          backgroundColor: '#1e1e32',
          border: '1px solid #3b3b5a',
          borderRadius: '8px',
          width: '520px',
          maxHeight: '80vh',
          display: 'flex',
          flexDirection: 'column',
          overflow: 'hidden'
        }}
        onClick={(e) => e.stopPropagation()}
      >
        <div style={{ backgroundColor: '#13131f', padding: '12px 16px', display: 'flex', alignItems: 'center' }}>
          <div style={{ color: '#4ade80', fontSize: '16px', fontWeight: 600, flex: 1 }}>
            ⚡ {skill.name} · {skill.category}
          </div>
          <div style={{ color: '#666', cursor: 'pointer' }} onClick={onClose}>
            ✕
          </div>
        </div>

        <div style={{ padding: '16px', overflowY: 'auto', flex: 1 }}>
          <div style={{ color: '#aaa', fontSize: '14px', marginBottom: '16px' }}>
            {skill.description}
          </div>

          <div style={{ marginBottom: '16px' }}>
            <div style={{ color: '#555', fontSize: '10px', textTransform: 'uppercase', marginBottom: '8px' }}>
              PARAMETERS
            </div>
            {skill.params.map((param) => (
              <div key={param.name} style={{ marginBottom: '12px' }}>
                <div style={{ color: '#aaa', fontSize: '12px', marginBottom: '4px' }}>
                  {param.name}
                  {param.unit && ` [${param.unit}]`}
                </div>
                <input
                  type="text"
                  value={paramValues[param.name]}
                  onChange={(e) => handleParamChange(param.name, e.target.value)}
                  style={{
                    backgroundColor: '#23233a',
                    border: '1px solid #3b3b5a',
                    color: '#ddd',
                    borderRadius: '4px',
                    padding: '6px 8px',
                    width: '100%',
                    fontSize: '14px'
                  }}
                />
              </div>
            ))}
          </div>

          <div style={{ marginBottom: '16px' }}>
            <div style={{ color: '#555', fontSize: '10px', textTransform: 'uppercase', marginBottom: '8px' }}>
              PROMPT PREVIEW
            </div>
            <textarea
              value={preview}
              onChange={handlePreviewChange}
              style={{
                backgroundColor: '#0e0e1c',
                border: '1px solid #2a2a48',
                color: '#a8d8a8',
                fontFamily: 'monospace',
                fontSize: '14px',
                padding: '12px',
                width: '100%',
                minHeight: '140px',
                resize: 'vertical'
              }}
            />
          </div>
        </div>

        <div style={{ padding: '12px 16px', display: 'flex', gap: '8px', borderTop: '1px solid #3b3b5a' }}>
          <button
            onClick={handleSave}
            style={{
              backgroundColor: '#2a2a42',
              color: '#aaa',
              border: 'none',
              borderRadius: '4px',
              padding: '8px 12px',
              fontSize: '14px',
              cursor: 'pointer'
            }}
          >
            💾 Save copy
          </button>
          <button
            onClick={onClose}
            style={{
              backgroundColor: 'transparent',
              color: '#666',
              border: '1px solid #333',
              borderRadius: '4px',
              padding: '8px 12px',
              fontSize: '14px',
              cursor: 'pointer'
            }}
          >
            ✕ Cancel
          </button>
          <button
            onClick={handleSubmit}
            style={{
              backgroundColor: '#4ade80',
              color: '#000',
              border: 'none',
              borderRadius: '4px',
              padding: '8px 12px',
              fontSize: '14px',
              fontWeight: 600,
              cursor: 'pointer',
              marginLeft: 'auto'
            }}
          >
            ▶ Submit
          </button>
        </div>
      </div>
    </div>
  );
}

// renderPrompt defined above