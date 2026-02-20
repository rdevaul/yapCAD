/**
 * DSL Editor Component - Monaco Editor for yapCAD DSL
 */

import { useState, useCallback, useRef, useEffect } from 'react';
import Editor from '@monaco-editor/react';
import * as monaco from 'monaco-editor';

interface DSLEditorProps {
  onEvaluate: (source: string) => Promise<void>;
  onSourceChange?: (source: string) => void;
  externalSource?: string | null; // When set, loads this source into the editor
  onExternalSourceConsumed?: () => void; // Called after consuming externalSource
  isEvaluating: boolean;
  evaluationError: string | null;
}

const styles = {
  container: {
    display: 'flex',
    flexDirection: 'column' as const,
    height: '100%',
    backgroundColor: '#1a1a2e',
    border: '1px solid #333',
  },
  toolbar: {
    display: 'flex',
    alignItems: 'center',
    padding: '8px 12px',
    backgroundColor: '#252540',
    borderBottom: '1px solid #333',
    gap: '12px',
  },
  evaluateButton: {
    padding: '6px 12px',
    fontSize: '13px',
    backgroundColor: '#4c9f38',
    color: '#fff',
    border: 'none',
    borderRadius: '4px',
    cursor: 'pointer',
    display: 'flex',
    alignItems: 'center',
    gap: '6px',
    transition: 'background-color 0.2s',
  },
  evaluateButtonDisabled: {
    backgroundColor: '#666',
    cursor: 'not-allowed',
  },
  statusText: {
    fontSize: '12px',
    color: '#888',
  },
  editorContainer: {
    flex: 1,
    position: 'relative' as const,
    minHeight: 0,
  },
  errorPanel: {
    backgroundColor: '#3d1a1a',
    border: '1px solid #cc3333',
    borderRadius: '4px',
    padding: '8px 12px',
    margin: '8px',
    color: '#ff6666',
    fontSize: '12px',
    fontFamily: 'monospace',
    whiteSpace: 'pre-wrap' as const,
  },
};

// yapCAD DSL language configuration
const DSL_LANGUAGE_CONFIG = {
  keywords: [
    'command', 'module', 'let', 'emit', 'require', 'for', 'if', 'else', 'elif', 'return',
    'true', 'false', 'nil', 'function', 'var', 'const'
  ],
  builtins: [
    // Basic shapes
    'box', 'sphere', 'cylinder', 'cone', 'torus', 'tube', 'conic_tube',
    'spherical_shell', 'dodecahedron',
    // Boolean operations
    'union', 'difference', 'intersection',
    // Transformations
    'translate', 'rotate', 'scale', 'mirror', 'scale_uniform',
    'rotate_2d', 'mirror_2d', 'mirror_y',
    // Advanced operations
    'extrude', 'revolve', 'loft', 'helical_extrude', 'bezier',
    // Path operations
    'path3d_eval', 'path3d_length', 'split_solid',
    // Text
    'text_solid', 'engrave_text', 'text_width',
    // Utilities
    'point', 'point2d', 'vector', 'vector2d', 'centroid', 'distance',
  ],
};

export function DSLEditor({ onEvaluate, onSourceChange, externalSource, onExternalSourceConsumed, isEvaluating, evaluationError }: DSLEditorProps) {
  const [source, setSource] = useState<string>('');
  const [hasUnsavedChanges, setHasUnsavedChanges] = useState(false);
  const editorRef = useRef<monaco.editor.IStandaloneCodeEditor | null>(null);
  const debounceTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Load saved source from localStorage on mount
  useEffect(() => {
    const saved = localStorage.getItem('yapcad-dsl-source');
    if (saved) {
      setSource(saved);
      onSourceChange?.(saved);
    } else {
      // Default content for new users
      const defaultContent = `# Welcome to yapCAD DSL Editor!
# Uses Python-style indentation (colon + indent, no braces)

command myBox() -> solid:
    let b = box(20, 20, 20)
    emit b

# Press Ctrl+Enter to evaluate, or click the "Evaluate" button`;
      setSource(defaultContent);
      onSourceChange?.(defaultContent);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Handle external source loading (from FileBar)
  useEffect(() => {
    if (externalSource != null) {
      setSource(externalSource);
      onSourceChange?.(externalSource);
      localStorage.setItem('yapcad-dsl-source', externalSource);
      setHasUnsavedChanges(false);
      onExternalSourceConsumed?.();
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [externalSource]);

  // Debounced auto-save
  const debouncedSave = useCallback((value: string) => {
    if (debounceTimeoutRef.current) {
      clearTimeout(debounceTimeoutRef.current);
    }
    
    debounceTimeoutRef.current = setTimeout(() => {
      localStorage.setItem('yapcad-dsl-source', value);
      setHasUnsavedChanges(false);
    }, 1000);
  }, []);

  const handleEditorChange = useCallback((value: string | undefined) => {
    const newValue = value || '';
    setSource(newValue);
    setHasUnsavedChanges(true);
    debouncedSave(newValue);
    onSourceChange?.(newValue);
  }, [debouncedSave, onSourceChange]);

  const handleEvaluate = useCallback(async () => {
    if (!source.trim() || isEvaluating) return;
    
    try {
      await onEvaluate(source);
    } catch (error) {
      console.error('Evaluation failed:', error);
    }
  }, [source, isEvaluating, onEvaluate]);

  const handleEditorDidMount = useCallback((editor: monaco.editor.IStandaloneCodeEditor, monaco: typeof import('monaco-editor')) => {
    editorRef.current = editor;
    
    // Configure yapCAD DSL language
    monaco.languages.register({ id: 'yapcad-dsl' });
    
    monaco.languages.setMonarchTokensProvider('yapcad-dsl', {
      tokenizer: {
        root: [
          // Keywords
          [/\b(?:command|module|let|emit|require|for|if|else|elif|return|true|false|nil)\b/, 'keyword'],
          [/\b(?:solid|float|int|bool|string|point|vector|path|surface|shell|transform)\b/, 'type.identifier'],
          
          // Built-in functions
          [/\b(?:box|sphere|cylinder|cone|torus|tube|conic_tube|spherical_shell|dodecahedron|union|difference|intersection|translate|rotate|scale|mirror|extrude|revolve|loft|helical_extrude|bezier|point|vector|centroid|distance|text_solid|engrave_text)\b/, 'type'],
          
          // Numbers
          [/\d*\.\d+([eE][\-+]?\d+)?/, 'number.float'],
          [/\d+/, 'number'],
          
          // Strings
          [/"([^"\\]|\\.)*$/, 'string.invalid'],
          [/"/, { token: 'string.quote', bracket: '@open', next: '@string' }],
          
          // Comments
          [/\/\/.*$/, 'comment'],
          [/\/\*/, 'comment', '@comment'],
          
          // Identifiers
          [/[a-zA-Z_]\w*/, 'identifier'],
          
          // Operators
          [/[{}()\[\]]/, '@brackets'],
          [/[;,.]/, 'delimiter'],
          [/[=!<>]=?/, 'operator'],
          [/[+\-*/]/, 'operator'],
        ],
        
        string: [
          [/[^\\"]+/, 'string'],
          [/\\./, 'string.escape.invalid'],
          [/"/, { token: 'string.quote', bracket: '@close', next: '@pop' }]
        ],
        
        comment: [
          [/[^\/*]+/, 'comment'],
          [/\/\*/, 'comment', '@push'],
          [/\*\//, 'comment', '@pop'],
          [/[\/*]/, 'comment']
        ],
      },
    });
    
    // Configure auto-completion
    monaco.languages.registerCompletionItemProvider('yapcad-dsl', {
      provideCompletionItems: (_, position) => {
        const suggestions: monaco.languages.CompletionItem[] = [];
        const range = {
          startLineNumber: position.lineNumber,
          endLineNumber: position.lineNumber,
          startColumn: position.column,
          endColumn: position.column,
        };
        
        // Add keywords
        DSL_LANGUAGE_CONFIG.keywords.forEach(keyword => {
          suggestions.push({
            label: keyword,
            kind: monaco.languages.CompletionItemKind.Keyword,
            insertText: keyword,
            detail: 'Keyword',
            range: range,
          });
        });
        
        // Add built-in functions
        DSL_LANGUAGE_CONFIG.builtins.forEach(builtin => {
          suggestions.push({
            label: builtin,
            kind: monaco.languages.CompletionItemKind.Function,
            insertText: `${builtin}()`,
            detail: 'Built-in function',
            insertTextRules: monaco.languages.CompletionItemInsertTextRule.InsertAsSnippet,
            range: range,
          });
        });
        
        return { suggestions };
      },
    });
    
    // Set theme
    monaco.editor.defineTheme('yapcad-dark', {
      base: 'vs-dark',
      inherit: true,
      rules: [
        { token: 'keyword', foreground: 'c678dd' },
        { token: 'type', foreground: '61dafb' },
        { token: 'string', foreground: '98c379' },
        { token: 'comment', foreground: '7c7c7c' },
        { token: 'number', foreground: 'd19a66' },
        { token: 'operator', foreground: '56b6c2' },
      ],
      colors: {
        'editor.background': '#1a1a2e',
        'editor.foreground': '#eeeeee',
        'editorLineNumber.foreground': '#666666',
        'editorCursor.foreground': '#eeeeee',
        'editor.selectionBackground': '#3d4454',
        'editor.lineHighlightBackground': '#252540',
      },
    });
    
    monaco.editor.setTheme('yapcad-dark');
    
    // Configure editor options
    editor.updateOptions({
      fontSize: 13,
      lineHeight: 20,
      tabSize: 2,
      insertSpaces: true,
      wordWrap: 'on',
      minimap: { enabled: false },
      folding: true,
      lineNumbers: 'on',
      renderWhitespace: 'selection',
      scrollBeyondLastLine: false,
    });
    
    // Keyboard shortcuts
    editor.addCommand(monaco.KeyMod.CtrlCmd | monaco.KeyCode.Enter, () => {
      handleEvaluate();
    });
  }, [handleEvaluate]);

  return (
    <div style={styles.container}>
      <div style={styles.toolbar}>
        <span style={styles.statusText}>
          {isEvaluating ? '⏳ Evaluating...' : hasUnsavedChanges ? '● Unsaved' : '✓ Saved'}
        </span>
        {source.trim() && (
          <span style={styles.statusText}>
            {source.split('\n').length} lines
          </span>
        )}
        <span style={{ ...styles.statusText, marginLeft: 'auto', fontSize: '10px', color: '#555' }}>
          Ctrl+Enter to evaluate
        </span>
      </div>
      
      <div style={styles.editorContainer}>
        <Editor
          height="100%"
          language="yapcad-dsl"
          value={source}
          onChange={handleEditorChange}
          onMount={handleEditorDidMount}
          loading={<div style={{ padding: '16px', color: '#888' }}>Loading editor...</div>}
        />
      </div>
      
      {evaluationError && (
        <div style={styles.errorPanel}>
          <strong>Evaluation Error:</strong>
          {evaluationError}
        </div>
      )}
    </div>
  );
}