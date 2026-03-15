/**
 * Multi-file tab state management for yapCAD workbench.
 *
 * Each tab represents an independent DSL editing context with:
 *   - filename (editable)
 *   - DSL source content
 *   - parsed command list
 *   - active command selection
 *   - parameter values per command
 *
 * State is persisted in sessionStorage under key `yapcad-tabs`.
 */

import { useState, useCallback, useRef, useEffect } from 'react';

// ── Types ─────────────────────────────────────────────────────────────────────

export interface CommandParam {
  name: string;
  type: string;
  default: number | string | boolean | null;
  ui_hint?: {
    min?: number;
    max?: number;
    step?: number;
    kind?: string;
  };
}

export interface CommandInfo {
  name: string;
  params: CommandParam[];
}

export interface TabState {
  id: string;
  filename: string;
  source: string;
  commands: CommandInfo[];
  selectedCommand: string;
  paramValues: Record<string, Record<string, unknown>>; // command -> param -> value
  isParsing: boolean;
  parseError: string;
}

interface SerializedTabState {
  id: string;
  filename: string;
  source: string;
  selectedCommand: string;
  paramValues: Record<string, Record<string, unknown>>;
}

const SESSION_KEY = 'yapcad-tabs';
const SESSION_ACTIVE_KEY = 'yapcad-active-tab';

let tabCounter = 0;

function generateTabId(): string {
  return `tab_${Date.now()}_${++tabCounter}`;
}

function createDefaultTab(index: number): TabState {
  const id = generateTabId();
  return {
    id,
    filename: `untitled_${index}`,
    source: `module untitled_${index}\n\n# Write your yapCAD DSL here\ncommand main() -> solid:\n    let b = box(20, 20, 20)\n    emit b\n`,
    commands: [],
    selectedCommand: '',
    paramValues: {},
    isParsing: false,
    parseError: '',
  };
}

// ── Persistence ───────────────────────────────────────────────────────────────

function saveTabs(tabs: TabState[], activeTabId: string): void {
  try {
    const serialized: SerializedTabState[] = tabs.map(t => ({
      id: t.id,
      filename: t.filename,
      source: t.source,
      selectedCommand: t.selectedCommand,
      paramValues: t.paramValues,
    }));
    sessionStorage.setItem(SESSION_KEY, JSON.stringify(serialized));
    sessionStorage.setItem(SESSION_ACTIVE_KEY, activeTabId);
  } catch {
    // sessionStorage full or unavailable
  }
}

function loadTabs(): { tabs: TabState[]; activeTabId: string } | null {
  try {
    const raw = sessionStorage.getItem(SESSION_KEY);
    const activeId = sessionStorage.getItem(SESSION_ACTIVE_KEY);
    if (!raw) return null;

    const serialized: SerializedTabState[] = JSON.parse(raw);
    if (!Array.isArray(serialized) || serialized.length === 0) return null;

    const tabs: TabState[] = serialized.map(s => ({
      id: s.id,
      filename: s.filename,
      source: s.source,
      commands: [],
      selectedCommand: s.selectedCommand,
      paramValues: s.paramValues,
      isParsing: false,
      parseError: '',
    }));

    return {
      tabs,
      activeTabId: activeId && tabs.some(t => t.id === activeId) ? activeId : tabs[0].id,
    };
  } catch {
    return null;
  }
}

// ── Hook ──────────────────────────────────────────────────────────────────────

interface UseTabStateOptions {
  backendUrl: string;
}

export function useTabState({ backendUrl }: UseTabStateOptions) {
  const [tabs, setTabs] = useState<TabState[]>(() => {
    const saved = loadTabs();
    if (saved) return saved.tabs;
    return [createDefaultTab(1)];
  });

  const [activeTabId, setActiveTabId] = useState<string>(() => {
    const saved = loadTabs();
    if (saved) return saved.activeTabId;
    return tabs[0]?.id ?? '';
  });

  // Track in-flight parse requests for cancellation
  const parseAbortRef = useRef<Map<string, AbortController>>(new Map());

  // Persist on change
  useEffect(() => {
    saveTabs(tabs, activeTabId);
  }, [tabs, activeTabId]);

  // Active tab helper
  const activeTab = tabs.find(t => t.id === activeTabId) ?? tabs[0];

  // ── Tab operations ──────────────────────────────────────────────────────────

  const addTab = useCallback(() => {
    const index = tabs.length + 1;
    const newTab = createDefaultTab(index);
    setTabs(prev => [...prev, newTab]);
    setActiveTabId(newTab.id);
    return newTab;
  }, [tabs.length]);

  const closeTab = useCallback((tabId: string) => {
    setTabs(prev => {
      if (prev.length <= 1) return prev; // Don't close last tab
      const filtered = prev.filter(t => t.id !== tabId);
      // If closing active tab, switch to adjacent
      if (activeTabId === tabId) {
        const idx = prev.findIndex(t => t.id === tabId);
        const newIdx = Math.min(idx, filtered.length - 1);
        setActiveTabId(filtered[newIdx].id);
      }
      return filtered;
    });
  }, [activeTabId]);

  const switchTab = useCallback((tabId: string) => {
    setActiveTabId(tabId);
  }, []);

  const renameTab = useCallback((tabId: string, newName: string) => {
    setTabs(prev => prev.map(t =>
      t.id === tabId ? { ...t, filename: newName } : t
    ));
  }, []);

  // ── Source updates ──────────────────────────────────────────────────────────

  const updateSource = useCallback((tabId: string, source: string) => {
    setTabs(prev => prev.map(t =>
      t.id === tabId ? { ...t, source } : t
    ));
  }, []);

  // ── Command parsing ─────────────────────────────────────────────────────────

  const parseCommands = useCallback(async (tabId: string, source: string) => {
    if (!source.trim() || !backendUrl) return;

    // Cancel previous parse for this tab
    const existing = parseAbortRef.current.get(tabId);
    if (existing) existing.abort();

    const controller = new AbortController();
    parseAbortRef.current.set(tabId, controller);

    setTabs(prev => prev.map(t =>
      t.id === tabId ? { ...t, isParsing: true, parseError: '' } : t
    ));

    try {
      const response = await fetch(`${backendUrl}/dsl/commands`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ source }),
        signal: controller.signal,
      });

      if (!response.ok) {
        const err = await response.text();
        setTabs(prev => prev.map(t =>
          t.id === tabId ? { ...t, isParsing: false, parseError: err } : t
        ));
        return;
      }

      const data = await response.json();
      const cmds: CommandInfo[] = data.commands || [];

      setTabs(prev => prev.map(t => {
        if (t.id !== tabId) return t;

        // Preserve selected command if still valid
        let selectedCommand = t.selectedCommand;
        if (!selectedCommand || !cmds.find(c => c.name === selectedCommand)) {
          selectedCommand = cmds.length > 0 ? cmds[0].name : '';
        }

        // Build default param values for new commands
        const paramValues = { ...t.paramValues };
        for (const cmd of cmds) {
          if (!paramValues[cmd.name]) {
            const defaults: Record<string, unknown> = {};
            for (const p of cmd.params) {
              defaults[p.name] = p.default;
            }
            paramValues[cmd.name] = defaults;
          }
        }

        return {
          ...t,
          commands: cmds,
          selectedCommand,
          paramValues,
          isParsing: false,
          parseError: '',
        };
      }));
    } catch (err) {
      if (err instanceof Error && err.name === 'AbortError') return;
      setTabs(prev => prev.map(t =>
        t.id === tabId ? {
          ...t,
          isParsing: false,
          parseError: err instanceof Error ? err.message : 'Parse failed',
        } : t
      ));
    } finally {
      parseAbortRef.current.delete(tabId);
    }
  }, [backendUrl]);

  // ── Command selection ───────────────────────────────────────────────────────

  const selectCommand = useCallback((tabId: string, commandName: string) => {
    setTabs(prev => prev.map(t => {
      if (t.id !== tabId) return t;

      const cmd = t.commands.find(c => c.name === commandName);
      if (!cmd) return { ...t, selectedCommand: commandName };

      // Initialize defaults if missing
      const paramValues = { ...t.paramValues };
      if (!paramValues[commandName]) {
        const defaults: Record<string, unknown> = {};
        for (const p of cmd.params) {
          defaults[p.name] = p.default;
        }
        paramValues[commandName] = defaults;
      }

      return { ...t, selectedCommand: commandName, paramValues };
    }));
  }, []);

  // ── Parameter updates ───────────────────────────────────────────────────────

  const updateParam = useCallback((tabId: string, command: string, paramName: string, value: unknown) => {
    setTabs(prev => prev.map(t => {
      if (t.id !== tabId) return t;
      const paramValues = { ...t.paramValues };
      paramValues[command] = { ...paramValues[command], [paramName]: value };
      return { ...t, paramValues };
    }));
  }, []);

  // ── External source loading (from FileBar) ─────────────────────────────────

  const loadExternalSource = useCallback((source: string) => {
    // Load into the active tab
    setTabs(prev => prev.map(t =>
      t.id === activeTabId ? { ...t, source } : t
    ));
  }, [activeTabId]);

  return {
    tabs,
    activeTab,
    activeTabId,
    addTab,
    closeTab,
    switchTab,
    renameTab,
    updateSource,
    parseCommands,
    selectCommand,
    updateParam,
    loadExternalSource,
  };
}
