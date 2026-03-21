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
    widget?: string;
    label?: string;
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

  // ── Undo stacks — refs so mutations don't cause re-renders ─────────────────
  // paramHistoryRef: tabId → stack of paramValues snapshots
  const paramHistoryRef = useRef<Map<string, Array<Record<string, Record<string, unknown>>>>>(new Map());
  // closedTabsRef: stack of recently closed TabState objects
  const closedTabsRef = useRef<TabState[]>([]);
  const MAX_HISTORY = 20;

  // undoCount: reactive counter so canUndo() triggers re-renders when stacks change
  const [undoCount, setUndoCount] = useState(0);

  // Persist on change
  useEffect(() => {
    saveTabs(tabs, activeTabId);
  }, [tabs, activeTabId]);

  // Active tab helper
  const activeTab = tabs.find(t => t.id === activeTabId) ?? tabs[0];

  // ── Undo helpers ────────────────────────────────────────────────────────────

  const pushParamHistory = useCallback((tabId: string, snapshot: Record<string, Record<string, unknown>>) => {
    const stack = paramHistoryRef.current.get(tabId) ?? [];
    stack.push(JSON.parse(JSON.stringify(snapshot))); // deep clone
    if (stack.length > MAX_HISTORY) stack.shift();
    paramHistoryRef.current.set(tabId, stack);
    setUndoCount(n => n + 1);
  }, []);

  // Public: let App.tsx push a snapshot before handleSetDefaults rewrites source
  const pushParamSnapshotForTab = useCallback((tabId: string) => {
    const tab = tabs.find(t => t.id === tabId);
    if (tab) pushParamHistory(tabId, tab.paramValues);
  }, [tabs, pushParamHistory]);

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
      const closing = prev.find(t => t.id === tabId);
      if (closing) {
        const stack = closedTabsRef.current;
        stack.push(closing);
        if (stack.length > 10) stack.shift();
        setUndoCount(n => n + 1);
      }
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
    setTabs(prev => {
      const currentTab = prev.find(t => t.id === tabId);
      if (currentTab) pushParamHistory(tabId, currentTab.paramValues);
      return prev.map(t => {
        if (t.id !== tabId) return t;
        const paramValues = { ...t.paramValues };
        paramValues[command] = { ...paramValues[command], [paramName]: value };
        for (const cmd of t.commands) {
          if (cmd.name === command) continue;
          if (cmd.params.some(p => p.name === paramName)) {
            paramValues[cmd.name] = { ...paramValues[cmd.name], [paramName]: value };
          }
        }
        return { ...t, paramValues };
      });
    });
  }, [pushParamHistory]);

  // Batch-update multiple params in one setState call — use this when updating
  // several params at once (e.g. all control points from a drag) to avoid
  // React batching issues where each call reads stale state.
  const batchUpdateParams = useCallback((tabId: string, command: string, updates: Record<string, unknown>) => {
    setTabs(prev => {
      const currentTab = prev.find(t => t.id === tabId);
      if (currentTab) pushParamHistory(tabId, currentTab.paramValues);
      return prev.map(t => {
        if (t.id !== tabId) return t;
        const paramValues = { ...t.paramValues };
        // Apply all updates to the target command at once
        paramValues[command] = { ...paramValues[command], ...updates };
        // Sync each updated param to sibling commands
        for (const [paramName, value] of Object.entries(updates)) {
          for (const cmd of t.commands) {
            if (cmd.name === command) continue;
            if (cmd.params.some(p => p.name === paramName)) {
              paramValues[cmd.name] = { ...paramValues[cmd.name], [paramName]: value };
            }
          }
        }
        return { ...t, paramValues };
      });
    });
  }, [pushParamHistory]);

  // ── Reset params to DSL defaults ───────────────────────────────────────────

  const resetParamDefaults = useCallback((tabId: string, command: string) => {
    setTabs(prev => {
      const currentTab = prev.find(t => t.id === tabId);
      if (currentTab) pushParamHistory(tabId, currentTab.paramValues);
      return prev.map(t => {
        if (t.id !== tabId) return t;
        const cmdDef = t.commands.find(c => c.name === command);
        if (!cmdDef) return t;
        const defaults: Record<string, unknown> = {};
        for (const p of cmdDef.params) {
          defaults[p.name] = p.default;
        }
        return {
          ...t,
          paramValues: {
            ...t.paramValues,
            [command]: defaults,
          },
        };
      });
    });
  }, [pushParamHistory]);

  // ── Undo ────────────────────────────────────────────────────────────────────

  const undoLastAction = useCallback((): boolean => {
    // First priority: undo param change on active tab
    if (activeTabId) {
      const stack = paramHistoryRef.current.get(activeTabId) ?? [];
      if (stack.length > 0) {
        const prev = stack.pop()!;
        setTabs(t => t.map(tab =>
          tab.id === activeTabId ? { ...tab, paramValues: prev } : tab
        ));
        // Keep undoCount accurate — decrement since we popped
        setUndoCount(n => Math.max(0, n - 1));
        return true;
      }
    }
    // Second priority: reopen last closed tab
    const stack = closedTabsRef.current;
    if (stack.length > 0) {
      const restored = stack.pop()!;
      setTabs(prev => [...prev, restored]);
      setActiveTabId(restored.id);
      setUndoCount(n => Math.max(0, n - 1));
      return true;
    }
    return false;
  }, [activeTabId]);

  const canUndo = useCallback((): boolean => {
    if (undoCount <= 0) return false;
    const paramStack = activeTabId ? (paramHistoryRef.current.get(activeTabId) ?? []) : [];
    return paramStack.length > 0 || closedTabsRef.current.length > 0;
  }, [activeTabId, undoCount]);

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
    batchUpdateParams,
    resetParamDefaults,
    loadExternalSource,
    undoLastAction,
    canUndo,
    pushParamSnapshotForTab,
    undoCount,
  };
}
