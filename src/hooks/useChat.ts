/**
 * Streaming chat hook for yapCAD workbench.
 *
 * Two modes:
 *   1. AGENT SESSION (default): Routes through a named OpenClaw agent session
 *      via x-openclaw-agent-id + x-openclaw-session-key headers.
 *      The session maintains its own history, memory, and tool access.
 *      Only the new user message is sent each turn.
 *
 *   2. STATELESS (legacy): Raw /v1/chat/completions with full local history.
 *      Used when sessionKey is empty/null.
 */

import { useState, useCallback, useRef } from 'react';

export interface ChatMessage {
  id: string;
  role: 'user' | 'assistant';
  content: string;
  imagePreview?: string;   // ADD: base64 data URL for display in chat history
  pending?: boolean;
  timestamp: number;
}

export const AVAILABLE_MODELS = [
  { id: 'anthropic/claude-opus-4-6',   label: 'Opus 4.6 (best)' },
  { id: 'anthropic/claude-sonnet-4-6', label: 'Sonnet 4.6 (fast)' },
  { id: 'qwen3:235b',                  label: 'Qwen3 235B (local)' },
  { id: 'qwen3-coder:480b',            label: 'Qwen3 Coder 480B (local)' },
] as const;

/**
 * Known DML team members → their OpenClaw agent ID.
 * Each agent has its own workspace, memory, and session history.
 * The agent ID determines which per-user context Jarvis loads.
 */
export const KNOWN_USERS: Record<string, string> = {
  rich:    'jarvis-rich',
  garrett: 'jarvis-garrett',
  jeremy:  'jarvis-jeremy',
  umair:   'jarvis-umair',
  yang:    'jarvis-yang',
};

export const DEFAULT_MODEL    = 'anthropic/claude-opus-4-6';
export const DEFAULT_USER     = 'rich';
export const DEFAULT_AGENT_ID = KNOWN_USERS[DEFAULT_USER];
export const DEFAULT_SESSION_KEY = `agent:${DEFAULT_AGENT_ID}:yapcad`;

/** Derive agent ID and session key from a username */
export function resolveUserSession(userName: string): { agentId: string; sessionKey: string } {
  const normalized = userName.trim().toLowerCase();
  const agentId = KNOWN_USERS[normalized] ?? `jarvis-${normalized}`;
  const sessionKey = `agent:${agentId}:yapcad`;
  return { agentId, sessionKey };
}

interface UseChatOptions {
  gatewayUrl: string;
  token: string;
  dslSource: string;
  model?: string;
  /** OpenClaw agent id to route to (e.g. "jarvis-rich"). Empty = stateless mode. */
  agentId?: string;
  /** Full session key (e.g. "agent:jarvis-rich:yapcad"). Empty = stateless mode. */
  sessionKey?: string;
}

interface UseChatReturn {
  messages: ChatMessage[];
  isStreaming: boolean;
  error: string | null;
  send: (text: string, image?: { base64: string; mimeType: string }) => Promise<void>;
  cancel: () => void;   // ADD
  clear: () => void;
}

// ── Stateless (legacy) system prompt ──────────────────────────────────────────
const STATELESS_SYSTEM_PROMPT = `You are Jarvis, an AI engineering assistant embedded in the yapCAD workbench.
You help users design 3D geometry using the yapCAD DSL language.

The yapCAD DSL supports:
- Primitives: box(w, d, h), cylinder(r, h), sphere(r), cone(r1, r2, h)
- Boolean ops: union(a, b), difference(a, b), intersection(a, b)
- Transforms: translate(solid, [x,y,z]), rotate(solid, angle, axis), scale(solid, [x,y,z])
- Revolution: lathe(profile, steps) — creates watertight revolution solid
- Assembly: assembly(name, parts), datum(name, type, pos, dir), mate(a, b, type)
- Parametric commands: command name(params...) -> type: ... emit result

When you provide DSL code, always wrap it in \`\`\`dsl fences so it can be applied to the editor.
Be concise. Focus on geometry and engineering. Suggest improvements proactively.`;

function buildStatelessSystemPrompt(dslSource: string): string {
  if (!dslSource.trim()) return STATELESS_SYSTEM_PROMPT;
  return `${STATELESS_SYSTEM_PROMPT}

Current DSL source in the editor:
\`\`\`dsl
${dslSource}
\`\`\``;
}

// ── Agent session: DSL context prepended to user turn ─────────────────────────
function buildAgentUserMessage(text: string, dslSource: string): string {
  if (!dslSource.trim()) return text;
  return `[yapCAD Editor Context]
\`\`\`dsl
${dslSource}
\`\`\`

${text}`;
}

function generateId(): string {
  return `${Date.now()}-${Math.random().toString(36).slice(2)}`;
}

export function useChat({
  gatewayUrl,
  token,
  dslSource,
  model = DEFAULT_MODEL,
  agentId = DEFAULT_AGENT_ID,
  sessionKey = DEFAULT_SESSION_KEY,
}: UseChatOptions): UseChatReturn {
  const [messages, setMessages] = useState<ChatMessage[]>([]);
  const [isStreaming, setIsStreaming] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const abortRef = useRef<AbortController | null>(null);

  const useAgentSession = !!(agentId && sessionKey);

  const send = useCallback(async (text: string, image?: { base64: string; mimeType: string }) => {
    if (!text.trim() || isStreaming) return;
    if (!token) {
      setError('No OpenClaw token configured. Add it in chat settings (⚙).');
      return;
    }

    if (abortRef.current) abortRef.current.abort();
    const controller = new AbortController();
    abortRef.current = controller;
    setError(null);

    const userMsg: ChatMessage = {
      id: generateId(),
      role: 'user',
      content: text,
      imagePreview: image ? `data:${image.mimeType};base64,${image.base64}` : undefined,
      timestamp: Date.now(),
    };
    const assistantId = generateId();
    const assistantMsg: ChatMessage = {
      id: assistantId,
      role: 'assistant',
      content: '',
      pending: true,
      timestamp: Date.now(),
    };

    setMessages(prev => [...prev, userMsg, assistantMsg]);
    setIsStreaming(true);

    // Build the messages payload
    let apiMessages: { role: string; content: string }[];

    if (useAgentSession) {
      // Agent session mode: only send the new user message.
      // DSL context is prepended to the content so the agent sees the editor state.
      // The session maintains history, memory, and tool access server-side.
      if (image) {
        apiMessages = [{
          role: 'user',
          content: [
            { type: 'image', source: { type: 'base64', media_type: image.mimeType, data: image.base64 } },
            { type: 'text', text: buildAgentUserMessage(text, dslSource) },
          ] as unknown as string,
        }];
      } else {
        apiMessages = [
          { role: 'user', content: buildAgentUserMessage(text, dslSource) },
        ];
      }
    } else {
      // Stateless mode: send full local history with system prompt.
      const contentItems: any[] = [];
      if (image) {
        contentItems.push({
          type: 'image',
          source: {
            type: 'base64',
            media_type: image.mimeType,
            data: image.base64,
          },
        });
      }
      contentItems.push({
        type: 'text',
        text: text,
      });

      apiMessages = [
        { role: 'system', content: buildStatelessSystemPrompt(dslSource) },
        ...messages.map(m => ({ role: m.role, content: m.content })),
        { role: 'user', content: JSON.stringify(contentItems) },
      ];
    }

    // Build headers
    const headers: Record<string, string> = {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${token}`,
    };

    if (useAgentSession) {
      headers['x-openclaw-agent-id']   = agentId;
      headers['x-openclaw-session-key'] = sessionKey;
    }

    try {
      const resp = await fetch(`${gatewayUrl}/v1/chat/completions`, {
        method: 'POST',
        headers,
        body: JSON.stringify({
          model,
          messages: apiMessages,
          stream: true,
          max_tokens: 4096,
        }),
        signal: controller.signal,
      });

      if (!resp.ok) {
        const errText = await resp.text();
        throw new Error(`Gateway error ${resp.status}: ${errText}`);
      }

      const reader = resp.body?.getReader();
      if (!reader) throw new Error('No response body');

      const decoder = new TextDecoder();
      let accumulated = '';
      let buffer = '';

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;

        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split('\n');
        buffer = lines.pop() ?? '';

        for (const line of lines) {
          if (!line.startsWith('data: ')) continue;
          const data = line.slice(6).trim();
          if (data === '[DONE]') break;

          try {
            const parsed = JSON.parse(data);
            const delta = parsed.choices?.[0]?.delta?.content ?? '';
            if (delta) {
              accumulated += delta;
              setMessages(prev =>
                prev.map(m =>
                  m.id === assistantId ? { ...m, content: accumulated } : m
                )
              );
            }
          } catch {
            // malformed SSE chunk, skip
          }
        }
      }

      setMessages(prev =>
        prev.map(m =>
          m.id === assistantId
            ? { ...m, content: accumulated, pending: false }
            : m
        )
      );
    } catch (err: unknown) {
      if (err instanceof Error && err.name === 'AbortError') return;
      const msg = err instanceof Error ? err.message : 'Unknown error';
      setError(msg);
      setMessages(prev => prev.filter(m => m.id !== assistantId));
    } finally {
      setIsStreaming(false);
      abortRef.current = null;
    }
  }, [messages, isStreaming, token, dslSource, gatewayUrl, model, agentId, sessionKey, useAgentSession]);

  const cancel = useCallback(() => {
    abortRef.current?.abort();
  }, []);

  const clear = useCallback(() => {
    setMessages([]);
    setError(null);
  }, []);

  return { messages, isStreaming, error, send, cancel, clear };
}