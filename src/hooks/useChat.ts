/**
 * Streaming chat hook for yapCAD workbench.
 * Connects to the OpenClaw gateway's /v1/chat/completions SSE endpoint.
 * DSL source is injected as context in the system prompt.
 */

import { useState, useCallback, useRef } from 'react';

export interface ChatMessage {
  id: string;
  role: 'user' | 'assistant';
  content: string;
  pending?: boolean;
  timestamp: number;
}

export const AVAILABLE_MODELS = [
  { id: 'anthropic/claude-opus-4-6',    label: 'Opus 4.6 (best)' },
  { id: 'anthropic/claude-sonnet-4-6',  label: 'Sonnet 4.6 (fast)' },
  { id: 'qwen3:235b',                   label: 'Qwen3 235B (local)' },
  { id: 'qwen3-coder:480b',             label: 'Qwen3 Coder 480B (local)' },
] as const;

export const DEFAULT_MODEL = 'anthropic/claude-opus-4-6';

interface UseChatOptions {
  gatewayUrl: string;
  token: string;
  dslSource: string;
  model?: string;
}

interface UseChatReturn {
  messages: ChatMessage[];
  isStreaming: boolean;
  error: string | null;
  send: (text: string) => Promise<void>;
  clear: () => void;
}

const SYSTEM_PROMPT = `You are Jarvis, an AI engineering assistant embedded in the yapCAD workbench.
You help users design 3D geometry using the yapCAD DSL language.

The yapCAD DSL supports:
- Primitives: box(w, d, h), cylinder(r, h), sphere(r), cone(r1, r2, h)
- Boolean ops: union(a, b), difference(a, b), intersection(a, b)
- Transforms: translate(solid, [x,y,z]), rotate(solid, angle, axis), scale(solid, [x,y,z])
- Revolution: lathe(profile, steps) — creates watertight revolution solid
- Assembly: assembly(name, parts), datum(name, type, pos, dir), mate(a, b, type)
- Variables and parametric commands: command name(params...) { ... }

When you provide DSL code, always wrap it in \`\`\`dsl fences so it can be applied to the editor.
Be concise. Focus on geometry and engineering. Suggest improvements proactively.`;

function buildSystemPrompt(dslSource: string): string {
  if (!dslSource.trim()) {
    return SYSTEM_PROMPT;
  }
  return `${SYSTEM_PROMPT}

Current DSL source in the editor:
\`\`\`dsl
${dslSource}
\`\`\`

You can suggest modifications or improvements to this code. The user can apply DSL blocks you provide.`;
}

function generateId(): string {
  return `${Date.now()}-${Math.random().toString(36).slice(2)}`;
}

export function useChat({ gatewayUrl, token, dslSource, model = DEFAULT_MODEL }: UseChatOptions): UseChatReturn {
  const [messages, setMessages] = useState<ChatMessage[]>([]);
  const [isStreaming, setIsStreaming] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const abortRef = useRef<AbortController | null>(null);

  const send = useCallback(async (text: string) => {
    if (!text.trim() || isStreaming) return;
    if (!token) {
      setError('No OpenClaw token configured. Add it in chat settings.');
      return;
    }

    // Abort any in-flight request
    if (abortRef.current) {
      abortRef.current.abort();
    }
    const controller = new AbortController();
    abortRef.current = controller;

    setError(null);

    // Add user message
    const userMsg: ChatMessage = {
      id: generateId(),
      role: 'user',
      content: text,
      timestamp: Date.now(),
    };

    // Add pending assistant message
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

    // Build messages array for API
    const apiMessages = [
      { role: 'system', content: buildSystemPrompt(dslSource) },
      ...messages.map(m => ({ role: m.role, content: m.content })),
      { role: 'user', content: text },
    ];

    try {
      const resp = await fetch(`${gatewayUrl}/v1/chat/completions`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`,
        },
        body: JSON.stringify({
          model,
          messages: apiMessages,
          stream: true,
          max_tokens: 2048,
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
                  m.id === assistantId
                    ? { ...m, content: accumulated }
                    : m
                )
              );
            }
          } catch {
            // malformed chunk, skip
          }
        }
      }

      // Finalize — clear pending flag
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
      // Remove the pending assistant message on error
      setMessages(prev => prev.filter(m => m.id !== assistantId));
    } finally {
      setIsStreaming(false);
      abortRef.current = null;
    }
  }, [messages, isStreaming, token, dslSource, gatewayUrl, model]);

  const clear = useCallback(() => {
    setMessages([]);
    setError(null);
  }, []);

  return { messages, isStreaming, error, send, clear };
}
