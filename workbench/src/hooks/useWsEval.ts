/**
 * WebSocket eval hook for yapCAD workbench.
 *
 * Maintains a single WebSocket connection to the DSL eval endpoint.
 * Falls back to REST `/dsl/eval` if WebSocket is unavailable.
 *
 * Features:
 *   - Request ID generation and in-flight cancellation
 *   - Latency tier tracking per command
 *   - Graceful reconnection
 */

import { useRef, useCallback, useEffect, useState } from 'react';

// ── Types ─────────────────────────────────────────────────────────────────────

export interface EvalResult {
  success: boolean;
  geometry: Record<string, unknown> | null;
  scalar_result: unknown;
  error_message: string | null;
  volume: number | null;
  elapsed_ms: number;
}

export type LatencyTier = 'fast' | 'medium' | 'slow';

interface PendingEval {
  requestId: string;
  resolve: (result: EvalResult) => void;
  reject: (err: Error) => void;
}

interface WsEvalMessage {
  type: 'eval';
  request_id: string;
  source: string;
  command: string;
  parameters: Record<string, unknown>;
  format: string;
}

interface WsCancelMessage {
  type: 'cancel';
  request_id: string;
}

interface WsResultMessage {
  type: 'result';
  request_id: string;
  success: boolean;
  geometry?: Record<string, unknown>;
  scalar_result?: unknown;
  error_message?: string | null;
  volume?: number | null;
  elapsed_ms?: number;
}

// ── Hook ──────────────────────────────────────────────────────────────────────

interface UseWsEvalOptions {
  /** Base URL for the backend (REST fallback). e.g. "http://localhost:8000" or "" for relative */
  backendUrl: string;
  /** WebSocket URL. If empty, derived from backendUrl. */
  wsUrl?: string;
}

export function useWsEval({ backendUrl, wsUrl }: UseWsEvalOptions) {
  const wsRef = useRef<WebSocket | null>(null);
  const pendingRef = useRef<PendingEval | null>(null);
  const latencyMapRef = useRef<Map<string, number>>(new Map());
  const [wsConnected, setWsConnected] = useState(false);
  const reconnectTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const mountedRef = useRef(true);

  // Derive WS URL from backend URL if not provided
  const resolvedWsUrl = wsUrl || (() => {
    // If backendUrl is relative (e.g. "/api"), use window.location
    if (backendUrl.startsWith('/') || backendUrl === '') {
      const proto = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
      return `${proto}//${window.location.host}/ws/dsl/eval`;
    }
    // Absolute URL — convert http(s) to ws(s)
    const url = new URL(backendUrl);
    const proto = url.protocol === 'https:' ? 'wss:' : 'ws:';
    return `${proto}//${url.host}/ws/dsl/eval`;
  })();

  const connectWs = useCallback(() => {
    if (wsRef.current?.readyState === WebSocket.OPEN || wsRef.current?.readyState === WebSocket.CONNECTING) {
      return;
    }

    try {
      const ws = new WebSocket(resolvedWsUrl);

      ws.onopen = () => {
        if (!mountedRef.current) { ws.close(); return; }
        setWsConnected(true);
        console.log('[useWsEval] WebSocket connected');
      };

      ws.onmessage = (event) => {
        try {
          const msg = JSON.parse(event.data) as WsResultMessage;
          if (msg.type === 'result' && pendingRef.current?.requestId === msg.request_id) {
            const pending = pendingRef.current;
            pendingRef.current = null;

            const result: EvalResult = {
              success: msg.success,
              geometry: msg.geometry ?? null,
              scalar_result: msg.scalar_result ?? null,
              error_message: msg.error_message ?? null,
              volume: msg.volume ?? null,
              elapsed_ms: msg.elapsed_ms ?? 0,
            };

            pending.resolve(result);
          }
        } catch {
          // Malformed message, ignore
        }
      };

      ws.onclose = () => {
        if (!mountedRef.current) return;
        setWsConnected(false);
        wsRef.current = null;
        // Reconnect after delay
        reconnectTimerRef.current = setTimeout(() => {
          if (mountedRef.current) connectWs();
        }, 3000);
      };

      ws.onerror = () => {
        // onclose will fire after onerror, which handles reconnect
        ws.close();
      };

      wsRef.current = ws;
    } catch {
      // WebSocket constructor failed — will fallback to REST
      console.warn('[useWsEval] WebSocket not available, using REST fallback');
    }
  }, [resolvedWsUrl]);

  // Connect on mount
  useEffect(() => {
    mountedRef.current = true;
    connectWs();
    return () => {
      mountedRef.current = false;
      if (reconnectTimerRef.current) clearTimeout(reconnectTimerRef.current);
      if (wsRef.current) {
        wsRef.current.onclose = null; // prevent reconnect on unmount
        wsRef.current.close();
        wsRef.current = null;
      }
      // Reject any pending eval
      if (pendingRef.current) {
        pendingRef.current.reject(new Error('Component unmounted'));
        pendingRef.current = null;
      }
    };
  }, [connectWs]);

  // Generate unique request IDs
  const nextRequestId = useCallback((): string => {
    return `eval_${Date.now()}_${Math.random().toString(36).slice(2, 8)}`;
  }, []);

  // Cancel in-flight eval
  const cancelPending = useCallback(() => {
    if (!pendingRef.current) return;

    const { requestId } = pendingRef.current;

    // Send cancel via WS if connected
    if (wsRef.current?.readyState === WebSocket.OPEN) {
      const cancelMsg: WsCancelMessage = { type: 'cancel', request_id: requestId };
      wsRef.current.send(JSON.stringify(cancelMsg));
    }

    pendingRef.current.reject(new Error('Cancelled'));
    pendingRef.current = null;
  }, []);

  // REST fallback eval
  const restEval = useCallback(async (
    source: string,
    command: string,
    parameters: Record<string, unknown>,
    _requestId: string
  ): Promise<EvalResult> => {
    const startTime = performance.now();

    const response = await fetch(`${backendUrl}/dsl/eval`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        source,
        command,
        parameters,
        format: 'json',
      }),
      signal: typeof AbortSignal.timeout === 'function' ? AbortSignal.timeout(300000) : undefined,
    });

    const elapsed = performance.now() - startTime;

    if (!response.ok) {
      const errorText = await response.text();
      return {
        success: false,
        geometry: null,
        scalar_result: null,
        error_message: `Server error (${response.status}): ${errorText}`,
        volume: null,
        elapsed_ms: elapsed,
      };
    }

    const data = await response.json();
    return {
      success: data.success ?? false,
      geometry: data.geometry ?? null,
      scalar_result: data.scalar_result ?? null,
      error_message: data.error_message ?? null,
      volume: data.volume ?? null,
      elapsed_ms: elapsed,
    };
  }, [backendUrl]);

  // Main eval function
  const evalDsl = useCallback(async (
    source: string,
    command: string,
    parameters: Record<string, unknown>
  ): Promise<EvalResult> => {
    // Cancel any in-flight eval
    cancelPending();

    const requestId = nextRequestId();

    // Try WebSocket first
    if (wsRef.current?.readyState === WebSocket.OPEN) {
      return new Promise<EvalResult>((resolve, reject) => {
        pendingRef.current = { requestId, resolve, reject };

        const msg: WsEvalMessage = {
          type: 'eval',
          request_id: requestId,
          source,
          command,
          parameters,
          format: 'json',
        };

        wsRef.current!.send(JSON.stringify(msg));

        // Timeout fallback — if WS doesn't respond in 300s, reject
        setTimeout(() => {
          if (pendingRef.current?.requestId === requestId) {
            pendingRef.current.reject(new Error('WebSocket eval timeout'));
            pendingRef.current = null;
          }
        }, 300000);
      });
    }

    // REST fallback
    return restEval(source, command, parameters, requestId);
  }, [cancelPending, nextRequestId, restEval]);

  // Track latency and expose tier
  const recordLatency = useCallback((command: string, elapsedMs: number) => {
    // Only record first call per command (as per spec)
    if (!latencyMapRef.current.has(command)) {
      latencyMapRef.current.set(command, elapsedMs);
    }
  }, []);

  const getLatencyTier = useCallback((command: string): LatencyTier => {
    const ms = latencyMapRef.current.get(command);
    if (ms === undefined) return 'medium'; // default before first measurement
    if (ms < 50) return 'fast';
    if (ms <= 500) return 'medium';
    return 'slow';
  }, []);

  return {
    eval: evalDsl,
    cancelPending,
    getLatencyTier,
    recordLatency,
    wsConnected,
  };
}
