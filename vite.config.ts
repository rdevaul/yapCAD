import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { localPackagesPlugin } from './vite-plugin-local-packages';

export default defineConfig({
  plugins: [
    react(),
    localPackagesPlugin(),
  ],
  server: {
    port: 5173,
    host: '0.0.0.0',
    proxy: {
      // Proxy OpenClaw gateway to avoid CORS issues in browser
      '/openclaw': {
        target: 'http://localhost:18789',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/openclaw/, ''),
      },
      // Proxy yapCAD backend WebSocket (DSL eval)
      '/ws': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        ws: true,
      },
      // Proxy yapCAD backend REST API
      '/dsl': {
        target: 'http://localhost:8000',
        changeOrigin: true,
      },
      // Proxy WBS mission control — allows remote clients (e.g. Garrett) to
      // reach the WBS server via the Vite host without needing direct access.
      // wbs:// paths in FileBar are fetched as /wbs/api/dsl/<rel>.
      // Agentic-1 runs on 8765; Hydra runs on 8766.
      // Change target here to match whichever project is active.
      '/wbs': {
        target: 'http://localhost:8766',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/wbs/, ''),
      },
    },
  },
  build: {
    outDir: 'dist',
    sourcemap: true,
  },
});
