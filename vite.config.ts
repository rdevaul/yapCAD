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
    },
  },
  build: {
    outDir: 'dist',
    sourcemap: true,
  },
});
