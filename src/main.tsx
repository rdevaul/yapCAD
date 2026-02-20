/**
 * yapCAD Web Viewer - Entry Point
 */

import React from 'react';
import ReactDOM from 'react-dom/client';
import { App } from './App';

// Global styles
const globalStyles = `
  * {
    box-sizing: border-box;
    margin: 0;
    padding: 0;
  }
  
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background-color: #0f0f1a;
    color: #eee;
    overflow: hidden;
  }
  
  ::-webkit-scrollbar {
    width: 8px;
    height: 8px;
  }
  
  ::-webkit-scrollbar-track {
    background: #1a1a2e;
  }
  
  ::-webkit-scrollbar-thumb {
    background: #333;
    border-radius: 4px;
  }
  
  ::-webkit-scrollbar-thumb:hover {
    background: #444;
  }

  @keyframes blink {
    0%, 100% { opacity: 1; }
    50% { opacity: 0; }
  }
`;

// Inject global styles
const styleSheet = document.createElement('style');
styleSheet.textContent = globalStyles;
document.head.appendChild(styleSheet);

// Seed chat config from env vars (only if not already set by user)
const envToken = import.meta.env.VITE_OPENCLAW_TOKEN as string | undefined;
const envGateway = import.meta.env.VITE_OPENCLAW_GATEWAY as string | undefined;
if (envToken && !localStorage.getItem('yapcad-chat-token')) {
  localStorage.setItem('yapcad-chat-token', envToken);
}
if (envGateway && !localStorage.getItem('yapcad-chat-gateway-url')) {
  localStorage.setItem('yapcad-chat-gateway-url', envGateway);
}

// Mount app
const root = ReactDOM.createRoot(document.getElementById('root')!);
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
