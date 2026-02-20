# yapCAD DSL Editor Integration

This document describes the new DSL code editor panel added to the yapCAD Web Viewer.

## Overview

The yapCAD Web Viewer now includes an integrated Monaco Editor panel that allows users to write and evaluate yapCAD DSL code directly in the browser. The evaluated geometry is loaded into the 3D viewer in real-time.

## Features Added

### 1. DSL Editor Panel
- **Location**: Right side of the interface (40% width, collapsible)
- **Editor**: Monaco Editor with yapCAD DSL syntax highlighting
- **Functionality**: 
  - Real-time syntax highlighting for yapCAD keywords and built-ins
  - Auto-completion for DSL functions (cube, sphere, union, etc.)
  - Auto-save to localStorage (1-second debounce)
  - Keyboard shortcut: Ctrl+Enter to evaluate
  - Error display with detailed messages

### 2. Chat Panel Placeholder  
- **Location**: Bottom half of right panel (resizable split)
- **Status**: Placeholder for Phase 3 implementation
- **Purpose**: Will provide AI assistance for DSL development

### 3. Session Management
- **Connection Status**: Green/red indicator showing backend connectivity
- **Backend URL**: Editable backend URL (saved to localStorage)
- **Default Backend**: http://localhost:8000
- **Authentication**: Placeholder for Phase 2 token-based auth

### 4. Layout Changes
- **Left Panel**: 3D viewport (60% width when editor is open)
- **Right Panel**: DSL editor + chat (40% width, collapsible) 
- **Resizable Split**: Drag handle between editor and chat panels
- **Responsive**: Smooth transitions when toggling panels

## Usage

### Writing DSL Code
1. Click "Show Editor" to open the DSL panel
2. Write yapCAD DSL code in the Monaco editor
3. Use auto-completion (Ctrl+Space) for built-in functions
4. Click "▶ Evaluate" or press Ctrl+Enter to run the code

### Example DSL Code
```yapcad-dsl
# Create a cube
let my_cube = cube(10, 10, 10);

# Create a sphere  
let my_sphere = sphere(5);

# Combine with union
let combined = union(my_cube, my_sphere);

# Emit the result
emit combined;
```

### Backend Communication
- **Endpoint**: POST `/dsl/eval`
- **Payload**: `{ "source": "<dsl_code>", "format": "json" }`
- **Response**: yapCAD geometry JSON format
- **Timeout**: 30 seconds for evaluation

## Technical Implementation

### Components Added
- `DSLEditor.tsx`: Monaco-based editor with yapCAD language support
- `ChatPanel.tsx`: Placeholder panel for future chat functionality
- `SessionBar.tsx`: Backend connection and authentication status
- `ResizableSplit.tsx`: Resizable panel splitter component

### Integration Points
- **Geometry Loading**: Reuses existing `yapcad-loader.ts` logic
- **Viewer Integration**: Loads evaluated geometry into YapCADViewer
- **State Management**: New state for editor, backend URL, evaluation status
- **Storage**: localStorage for editor content, backend URL, panel sizes

### Styling
- **Theme**: Consistent dark theme matching the existing viewer  
- **Colors**: #1a1a2e background, #252540 panels, #4c9f38 accents
- **Layout**: All styles are inline (no CSS files)

## Future Phases

### Phase 2: Authentication
- Token-based authentication system
- Secure backend communication
- User session management

### Phase 3: AI Chat
- WebSocket-based chat with AI assistant
- DSL syntax help and code generation
- Interactive debugging and suggestions
- Context-aware geometry design guidance

## Verification

The implementation has been tested with:
- ✅ Successful TypeScript compilation (`npm run build`)
- ✅ Development server startup (`npm run dev`)  
- ✅ Monaco Editor integration
- ✅ Component rendering and layout
- ✅ localStorage persistence
- ✅ Resizable panels

## Backend Requirements

For full functionality, the backend service must implement:
- `GET /health` - Health check endpoint
- `POST /dsl/eval` - DSL evaluation endpoint

The backend should accept DSL source code and return yapCAD geometry JSON format compatible with the existing loader.