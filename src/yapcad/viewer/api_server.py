"""
REST API and WebSocket Server
=============================

Flask-based REST API and SocketIO WebSocket server for remote viewer control.

This module provides a complete HTTP API for controlling the VTK viewer
programmatically, along with WebSocket support for real-time event streaming.

Example Usage
-------------

Basic server startup::

    from yapcad.viewer import VTKViewer, ViewerAPIServer, ViewerConfig

    config = ViewerConfig(stl_dir="/path/to/stls")
    viewer = VTKViewer(config)
    server = ViewerAPIServer(viewer)
    server.run(port=5000)

With custom configuration::

    server = ViewerAPIServer(viewer, cors_origins=["http://localhost:3000"])
    server.run(host="0.0.0.0", port=8080, debug=True)

REST API Endpoints
------------------

**Status and Information:**

- ``GET /status`` - Get viewer status
- ``GET /parts`` - List loaded parts
- ``GET /config`` - Get viewer configuration

**Part Management:**

- ``POST /load`` - Load parts from JSON or part list
- ``POST /reload`` - Reload positions from JSON file
- ``POST /clear`` - Clear all parts

**Rendering Controls:**

- ``POST /xray`` - Toggle or set x-ray mode
- ``POST /highlight`` - Highlight specific parts

**Camera Controls:**

- ``POST /camera/reset`` - Reset camera to default
- ``POST /camera/set`` - Set camera position
- ``POST /camera/orbit`` - Orbit camera around focal point
- ``POST /camera/focus`` - Focus on a specific part

**Screenshot:**

- ``GET /screenshot`` - Capture and return screenshot
- ``POST /screenshot`` - Capture and save screenshot

WebSocket Events
----------------

**Server-to-Client:**

- ``status`` - Viewer status updates
- ``part_loaded`` - Part was loaded
- ``part_removed`` - Part was removed
- ``parts_cleared`` - All parts cleared
- ``xray_changed`` - X-ray mode changed
- ``camera_changed`` - Camera position changed
- ``screenshot_saved`` - Screenshot captured
- ``selection_changed`` - Selection/highlight changed
- ``error`` - Error occurred

**Client-to-Server:**

- ``command`` - Execute a viewer command
- ``subscribe`` - Subscribe to event types
- ``unsubscribe`` - Unsubscribe from event types

Dependencies
------------

Requires Flask and related packages::

    pip install flask flask-socketio flask-cors
"""

import io
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from flask import Flask, jsonify, request, send_file
from flask_cors import CORS
from flask_socketio import SocketIO, emit

from .config import ViewerConfig
from .vtk_viewer import VTKViewer
from .events import (
    EventType,
    CommandMessage,
    CommandType,
    StatusEvent,
    ErrorEvent,
)


class ViewerAPIServer:
    """
    REST API and WebSocket server for the VTK viewer.

    This class wraps a VTKViewer instance and provides HTTP endpoints
    and WebSocket events for remote control and monitoring.

    Parameters
    ----------
    viewer : VTKViewer
        The viewer instance to control.
    cors_origins : list of str, optional
        Allowed CORS origins. Default: ["*"] (all origins).
    async_mode : str, optional
        SocketIO async mode. Default: "threading".

    Attributes
    ----------
    viewer : VTKViewer
        The controlled viewer instance.
    app : Flask
        The Flask application.
    socketio : SocketIO
        The SocketIO instance.

    Examples
    --------
    Basic usage::

        viewer = VTKViewer(config)
        server = ViewerAPIServer(viewer)
        server.run(port=5000)

    With pre-configured Flask app::

        app = Flask(__name__)
        viewer = VTKViewer(config)
        server = ViewerAPIServer(viewer)
        # Add custom routes to server.app
        server.run()
    """

    def __init__(
        self,
        viewer: VTKViewer,
        cors_origins: Optional[List[str]] = None,
        async_mode: str = "threading"
    ):
        """
        Initialize the API server.

        Parameters
        ----------
        viewer : VTKViewer
            Viewer instance to control.
        cors_origins : list of str, optional
            CORS allowed origins.
        async_mode : str, optional
            SocketIO async mode.
        """
        self.viewer = viewer
        self._subscriptions: Dict[str, Set[str]] = {}  # sid -> set of event types
        self._viewer_thread: Optional[threading.Thread] = None
        self._running = False

        # Create Flask app
        self.app = Flask(__name__)

        # Enable CORS
        if cors_origins is None:
            cors_origins = ["*"]
        CORS(self.app, origins=cors_origins)

        # Create SocketIO
        self.socketio = SocketIO(
            self.app,
            cors_allowed_origins=cors_origins,
            async_mode=async_mode
        )

        # Register routes
        self._register_routes()
        self._register_socket_events()

        # Connect viewer events to WebSocket broadcasts
        self._connect_viewer_events()

    def _register_routes(self):
        """Register all REST API routes."""

        # =====================================================================
        # Status and Information
        # =====================================================================

        @self.app.route("/status", methods=["GET"])
        def get_status():
            """
            Get viewer status.

            ---
            tags:
              - Status
            responses:
              200:
                description: Current viewer status
                content:
                  application/json:
                    schema:
                      type: object
                      properties:
                        part_count:
                          type: integer
                          description: Number of loaded parts
                        xray_enabled:
                          type: boolean
                          description: Whether x-ray mode is active
                        viewport:
                          type: string
                          description: Current viewport name
                        camera_position:
                          type: array
                          items:
                            type: number
                          description: Camera position [x, y, z]
                        camera_focal_point:
                          type: array
                          items:
                            type: number
                          description: Camera focal point [x, y, z]
            """
            status = self.viewer.get_status()
            return jsonify(status.to_dict())

        @self.app.route("/parts", methods=["GET"])
        def get_parts():
            """
            List loaded parts.

            ---
            tags:
              - Parts
            responses:
              200:
                description: List of loaded parts with bounds
                content:
                  application/json:
                    schema:
                      type: object
                      properties:
                        parts:
                          type: array
                          items:
                            type: object
                            properties:
                              name:
                                type: string
                              bounds:
                                type: array
                                items:
                                  type: number
            """
            parts = []
            for name in self.viewer.get_part_names():
                bounds = self.viewer.get_part_bounds(name)
                parts.append({
                    "name": name,
                    "bounds": list(bounds) if bounds else None
                })
            return jsonify({"parts": parts, "count": len(parts)})

        @self.app.route("/config", methods=["GET"])
        def get_config():
            """
            Get viewer configuration.

            ---
            tags:
              - Configuration
            responses:
              200:
                description: Current viewer configuration
            """
            return jsonify(self.viewer.config.to_dict())

        # =====================================================================
        # Part Management
        # =====================================================================

        @self.app.route("/load", methods=["POST"])
        def load_parts():
            """
            Load parts into the viewer.

            ---
            tags:
              - Parts
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    properties:
                      parts:
                        type: array
                        items:
                          type: array
                          description: "[stl_file, part_name, color_name]"
                        description: List of part specifications
                      positions_file:
                        type: string
                        description: Path to positions JSON file (alternative to parts)
                      clear:
                        type: boolean
                        default: true
                        description: Clear existing parts first
            responses:
              200:
                description: Parts loaded successfully
                content:
                  application/json:
                    schema:
                      type: object
                      properties:
                        loaded:
                          type: integer
                          description: Number of parts loaded
            """
            data = request.get_json() or {}
            clear = data.get("clear", True)

            try:
                if "positions_file" in data:
                    loaded = self.viewer.load_from_json(data["positions_file"])
                elif "parts" in data:
                    loaded = self.viewer.load_parts(data["parts"], clear=clear)
                else:
                    # Try to reload from configured positions file
                    loaded = self.viewer.load_from_json()

                return jsonify({"loaded": loaded, "success": True})
            except Exception as e:
                return jsonify({"error": str(e), "success": False}), 400

        @self.app.route("/reload", methods=["POST"])
        def reload_positions():
            """
            Reload positions from JSON file.

            Updates transforms without reloading STL files.

            ---
            tags:
              - Parts
            responses:
              200:
                description: Positions reloaded successfully
            """
            try:
                self.viewer.reload_positions()
                return jsonify({"success": True})
            except Exception as e:
                return jsonify({"error": str(e), "success": False}), 400

        @self.app.route("/clear", methods=["POST"])
        def clear_parts():
            """
            Clear all loaded parts.

            ---
            tags:
              - Parts
            responses:
              200:
                description: Parts cleared successfully
            """
            self.viewer.clear_parts()
            return jsonify({"success": True})

        # =====================================================================
        # Rendering Controls
        # =====================================================================

        @self.app.route("/xray", methods=["POST"])
        def set_xray():
            """
            Toggle or set x-ray mode.

            ---
            tags:
              - Rendering
            requestBody:
              content:
                application/json:
                  schema:
                    type: object
                    properties:
                      enabled:
                        type: boolean
                        description: Set x-ray state (omit to toggle)
                      opacity:
                        type: number
                        description: X-ray opacity (0-1)
            responses:
              200:
                description: X-ray mode updated
                content:
                  application/json:
                    schema:
                      type: object
                      properties:
                        enabled:
                          type: boolean
                        opacity:
                          type: number
            """
            data = request.get_json() or {}

            if "enabled" in data:
                self.viewer.set_xray(
                    data["enabled"],
                    opacity=data.get("opacity")
                )
            else:
                self.viewer.toggle_xray()

            return jsonify({
                "enabled": self.viewer.xray_mode,
                "opacity": (
                    self.viewer.config.xray_opacity
                    if self.viewer.xray_mode
                    else self.viewer.config.default_opacity
                )
            })

        @self.app.route("/highlight", methods=["POST"])
        def highlight_parts():
            """
            Highlight specific parts.

            ---
            tags:
              - Rendering
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    properties:
                      parts:
                        type: array
                        items:
                          type: string
                        description: Part names to highlight
                      dim_opacity:
                        type: number
                        default: 0.15
                        description: Opacity for non-highlighted parts
            responses:
              200:
                description: Parts highlighted
            """
            data = request.get_json() or {}
            parts = data.get("parts", [])
            dim_opacity = data.get("dim_opacity", 0.15)

            self.viewer.highlight_parts(parts, dim_opacity)
            return jsonify({"highlighted": parts})

        # =====================================================================
        # Camera Controls
        # =====================================================================

        @self.app.route("/camera/reset", methods=["POST"])
        def reset_camera():
            """
            Reset camera to default position.

            ---
            tags:
              - Camera
            responses:
              200:
                description: Camera reset successfully
            """
            self.viewer.reset_camera()
            return jsonify({"success": True})

        @self.app.route("/camera/set", methods=["POST"])
        def set_camera():
            """
            Set camera position and focal point.

            ---
            tags:
              - Camera
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    required:
                      - position
                      - focal_point
                    properties:
                      position:
                        type: array
                        items:
                          type: number
                        description: Camera position [x, y, z]
                      focal_point:
                        type: array
                        items:
                          type: number
                        description: Focal point [x, y, z]
                      viewport:
                        type: string
                        default: ISO
                        description: Viewport to modify
                      view_up:
                        type: array
                        items:
                          type: number
                        default: [0, 0, 1]
                        description: Up direction vector
            responses:
              200:
                description: Camera position set
            """
            data = request.get_json() or {}

            try:
                self.viewer.set_camera(
                    position=tuple(data["position"]),
                    focal_point=tuple(data["focal_point"]),
                    viewport=data.get("viewport", "ISO"),
                    view_up=tuple(data.get("view_up", [0, 0, 1]))
                )
                return jsonify({"success": True})
            except KeyError as e:
                return jsonify({"error": f"Missing required field: {e}"}), 400
            except Exception as e:
                return jsonify({"error": str(e)}), 400

        @self.app.route("/camera/orbit", methods=["POST"])
        def orbit_camera():
            """
            Orbit camera using spherical coordinates.

            ---
            tags:
              - Camera
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    required:
                      - azimuth
                      - elevation
                      - distance
                    properties:
                      azimuth:
                        type: number
                        description: Azimuth angle in degrees
                      elevation:
                        type: number
                        description: Elevation angle in degrees
                      distance:
                        type: number
                        description: Distance from focal point
                      viewport:
                        type: string
                        default: ISO
                        description: Viewport to modify
            responses:
              200:
                description: Camera orbited successfully
            """
            data = request.get_json() or {}

            try:
                self.viewer.orbit_camera(
                    azimuth=data["azimuth"],
                    elevation=data["elevation"],
                    distance=data["distance"],
                    viewport=data.get("viewport", "ISO")
                )
                return jsonify({"success": True})
            except KeyError as e:
                return jsonify({"error": f"Missing required field: {e}"}), 400
            except Exception as e:
                return jsonify({"error": str(e)}), 400

        @self.app.route("/camera/focus", methods=["POST"])
        def focus_camera():
            """
            Focus camera on a specific part.

            ---
            tags:
              - Camera
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    required:
                      - part
                    properties:
                      part:
                        type: string
                        description: Part name to focus on
            responses:
              200:
                description: Camera focused on part
              404:
                description: Part not found
            """
            data = request.get_json() or {}

            try:
                self.viewer.focus_part(data["part"])
                return jsonify({"success": True, "focused_on": data["part"]})
            except KeyError:
                return jsonify({"error": f"Part not found: {data.get('part')}"}), 404
            except Exception as e:
                return jsonify({"error": str(e)}), 400

        @self.app.route("/camera", methods=["GET"])
        def get_camera():
            """
            Get current camera state.

            ---
            tags:
              - Camera
            parameters:
              - name: viewport
                in: query
                schema:
                  type: string
                  default: ISO
                description: Viewport to query
            responses:
              200:
                description: Camera state
            """
            viewport = request.args.get("viewport", "ISO")
            try:
                state = self.viewer.get_camera_state(viewport)
                return jsonify(state)
            except KeyError:
                return jsonify({"error": f"Unknown viewport: {viewport}"}), 400

        # =====================================================================
        # Screenshot
        # =====================================================================

        @self.app.route("/screenshot", methods=["GET"])
        def get_screenshot():
            """
            Capture and return screenshot as PNG.

            ---
            tags:
              - Screenshot
            responses:
              200:
                description: PNG image
                content:
                  image/png:
                    schema:
                      type: string
                      format: binary
            """
            # Capture to temp file and return
            import tempfile
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
                filepath = self.viewer.screenshot(f.name)
                return send_file(
                    filepath,
                    mimetype="image/png",
                    as_attachment=False
                )

        @self.app.route("/screenshot", methods=["POST"])
        def save_screenshot():
            """
            Capture and save screenshot.

            ---
            tags:
              - Screenshot
            requestBody:
              content:
                application/json:
                  schema:
                    type: object
                    properties:
                      filename:
                        type: string
                        description: Output filename
                      directory:
                        type: string
                        description: Output directory (optional)
            responses:
              200:
                description: Screenshot saved
                content:
                  application/json:
                    schema:
                      type: object
                      properties:
                        filepath:
                          type: string
                          description: Full path to saved file
            """
            data = request.get_json() or {}
            filename = data.get("filename")
            directory = data.get("directory")

            filepath = self.viewer.screenshot(filename, directory)
            return jsonify({"filepath": filepath, "success": True})

        # =====================================================================
        # Part Movement
        # =====================================================================

        @self.app.route("/move", methods=["POST"])
        def move_part():
            """
            Move a part by relative offset.

            ---
            tags:
              - Parts
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    required:
                      - part
                      - dx
                      - dy
                      - dz
                    properties:
                      part:
                        type: string
                        description: Part name
                      dx:
                        type: number
                        description: X offset
                      dy:
                        type: number
                        description: Y offset
                      dz:
                        type: number
                        description: Z offset
            responses:
              200:
                description: Part moved
              404:
                description: Part not found
            """
            data = request.get_json() or {}

            try:
                self.viewer.move_part(
                    data["part"],
                    data["dx"],
                    data["dy"],
                    data["dz"]
                )
                return jsonify({"success": True})
            except KeyError as e:
                if "part" in str(e):
                    return jsonify({"error": f"Part not found: {data.get('part')}"}), 404
                return jsonify({"error": f"Missing required field: {e}"}), 400

        @self.app.route("/moveto", methods=["POST"])
        def move_part_to():
            """
            Move a part to absolute position.

            ---
            tags:
              - Parts
            requestBody:
              required: true
              content:
                application/json:
                  schema:
                    type: object
                    required:
                      - part
                      - x
                      - y
                      - z
                    properties:
                      part:
                        type: string
                        description: Part name
                      x:
                        type: number
                        description: Target X coordinate
                      y:
                        type: number
                        description: Target Y coordinate
                      z:
                        type: number
                        description: Target Z coordinate
            responses:
              200:
                description: Part moved
              404:
                description: Part not found
            """
            data = request.get_json() or {}

            try:
                self.viewer.move_part_to(
                    data["part"],
                    data["x"],
                    data["y"],
                    data["z"]
                )
                return jsonify({"success": True})
            except KeyError as e:
                if "part" in str(e):
                    return jsonify({"error": f"Part not found: {data.get('part')}"}), 404
                return jsonify({"error": f"Missing required field: {e}"}), 400

    def _register_socket_events(self):
        """Register WebSocket event handlers."""

        @self.socketio.on("connect")
        def handle_connect():
            """Handle client connection."""
            sid = request.sid
            self._subscriptions[sid] = set()
            # Send initial status
            status = self.viewer.get_status()
            emit(EventType.STATUS.value, status.to_dict())

        @self.socketio.on("disconnect")
        def handle_disconnect():
            """Handle client disconnection."""
            sid = request.sid
            if sid in self._subscriptions:
                del self._subscriptions[sid]

        @self.socketio.on(EventType.SUBSCRIBE.value)
        def handle_subscribe(data):
            """
            Handle subscription request.

            Data format::

                {
                    "events": ["part_loaded", "xray_changed"]
                }
            """
            sid = request.sid
            events = data.get("events", [])
            if sid not in self._subscriptions:
                self._subscriptions[sid] = set()
            self._subscriptions[sid].update(events)
            emit("subscribed", {"events": list(self._subscriptions[sid])})

        @self.socketio.on(EventType.UNSUBSCRIBE.value)
        def handle_unsubscribe(data):
            """
            Handle unsubscription request.

            Data format::

                {
                    "events": ["part_loaded"]
                }
            """
            sid = request.sid
            events = data.get("events", [])
            if sid in self._subscriptions:
                self._subscriptions[sid] -= set(events)
            emit("unsubscribed", {"events": events})

        @self.socketio.on(EventType.COMMAND.value)
        def handle_command(data):
            """
            Handle command from client.

            Data format::

                {
                    "command": "xray",
                    "params": {"enabled": true}
                }
            """
            try:
                cmd = CommandMessage.from_dict(data)
                self._execute_command(cmd)
                emit("command_result", {"success": True, "command": cmd.command})
            except Exception as e:
                emit("command_result", {
                    "success": False,
                    "command": data.get("command"),
                    "error": str(e)
                })

    def _execute_command(self, cmd: CommandMessage):
        """
        Execute a command message.

        Parameters
        ----------
        cmd : CommandMessage
            Command to execute.
        """
        params = cmd.params

        if cmd.command == CommandType.LOAD.value:
            if "positions_file" in params:
                self.viewer.load_from_json(params["positions_file"])
            elif "parts" in params:
                self.viewer.load_parts(params["parts"], clear=params.get("clear", True))
            else:
                self.viewer.load_from_json()

        elif cmd.command == CommandType.CLEAR.value:
            self.viewer.clear_parts()

        elif cmd.command == CommandType.RELOAD.value:
            self.viewer.reload_positions()

        elif cmd.command == CommandType.XRAY.value:
            if "enabled" in params:
                self.viewer.set_xray(params["enabled"], params.get("opacity"))
            else:
                self.viewer.toggle_xray()

        elif cmd.command == CommandType.SCREENSHOT.value:
            self.viewer.screenshot(
                params.get("filename"),
                params.get("directory")
            )

        elif cmd.command == CommandType.FOCUS.value:
            self.viewer.focus_part(params["part"])

        elif cmd.command == CommandType.HIGHLIGHT.value:
            self.viewer.highlight_parts(
                params.get("parts", []),
                params.get("dim_opacity", 0.15)
            )

        elif cmd.command == CommandType.CAMERA_RESET.value:
            self.viewer.reset_camera()

        elif cmd.command == CommandType.CAMERA_SET.value:
            self.viewer.set_camera(
                tuple(params["position"]),
                tuple(params["focal_point"]),
                params.get("viewport", "ISO"),
                tuple(params.get("view_up", [0, 0, 1]))
            )

        elif cmd.command == CommandType.CAMERA_ORBIT.value:
            self.viewer.orbit_camera(
                params["azimuth"],
                params["elevation"],
                params["distance"],
                params.get("viewport", "ISO")
            )

        elif cmd.command == CommandType.MOVE.value:
            self.viewer.move_part(
                params["part"],
                params["dx"],
                params["dy"],
                params["dz"]
            )

        elif cmd.command == CommandType.MOVE_TO.value:
            self.viewer.move_part_to(
                params["part"],
                params["x"],
                params["y"],
                params["z"]
            )

    def _connect_viewer_events(self):
        """Connect viewer events to WebSocket broadcasts."""

        def broadcast_event(event):
            """Broadcast event to subscribed clients."""
            event_type = event.event_type
            data = event.to_dict()

            # Broadcast to all clients (or filter by subscription)
            self.socketio.emit(event_type, data)

        # Register for all event types
        for event_type in [
            EventType.PART_LOADED.value,
            EventType.PART_REMOVED.value,
            EventType.PARTS_CLEARED.value,
            EventType.XRAY_CHANGED.value,
            EventType.CAMERA_CHANGED.value,
            EventType.SCREENSHOT_SAVED.value,
            EventType.SELECTION_CHANGED.value,
            EventType.ERROR.value,
        ]:
            self.viewer.on_event(event_type, broadcast_event)

    def run(
        self,
        host: str = "127.0.0.1",
        port: int = 5000,
        debug: bool = False,
        use_reloader: bool = False,
        viewer_in_background: bool = True
    ):
        """
        Start the API server.

        Parameters
        ----------
        host : str, optional
            Host to bind to. Default: "127.0.0.1".
        port : int, optional
            Port to listen on. Default: 5000.
        debug : bool, optional
            Enable Flask debug mode. Default: False.
        use_reloader : bool, optional
            Enable auto-reload on code changes. Default: False.
        viewer_in_background : bool, optional
            Run viewer rendering in background thread. Default: True.

        Notes
        -----
        This method blocks. The viewer render loop runs in a background
        thread if viewer_in_background is True.

        Examples
        --------
        Start server on default port::

            server.run()

        Start on custom host/port::

            server.run(host="0.0.0.0", port=8080)

        Start with debug mode::

            server.run(debug=True)
        """
        self._running = True

        # Start viewer in background thread if requested
        if viewer_in_background:
            self.viewer.start_offscreen()

        # Run Flask/SocketIO
        self.socketio.run(
            self.app,
            host=host,
            port=port,
            debug=debug,
            use_reloader=use_reloader
        )

    def stop(self):
        """Stop the API server."""
        self._running = False
        # SocketIO doesn't have a clean stop method in threading mode
        # The server will stop when the main thread exits
