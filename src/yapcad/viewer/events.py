"""
WebSocket Event Types
=====================

Defines event types for WebSocket communication between the viewer server
and connected clients.

Event Categories
----------------

**Client-to-Server Events:**
    Commands sent from clients to control the viewer.

**Server-to-Client Events:**
    Notifications sent from the viewer to connected clients.

Example Usage
-------------

Sending events from server::

    from yapcad.viewer.events import ServerEvent, PartLoadedEvent

    event = PartLoadedEvent(
        part_name="AXIS1_SUN_GEAR",
        stl_file="gearbox/axis1_sun_gear.stl",
        bounds=[-10, 10, -10, 10, 0, 20]
    )
    socketio.emit(event.event_type, event.to_dict())

Handling events on client::

    socket.on('part_loaded', (data) => {
        console.log(`Loaded: ${data.part_name}`);
    });
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


class EventType(str, Enum):
    """
    Enumeration of all WebSocket event types.

    Client-to-Server
    ----------------
    COMMAND
        Execute a viewer command.
    SUBSCRIBE
        Subscribe to specific event types.
    UNSUBSCRIBE
        Unsubscribe from event types.

    Server-to-Client
    ----------------
    STATUS
        Viewer status update.
    PART_LOADED
        Part was loaded into the scene.
    PART_REMOVED
        Part was removed from the scene.
    PARTS_CLEARED
        All parts were cleared.
    XRAY_CHANGED
        X-ray mode was toggled.
    CAMERA_CHANGED
        Camera position/orientation changed.
    SCREENSHOT_SAVED
        Screenshot was captured and saved.
    SELECTION_CHANGED
        Part selection changed.
    ERROR
        An error occurred.
    """
    # Client-to-Server
    COMMAND = "command"
    SUBSCRIBE = "subscribe"
    UNSUBSCRIBE = "unsubscribe"

    # Server-to-Client
    STATUS = "status"
    PART_LOADED = "part_loaded"
    PART_REMOVED = "part_removed"
    PARTS_CLEARED = "parts_cleared"
    XRAY_CHANGED = "xray_changed"
    CAMERA_CHANGED = "camera_changed"
    SCREENSHOT_SAVED = "screenshot_saved"
    SELECTION_CHANGED = "selection_changed"
    ERROR = "error"


@dataclass
class BaseEvent:
    """
    Base class for all WebSocket events.

    Attributes
    ----------
    event_type : str
        The event type identifier.
    """

    event_type: str = field(init=False)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert event to dictionary for JSON serialization.

        Returns
        -------
        dict
            Event data as dictionary.
        """
        return {
            k: v for k, v in self.__dict__.items()
            if not k.startswith("_") and k != "event_type"
        }


# =============================================================================
# Server-to-Client Events
# =============================================================================


@dataclass
class StatusEvent(BaseEvent):
    """
    Viewer status update event.

    Sent periodically or on request to provide current viewer state.

    Attributes
    ----------
    part_count : int
        Number of parts currently loaded.
    xray_enabled : bool
        Whether x-ray mode is active.
    viewport : str
        Current active viewport name.
    camera_position : tuple of float
        Camera position (x, y, z).
    camera_focal_point : tuple of float
        Camera focal point (x, y, z).

    WebSocket Event
    ---------------
    Event name: ``status``

    Payload::

        {
            "part_count": 42,
            "xray_enabled": false,
            "viewport": "ISO",
            "camera_position": [100.0, 100.0, 100.0],
            "camera_focal_point": [0.0, 0.0, 50.0]
        }
    """

    part_count: int
    xray_enabled: bool
    viewport: str = "ISO"
    camera_position: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    camera_focal_point: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self):
        self.event_type = EventType.STATUS.value


@dataclass
class PartLoadedEvent(BaseEvent):
    """
    Event emitted when a part is loaded.

    Attributes
    ----------
    part_name : str
        Name/identifier of the loaded part.
    stl_file : str
        Path to the STL file.
    bounds : list of float
        Bounding box [xmin, xmax, ymin, ymax, zmin, zmax].
    color : tuple of float, optional
        RGB color (0-1 range).

    WebSocket Event
    ---------------
    Event name: ``part_loaded``

    Payload::

        {
            "part_name": "AXIS1_SUN_GEAR",
            "stl_file": "gearbox/axis1_sun_gear.stl",
            "bounds": [-10.0, 10.0, -10.0, 10.0, 0.0, 20.0],
            "color": [0.8, 0.5, 0.2]
        }
    """

    part_name: str
    stl_file: str
    bounds: List[float] = field(default_factory=list)
    color: Optional[Tuple[float, float, float]] = None

    def __post_init__(self):
        self.event_type = EventType.PART_LOADED.value


@dataclass
class PartRemovedEvent(BaseEvent):
    """
    Event emitted when a part is removed.

    Attributes
    ----------
    part_name : str
        Name of the removed part.

    WebSocket Event
    ---------------
    Event name: ``part_removed``

    Payload::

        {
            "part_name": "AXIS1_SUN_GEAR"
        }
    """

    part_name: str

    def __post_init__(self):
        self.event_type = EventType.PART_REMOVED.value


@dataclass
class PartsClearedEvent(BaseEvent):
    """
    Event emitted when all parts are cleared.

    Attributes
    ----------
    previous_count : int
        Number of parts that were cleared.

    WebSocket Event
    ---------------
    Event name: ``parts_cleared``

    Payload::

        {
            "previous_count": 42
        }
    """

    previous_count: int

    def __post_init__(self):
        self.event_type = EventType.PARTS_CLEARED.value


@dataclass
class XrayChangedEvent(BaseEvent):
    """
    Event emitted when x-ray mode changes.

    Attributes
    ----------
    enabled : bool
        Whether x-ray mode is now enabled.
    opacity : float
        Current opacity value.

    WebSocket Event
    ---------------
    Event name: ``xray_changed``

    Payload::

        {
            "enabled": true,
            "opacity": 0.4
        }
    """

    enabled: bool
    opacity: float

    def __post_init__(self):
        self.event_type = EventType.XRAY_CHANGED.value


@dataclass
class CameraChangedEvent(BaseEvent):
    """
    Event emitted when camera position changes.

    Attributes
    ----------
    viewport : str
        Viewport name that changed.
    position : tuple of float
        New camera position (x, y, z).
    focal_point : tuple of float
        New focal point (x, y, z).
    view_up : tuple of float
        View up vector.

    WebSocket Event
    ---------------
    Event name: ``camera_changed``

    Payload::

        {
            "viewport": "ISO",
            "position": [100.0, 100.0, 100.0],
            "focal_point": [0.0, 0.0, 50.0],
            "view_up": [0.0, 0.0, 1.0]
        }
    """

    viewport: str
    position: Tuple[float, float, float]
    focal_point: Tuple[float, float, float]
    view_up: Tuple[float, float, float] = (0.0, 0.0, 1.0)

    def __post_init__(self):
        self.event_type = EventType.CAMERA_CHANGED.value


@dataclass
class ScreenshotSavedEvent(BaseEvent):
    """
    Event emitted when a screenshot is saved.

    Attributes
    ----------
    filename : str
        Name of the saved file.
    filepath : str
        Full path to the saved file.
    size : tuple of int
        Image dimensions (width, height).

    WebSocket Event
    ---------------
    Event name: ``screenshot_saved``

    Payload::

        {
            "filename": "screenshot_001.png",
            "filepath": "/path/to/renders/screenshot_001.png",
            "size": [1600, 1000]
        }
    """

    filename: str
    filepath: str
    size: Tuple[int, int] = (0, 0)

    def __post_init__(self):
        self.event_type = EventType.SCREENSHOT_SAVED.value


@dataclass
class SelectionChangedEvent(BaseEvent):
    """
    Event emitted when part selection changes.

    Attributes
    ----------
    selected_parts : list of str
        Names of currently selected parts.
    highlighted_parts : list of str
        Names of currently highlighted parts.

    WebSocket Event
    ---------------
    Event name: ``selection_changed``

    Payload::

        {
            "selected_parts": ["AXIS1_SUN_GEAR"],
            "highlighted_parts": ["AXIS1_SUN_GEAR", "AXIS1_RING_HOUSING"]
        }
    """

    selected_parts: List[str] = field(default_factory=list)
    highlighted_parts: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.event_type = EventType.SELECTION_CHANGED.value


@dataclass
class ErrorEvent(BaseEvent):
    """
    Event emitted when an error occurs.

    Attributes
    ----------
    message : str
        Error message.
    code : str, optional
        Error code for programmatic handling.
    details : dict, optional
        Additional error details.

    WebSocket Event
    ---------------
    Event name: ``error``

    Payload::

        {
            "message": "Failed to load STL file",
            "code": "LOAD_ERROR",
            "details": {"file": "missing.stl"}
        }
    """

    message: str
    code: Optional[str] = None
    details: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        self.event_type = EventType.ERROR.value


# =============================================================================
# Client-to-Server Command Types
# =============================================================================


class CommandType(str, Enum):
    """
    Types of commands that can be sent to the viewer.

    Commands
    --------
    LOAD
        Load parts into the viewer.
    CLEAR
        Clear all loaded parts.
    RELOAD
        Reload positions from JSON file.
    XRAY
        Toggle or set x-ray mode.
    SCREENSHOT
        Capture a screenshot.
    FOCUS
        Focus camera on a specific part.
    HIGHLIGHT
        Highlight specific parts.
    CAMERA_RESET
        Reset camera to default position.
    CAMERA_SET
        Set camera position directly.
    CAMERA_ORBIT
        Orbit camera around focal point.
    MOVE
        Move a part by relative offset.
    MOVE_TO
        Move a part to absolute position.
    """
    LOAD = "load"
    CLEAR = "clear"
    RELOAD = "reload"
    XRAY = "xray"
    SCREENSHOT = "screenshot"
    FOCUS = "focus"
    HIGHLIGHT = "highlight"
    CAMERA_RESET = "camera_reset"
    CAMERA_SET = "camera_set"
    CAMERA_ORBIT = "camera_orbit"
    MOVE = "move"
    MOVE_TO = "move_to"


@dataclass
class CommandMessage:
    """
    Command message sent from client to server.

    Attributes
    ----------
    command : str
        Command type (from CommandType enum).
    params : dict, optional
        Command parameters.

    Example
    -------
    Load command::

        {
            "command": "load",
            "params": {
                "parts": ["AXIS1_SUN_GEAR", "AXIS1_RING_HOUSING"],
                "clear": true
            }
        }

    Screenshot command::

        {
            "command": "screenshot",
            "params": {
                "filename": "my_screenshot.png"
            }
        }
    """

    command: str
    params: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CommandMessage":
        """
        Create command from dictionary.

        Parameters
        ----------
        data : dict
            Command data with 'command' and optional 'params'.

        Returns
        -------
        CommandMessage
            Parsed command message.
        """
        return cls(
            command=data.get("command", ""),
            params=data.get("params", {})
        )
