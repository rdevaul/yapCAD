"""
Viewer Configuration
====================

Configuration dataclass for the VTK viewer and API server.

Example Usage
-------------

Basic configuration::

    from yapcad.viewer import ViewerConfig

    config = ViewerConfig(
        stl_dir="/path/to/stl/files",
        positions_file="/path/to/positions.json",
        window_size=(1600, 1000),
        background_color=(0.12, 0.12, 0.15)
    )

Configuration with file-based commands (backward compatibility)::

    config = ViewerConfig(
        stl_dir="/path/to/stl/files",
        command_file="/path/to/viewer_cmd.txt",
        screenshot_dir="/path/to/screenshots"
    )
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, Dict, Any


@dataclass
class ViewerConfig:
    """
    Configuration for the VTK assembly viewer.

    This dataclass holds all configuration options for the viewer including
    file paths, window settings, rendering options, and color schemes.

    Parameters
    ----------
    stl_dir : str or Path
        Directory containing STL files to load.
    positions_file : str or Path, optional
        JSON file containing part positions and transforms.
        Expected format::

            {
                "PART_NAME": {
                    "transform_matrix": [[...], [...], [...], [...]],
                    "stl_file": "optional/path/to/file.stl",
                    "color": [r, g, b]  # optional, 0-1 range
                },
                ...
            }

    command_file : str or Path, optional
        Path to file-based command interface for backward compatibility.
        Commands are read from this file and executed.
    screenshot_dir : str or Path, optional
        Directory for saving screenshots. Defaults to stl_dir/renders.
    window_size : tuple of int, optional
        Window dimensions (width, height). Default: (1600, 1000).
    window_title : str, optional
        Window title. Default: "yapCAD Assembly Viewer".
    background_color : tuple of float, optional
        Background color as RGB (0-1 range). Default: (0.12, 0.12, 0.15).
    default_opacity : float, optional
        Default opacity for parts. Default: 1.0.
    xray_opacity : float, optional
        Opacity when x-ray mode is enabled. Default: 0.4.
    default_color : tuple of float, optional
        Default part color as RGB (0-1 range). Default: (0.4, 0.4, 0.4).
    color_palette : dict, optional
        Named color palette for automatic color assignment.
        Keys are color names, values are RGB tuples (0-1 range).
    enable_multi_viewport : bool, optional
        Enable 4-viewport mode (ISO, TOP, FRONT, SIDE). Default: True.
    command_poll_interval_ms : int, optional
        Interval in milliseconds for polling command file. Default: 100.

    Attributes
    ----------
    stl_dir : Path
        Resolved STL directory path.
    positions_file : Path or None
        Resolved positions file path.
    command_file : Path or None
        Resolved command file path.
    screenshot_dir : Path
        Resolved screenshot directory path.

    Examples
    --------
    Minimal configuration::

        config = ViewerConfig(stl_dir="/my/stls")

    Full configuration::

        config = ViewerConfig(
            stl_dir="/my/stls",
            positions_file="/my/positions.json",
            window_size=(1920, 1080),
            background_color=(0.1, 0.1, 0.12),
            color_palette={
                "primary": (0.2, 0.4, 0.8),
                "secondary": (0.8, 0.4, 0.2),
                "highlight": (1.0, 0.9, 0.2)
            }
        )
    """

    stl_dir: str
    positions_file: Optional[str] = None
    command_file: Optional[str] = None
    screenshot_dir: Optional[str] = None

    # Window settings
    window_size: Tuple[int, int] = (1600, 1000)
    window_title: str = "yapCAD Assembly Viewer"
    background_color: Tuple[float, float, float] = (0.12, 0.12, 0.15)

    # Rendering settings
    default_opacity: float = 1.0
    xray_opacity: float = 0.4
    default_color: Tuple[float, float, float] = (0.4, 0.4, 0.4)

    # Color palette for automatic assignment
    color_palette: Dict[str, Tuple[float, float, float]] = field(
        default_factory=lambda: {
            "gray": (0.4, 0.4, 0.4),
            "red": (0.8, 0.2, 0.2),
            "green": (0.2, 0.8, 0.2),
            "blue": (0.2, 0.4, 0.8),
            "yellow": (0.8, 0.8, 0.2),
            "orange": (0.8, 0.5, 0.2),
            "purple": (0.6, 0.2, 0.8),
            "cyan": (0.2, 0.7, 0.8),
            "brown": (0.5, 0.3, 0.15),
            "dark": (0.15, 0.15, 0.15),
        }
    )

    # Viewport settings
    enable_multi_viewport: bool = True

    # Command file polling
    command_poll_interval_ms: int = 100

    def __post_init__(self):
        """Validate and convert paths."""
        self._stl_dir = Path(self.stl_dir)
        self._positions_file = Path(self.positions_file) if self.positions_file else None
        self._command_file = Path(self.command_file) if self.command_file else None

        if self.screenshot_dir:
            self._screenshot_dir = Path(self.screenshot_dir)
        else:
            self._screenshot_dir = self._stl_dir / "renders"

    @property
    def stl_path(self) -> Path:
        """Get the resolved STL directory path."""
        return self._stl_dir

    @property
    def positions_path(self) -> Optional[Path]:
        """Get the resolved positions file path."""
        return self._positions_file

    @property
    def command_path(self) -> Optional[Path]:
        """Get the resolved command file path."""
        return self._command_file

    @property
    def screenshot_path(self) -> Path:
        """Get the resolved screenshot directory path."""
        return self._screenshot_dir

    def get_color(self, name: str) -> Tuple[float, float, float]:
        """
        Get a color from the palette by name.

        Parameters
        ----------
        name : str
            Color name (case-insensitive).

        Returns
        -------
        tuple of float
            RGB color tuple (0-1 range).

        Examples
        --------
        >>> config = ViewerConfig(stl_dir="/tmp")
        >>> config.get_color("blue")
        (0.2, 0.4, 0.8)
        >>> config.get_color("unknown")  # Returns default
        (0.4, 0.4, 0.4)
        """
        return self.color_palette.get(name.lower(), self.default_color)

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.

        Returns
        -------
        dict
            Configuration as a dictionary suitable for JSON serialization.
        """
        return {
            "stl_dir": str(self.stl_dir),
            "positions_file": str(self.positions_file) if self.positions_file else None,
            "command_file": str(self.command_file) if self.command_file else None,
            "screenshot_dir": str(self._screenshot_dir),
            "window_size": list(self.window_size),
            "window_title": self.window_title,
            "background_color": list(self.background_color),
            "default_opacity": self.default_opacity,
            "xray_opacity": self.xray_opacity,
            "default_color": list(self.default_color),
            "color_palette": {k: list(v) for k, v in self.color_palette.items()},
            "enable_multi_viewport": self.enable_multi_viewport,
            "command_poll_interval_ms": self.command_poll_interval_ms,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ViewerConfig":
        """
        Create configuration from dictionary.

        Parameters
        ----------
        data : dict
            Configuration dictionary.

        Returns
        -------
        ViewerConfig
            New configuration instance.
        """
        # Convert lists back to tuples where needed
        if "window_size" in data:
            data["window_size"] = tuple(data["window_size"])
        if "background_color" in data:
            data["background_color"] = tuple(data["background_color"])
        if "default_color" in data:
            data["default_color"] = tuple(data["default_color"])
        if "color_palette" in data:
            data["color_palette"] = {
                k: tuple(v) for k, v in data["color_palette"].items()
            }
        return cls(**data)
