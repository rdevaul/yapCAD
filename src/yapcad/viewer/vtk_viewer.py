"""
VTK-Based Assembly Viewer
=========================

Core VTK rendering engine for the yapCAD assembly viewer.

This module provides a general-purpose STL viewer with multi-viewport support,
transform application, and real-time rendering capabilities. It is designed
to work with any assembly defined by JSON position files and STL directories.

Example Usage
-------------

Basic usage::

    from yapcad.viewer import VTKViewer, ViewerConfig

    config = ViewerConfig(
        stl_dir="/path/to/stls",
        positions_file="/path/to/positions.json"
    )
    viewer = VTKViewer(config)
    viewer.load_from_json()
    viewer.start()

Manual part loading::

    viewer = VTKViewer(config)
    viewer.add_part(
        name="my_part",
        stl_path="/path/to/part.stl",
        transform=[[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]],
        color=(0.8, 0.4, 0.2)
    )
    viewer.start()

With file-based commands (backward compatibility)::

    config = ViewerConfig(
        stl_dir="/path/to/stls",
        command_file="/path/to/viewer_cmd.txt"
    )
    viewer = VTKViewer(config)
    viewer.start()  # Will poll command file for instructions

Dependencies
------------

Requires VTK >= 9.0::

    pip install vtk
"""

import json
import math
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import vtk

from .config import ViewerConfig
from .events import (
    PartLoadedEvent,
    PartRemovedEvent,
    PartsClearedEvent,
    XrayChangedEvent,
    CameraChangedEvent,
    ScreenshotSavedEvent,
    SelectionChangedEvent,
    StatusEvent,
    ErrorEvent,
)


# Type aliases
Matrix4x4 = List[List[float]]
RGB = Tuple[float, float, float]
EventCallback = Callable[[Any], None]


class VTKViewer:
    """
    VTK-based multi-viewport assembly viewer.

    This class provides a complete 3D viewer for STL assemblies with support
    for multiple viewports, transform application, x-ray mode, and various
    camera controls.

    Parameters
    ----------
    config : ViewerConfig
        Viewer configuration object.

    Attributes
    ----------
    config : ViewerConfig
        Current configuration.
    actors : dict
        Mapping of part names to VTK actors.
    positions : dict
        Current part positions loaded from JSON.
    xray_mode : bool
        Whether x-ray mode is currently enabled.

    Examples
    --------
    Create viewer and load from JSON::

        config = ViewerConfig(
            stl_dir="/my/stls",
            positions_file="/my/positions.json"
        )
        viewer = VTKViewer(config)
        viewer.load_from_json()
        viewer.start()

    Create viewer and add parts manually::

        viewer = VTKViewer(ViewerConfig(stl_dir="/my/stls"))
        viewer.add_part("gear", "/my/stls/gear.stl", transform, (0.8, 0.5, 0.2))
        viewer.add_part("housing", "/my/stls/housing.stl", transform, (0.4, 0.4, 0.4))
        viewer.start()
    """

    def __init__(self, config: ViewerConfig):
        """
        Initialize the VTK viewer.

        Parameters
        ----------
        config : ViewerConfig
            Viewer configuration.
        """
        self.config = config
        self.positions: Dict[str, Any] = {}
        self.actors: Dict[str, vtk.vtkActor] = {}
        self.actor_transforms: Dict[str, vtk.vtkTransform] = {}
        self.xray_mode = False
        self._color_index = 0
        self._event_callbacks: Dict[str, List[EventCallback]] = {}

        # Initialize VTK components
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.SetSize(*config.window_size)
        self.render_window.SetWindowName(config.window_title)

        self.renderers: Dict[str, vtk.vtkRenderer] = {}
        self._setup_viewports()

        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)

        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        # Set up keyboard handler
        self.interactor.AddObserver("KeyPressEvent", self._on_key)

        # Set up command file polling if configured
        self._last_cmd_mtime = 0
        if config.command_path:
            self.interactor.AddObserver("TimerEvent", self._on_timer)

        # Load positions if configured
        if config.positions_path and config.positions_path.exists():
            self._load_positions_file()

    def _setup_viewports(self):
        """Set up the multi-viewport layout."""
        if self.config.enable_multi_viewport:
            viewports = {
                "ISO": (0.5, 0.5, 1.0, 1.0),
                "TOP": (0.0, 0.5, 0.5, 1.0),
                "FRONT": (0.0, 0.0, 0.5, 0.5),
                "SIDE": (0.5, 0.0, 1.0, 0.5),
            }
        else:
            viewports = {
                "ISO": (0.0, 0.0, 1.0, 1.0),
            }

        bg = self.config.background_color

        for name, (x0, y0, x1, y1) in viewports.items():
            renderer = vtk.vtkRenderer()
            renderer.SetViewport(x0, y0, x1, y1)
            renderer.SetBackground(*bg)
            self.render_window.AddRenderer(renderer)
            self.renderers[name] = renderer

            # Add viewport label
            label = vtk.vtkTextActor()
            label.SetInput(name)
            label.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
            label.GetPositionCoordinate().SetValue(0.02, 0.92)
            label.GetTextProperty().SetFontSize(16)
            label.GetTextProperty().SetColor(0.8, 0.8, 0.8)
            label.GetTextProperty().SetBold(True)
            renderer.AddActor2D(label)

        self._reset_cameras()

    def _load_positions_file(self):
        """Load positions from the configured JSON file."""
        if self.config.positions_path and self.config.positions_path.exists():
            with open(self.config.positions_path) as f:
                self.positions = json.load(f)

    def _get_next_color(self) -> RGB:
        """
        Get the next color from the palette for automatic assignment.

        Returns
        -------
        tuple of float
            RGB color (0-1 range).
        """
        colors = list(self.config.color_palette.values())
        if not colors:
            return self.config.default_color
        color = colors[self._color_index % len(colors)]
        self._color_index += 1
        return color

    def _get_scene_bounds(self) -> Optional[Tuple[float, ...]]:
        """
        Calculate the bounding box of all loaded actors.

        Returns
        -------
        tuple of float or None
            Bounding box (xmin, xmax, ymin, ymax, zmin, zmax) or None if empty.
        """
        if not self.actors:
            return None
        all_bounds = [a.GetBounds() for a in self.actors.values()]
        return (
            min(b[0] for b in all_bounds),
            max(b[1] for b in all_bounds),
            min(b[2] for b in all_bounds),
            max(b[3] for b in all_bounds),
            min(b[4] for b in all_bounds),
            max(b[5] for b in all_bounds),
        )

    def _reset_cameras(self):
        """Reset all cameras to view the entire scene."""
        bounds = self._get_scene_bounds()
        if bounds is None:
            bounds = (-100, 100, -100, 100, 0, 200)

        cx = (bounds[0] + bounds[1]) / 2
        cy = (bounds[2] + bounds[3]) / 2
        cz = (bounds[4] + bounds[5]) / 2
        size = max(
            bounds[1] - bounds[0],
            bounds[3] - bounds[2],
            bounds[5] - bounds[4]
        )

        for name, renderer in self.renderers.items():
            camera = renderer.GetActiveCamera()
            camera.SetParallelProjection(name != "ISO")

            if name == "ISO":
                camera.SetPosition(cx + size, cy + size * 0.8, cz + size)
                camera.SetViewUp(0, 0, 1)
            elif name == "TOP":
                camera.SetPosition(cx, cy, cz + size * 1.5)
                camera.SetViewUp(0, 1, 0)
            elif name == "FRONT":
                camera.SetPosition(cx, cy - size * 1.5, cz)
                camera.SetViewUp(0, 0, 1)
            elif name == "SIDE":
                camera.SetPosition(cx + size * 1.5, cy, cz)
                camera.SetViewUp(0, 0, 1)

            camera.SetFocalPoint(cx, cy, cz)
            renderer.ResetCamera()

    def _load_stl(self, stl_path: Union[str, Path]) -> Optional[vtk.vtkPolyData]:
        """
        Load an STL file.

        Parameters
        ----------
        stl_path : str or Path
            Path to the STL file.

        Returns
        -------
        vtkPolyData or None
            Loaded mesh data or None if file doesn't exist.
        """
        path = Path(stl_path)
        if not path.exists():
            return None
        reader = vtk.vtkSTLReader()
        reader.SetFileName(str(path))
        reader.Update()
        return reader.GetOutput()

    def _emit_event(self, event):
        """
        Emit an event to all registered callbacks.

        Parameters
        ----------
        event : BaseEvent
            Event to emit.
        """
        callbacks = self._event_callbacks.get(event.event_type, [])
        for callback in callbacks:
            try:
                callback(event)
            except Exception as e:
                print(f"Event callback error: {e}")

    # =========================================================================
    # Public API - Part Management
    # =========================================================================

    def add_part(
        self,
        name: str,
        stl_path: Union[str, Path],
        transform: Optional[Matrix4x4] = None,
        color: Optional[RGB] = None
    ) -> bool:
        """
        Add a single part to the viewer.

        Parameters
        ----------
        name : str
            Unique identifier for the part.
        stl_path : str or Path
            Path to the STL file.
        transform : list of list of float, optional
            4x4 transformation matrix. Identity matrix if not provided.
        color : tuple of float, optional
            RGB color (0-1 range). Auto-assigned if not provided.

        Returns
        -------
        bool
            True if part was loaded successfully.

        Examples
        --------
        Add a part with transform::

            transform = [
                [1, 0, 0, 10],
                [0, 1, 0, 20],
                [0, 0, 1, 30],
                [0, 0, 0, 1]
            ]
            viewer.add_part("my_gear", "/path/to/gear.stl", transform, (0.8, 0.5, 0.2))

        Add a part at origin with auto color::

            viewer.add_part("housing", "/path/to/housing.stl")
        """
        polydata = self._load_stl(stl_path)
        if polydata is None:
            self._emit_event(ErrorEvent(
                message=f"Failed to load STL: {stl_path}",
                code="LOAD_ERROR",
                details={"part": name, "path": str(stl_path)}
            ))
            return False

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Apply color
        if color is None:
            color = self._get_next_color()
        actor.GetProperty().SetColor(*color)
        actor.GetProperty().SetInterpolationToPhong()

        # Apply transform
        if transform is not None:
            vtk_transform = vtk.vtkTransform()
            vtk_matrix = vtk.vtkMatrix4x4()
            for i in range(4):
                for j in range(4):
                    vtk_matrix.SetElement(i, j, transform[i][j])
            vtk_transform.SetMatrix(vtk_matrix)
            actor.SetUserTransform(vtk_transform)
            self.actor_transforms[name] = vtk_transform

        # Apply current x-ray mode
        if self.xray_mode:
            actor.GetProperty().SetOpacity(self.config.xray_opacity)

        # Add to all renderers
        for renderer in self.renderers.values():
            renderer.AddActor(actor)

        self.actors[name] = actor

        # Emit event
        bounds = list(actor.GetBounds())
        self._emit_event(PartLoadedEvent(
            part_name=name,
            stl_file=str(stl_path),
            bounds=bounds,
            color=color
        ))

        return True

    def remove_part(self, name: str) -> bool:
        """
        Remove a part from the viewer.

        Parameters
        ----------
        name : str
            Name of the part to remove.

        Returns
        -------
        bool
            True if part was found and removed.
        """
        if name not in self.actors:
            return False

        actor = self.actors.pop(name)
        for renderer in self.renderers.values():
            renderer.RemoveActor(actor)

        if name in self.actor_transforms:
            del self.actor_transforms[name]

        self._emit_event(PartRemovedEvent(part_name=name))
        return True

    def clear_parts(self):
        """Remove all parts from the viewer."""
        count = len(self.actors)

        for actor in self.actors.values():
            for renderer in self.renderers.values():
                renderer.RemoveActor(actor)

        self.actors.clear()
        self.actor_transforms.clear()
        self._color_index = 0

        self._emit_event(PartsClearedEvent(previous_count=count))

    def load_parts(
        self,
        part_specs: List[Tuple[str, str, Optional[str]]],
        clear: bool = True
    ) -> int:
        """
        Load multiple parts from specifications.

        Parameters
        ----------
        part_specs : list of tuple
            List of (stl_file, part_name, color_name) tuples.
            stl_file is relative to stl_dir.
            color_name is optional; if None, auto-assigned.
        clear : bool, optional
            Clear existing parts first. Default: True.

        Returns
        -------
        int
            Number of parts successfully loaded.

        Examples
        --------
        Load multiple parts::

            parts = [
                ("gearbox/sun_gear.stl", "SUN_GEAR", "orange"),
                ("gearbox/ring_housing.stl", "RING_HOUSING", "brown"),
                ("motors/servo.stl", "SERVO", "dark"),
            ]
            viewer.load_parts(parts)
        """
        if clear:
            self.clear_parts()

        loaded = 0
        for spec in part_specs:
            stl_file = spec[0]
            part_name = spec[1]
            color_name = spec[2] if len(spec) > 2 else None

            stl_path = self.config.stl_path / stl_file
            color = self.config.get_color(color_name) if color_name else None

            # Get transform from positions if available
            transform = None
            if part_name in self.positions:
                transform = self.positions[part_name].get("transform_matrix")
                # Override color from positions if specified
                if "color" in self.positions[part_name] and color is None:
                    color = tuple(self.positions[part_name]["color"])

            if self.add_part(part_name, stl_path, transform, color):
                loaded += 1

        self._reset_cameras()
        return loaded

    def load_from_json(self, positions_file: Optional[Union[str, Path]] = None) -> int:
        """
        Load all parts defined in a positions JSON file.

        The JSON file should have the format::

            {
                "PART_NAME": {
                    "transform_matrix": [[...], [...], [...], [...]],
                    "stl_file": "path/to/file.stl",  // optional
                    "color": [r, g, b]  // optional, 0-1 range
                },
                ...
            }

        Parameters
        ----------
        positions_file : str or Path, optional
            Path to positions JSON file. Uses config.positions_file if not provided.

        Returns
        -------
        int
            Number of parts successfully loaded.

        Raises
        ------
        FileNotFoundError
            If positions file doesn't exist.
        """
        if positions_file:
            path = Path(positions_file)
        elif self.config.positions_path:
            path = self.config.positions_path
        else:
            raise ValueError("No positions file specified")

        if not path.exists():
            raise FileNotFoundError(f"Positions file not found: {path}")

        with open(path) as f:
            self.positions = json.load(f)

        self.clear_parts()
        loaded = 0

        for part_name, part_data in self.positions.items():
            # Get STL path
            if "stl_file" in part_data:
                stl_path = self.config.stl_path / part_data["stl_file"]
            else:
                # Try to find STL by part name
                possible_paths = [
                    self.config.stl_path / f"{part_name}.stl",
                    self.config.stl_path / f"{part_name.lower()}.stl",
                ]
                stl_path = None
                for p in possible_paths:
                    if p.exists():
                        stl_path = p
                        break

                if stl_path is None:
                    continue

            transform = part_data.get("transform_matrix")
            color = tuple(part_data["color"]) if "color" in part_data else None

            if self.add_part(part_name, stl_path, transform, color):
                loaded += 1

        self._reset_cameras()
        return loaded

    def reload_positions(self):
        """
        Reload positions from the JSON file and update transforms.

        This updates transforms without reloading STL files, which is faster
        for iterating on kinematic positions.
        """
        self._load_positions_file()

        for part_name, part_data in self.positions.items():
            if part_name not in self.actor_transforms:
                continue

            transform = part_data.get("transform_matrix")
            if transform is None:
                continue

            vtk_transform = self.actor_transforms[part_name]
            vtk_matrix = vtk.vtkMatrix4x4()
            for i in range(4):
                for j in range(4):
                    vtk_matrix.SetElement(i, j, transform[i][j])
            vtk_transform.SetMatrix(vtk_matrix)

        self.render_window.Render()

    def get_part_names(self) -> List[str]:
        """
        Get list of loaded part names.

        Returns
        -------
        list of str
            Names of all loaded parts.
        """
        return list(self.actors.keys())

    def get_part_bounds(self, name: str) -> Optional[Tuple[float, ...]]:
        """
        Get bounding box of a specific part.

        Parameters
        ----------
        name : str
            Part name.

        Returns
        -------
        tuple of float or None
            Bounding box (xmin, xmax, ymin, ymax, zmin, zmax) or None.
        """
        if name not in self.actors:
            return None
        return self.actors[name].GetBounds()

    # =========================================================================
    # Public API - Rendering Controls
    # =========================================================================

    def set_xray(self, enabled: bool, opacity: Optional[float] = None):
        """
        Enable or disable x-ray mode.

        Parameters
        ----------
        enabled : bool
            Whether to enable x-ray mode.
        opacity : float, optional
            Opacity value when x-ray is enabled. Uses config default if not provided.

        Examples
        --------
        Enable x-ray::

            viewer.set_xray(True)

        Enable with custom opacity::

            viewer.set_xray(True, opacity=0.3)

        Disable::

            viewer.set_xray(False)
        """
        self.xray_mode = enabled
        if opacity is None:
            opacity = self.config.xray_opacity if enabled else self.config.default_opacity
        else:
            opacity = opacity if enabled else self.config.default_opacity

        for actor in self.actors.values():
            actor.GetProperty().SetOpacity(opacity)

        self.render_window.Render()
        self._emit_event(XrayChangedEvent(enabled=enabled, opacity=opacity))

    def toggle_xray(self):
        """Toggle x-ray mode on/off."""
        self.set_xray(not self.xray_mode)

    def highlight_parts(self, names: List[str], dim_opacity: float = 0.15):
        """
        Highlight specific parts by dimming others.

        Parameters
        ----------
        names : list of str
            Part names to highlight.
        dim_opacity : float, optional
            Opacity for non-highlighted parts. Default: 0.15.

        Examples
        --------
        Highlight gears::

            viewer.highlight_parts(["SUN_GEAR", "RING_HOUSING"])

        Clear highlighting (show all)::

            viewer.highlight_parts([])  # Or pass all part names
        """
        for name, actor in self.actors.items():
            if not names or name in names:
                actor.GetProperty().SetOpacity(1.0)
                actor.GetProperty().SetAmbient(0.3)
            else:
                actor.GetProperty().SetOpacity(dim_opacity)
                actor.GetProperty().SetAmbient(0.0)

        self.render_window.Render()
        self._emit_event(SelectionChangedEvent(highlighted_parts=names))

    def focus_part(self, name: str):
        """
        Focus all cameras on a specific part.

        Parameters
        ----------
        name : str
            Name of the part to focus on.

        Raises
        ------
        KeyError
            If part name is not found.
        """
        if name not in self.actors:
            raise KeyError(f"Part not found: {name}")

        actor = self.actors[name]
        bounds = actor.GetBounds()
        cx = (bounds[0] + bounds[1]) / 2
        cy = (bounds[2] + bounds[3]) / 2
        cz = (bounds[4] + bounds[5]) / 2
        size = max(
            bounds[1] - bounds[0],
            bounds[3] - bounds[2],
            bounds[5] - bounds[4]
        )

        for vp_name, renderer in self.renderers.items():
            camera = renderer.GetActiveCamera()
            camera.SetFocalPoint(cx, cy, cz)

            if vp_name == "ISO":
                camera.SetPosition(cx + size * 1.5, cy + size * 1.5, cz + size * 1.5)
            elif vp_name == "TOP":
                camera.SetPosition(cx, cy, cz + size * 2)
            elif vp_name == "FRONT":
                camera.SetPosition(cx, cy - size * 2, cz)
            elif vp_name == "SIDE":
                camera.SetPosition(cx + size * 2, cy, cz)

            renderer.ResetCameraClippingRange()

        self.render_window.Render()

    # =========================================================================
    # Public API - Camera Controls
    # =========================================================================

    def reset_camera(self):
        """Reset all cameras to view the entire scene."""
        self._reset_cameras()
        self.render_window.Render()

    def set_camera(
        self,
        position: Tuple[float, float, float],
        focal_point: Tuple[float, float, float],
        viewport: str = "ISO",
        view_up: Tuple[float, float, float] = (0, 0, 1)
    ):
        """
        Set camera position and orientation.

        Parameters
        ----------
        position : tuple of float
            Camera position (x, y, z).
        focal_point : tuple of float
            Point camera is looking at (x, y, z).
        viewport : str, optional
            Viewport to modify. Default: "ISO".
        view_up : tuple of float, optional
            Up direction vector. Default: (0, 0, 1).

        Examples
        --------
        Set camera for isometric view::

            viewer.set_camera(
                position=(200, 200, 150),
                focal_point=(0, 0, 50)
            )
        """
        if viewport not in self.renderers:
            raise KeyError(f"Unknown viewport: {viewport}")

        camera = self.renderers[viewport].GetActiveCamera()
        camera.SetPosition(*position)
        camera.SetFocalPoint(*focal_point)
        camera.SetViewUp(*view_up)
        self.renderers[viewport].ResetCameraClippingRange()
        self.render_window.Render()

        self._emit_event(CameraChangedEvent(
            viewport=viewport,
            position=position,
            focal_point=focal_point,
            view_up=view_up
        ))

    def orbit_camera(
        self,
        azimuth: float,
        elevation: float,
        distance: float,
        viewport: str = "ISO"
    ):
        """
        Position camera using spherical coordinates around focal point.

        Parameters
        ----------
        azimuth : float
            Azimuth angle in degrees (rotation around Z axis).
        elevation : float
            Elevation angle in degrees (angle from XY plane).
        distance : float
            Distance from focal point.
        viewport : str, optional
            Viewport to modify. Default: "ISO".

        Examples
        --------
        Orbit to 45 degrees azimuth, 30 degrees elevation::

            viewer.orbit_camera(azimuth=45, elevation=30, distance=500)
        """
        if viewport not in self.renderers:
            raise KeyError(f"Unknown viewport: {viewport}")

        camera = self.renderers[viewport].GetActiveCamera()
        fx, fy, fz = camera.GetFocalPoint()

        az_rad = math.radians(azimuth)
        el_rad = math.radians(elevation)

        x = fx + distance * math.cos(el_rad) * math.cos(az_rad)
        y = fy + distance * math.cos(el_rad) * math.sin(az_rad)
        z = fz + distance * math.sin(el_rad)

        camera.SetPosition(x, y, z)
        camera.SetViewUp(0, 0, 1)
        self.renderers[viewport].ResetCameraClippingRange()
        self.render_window.Render()

        self._emit_event(CameraChangedEvent(
            viewport=viewport,
            position=(x, y, z),
            focal_point=(fx, fy, fz),
            view_up=(0, 0, 1)
        ))

    def get_camera_state(self, viewport: str = "ISO") -> Dict[str, Any]:
        """
        Get current camera state.

        Parameters
        ----------
        viewport : str, optional
            Viewport to query. Default: "ISO".

        Returns
        -------
        dict
            Camera state with position, focal_point, view_up.
        """
        if viewport not in self.renderers:
            raise KeyError(f"Unknown viewport: {viewport}")

        camera = self.renderers[viewport].GetActiveCamera()
        return {
            "position": camera.GetPosition(),
            "focal_point": camera.GetFocalPoint(),
            "view_up": camera.GetViewUp(),
        }

    # =========================================================================
    # Public API - Screenshot
    # =========================================================================

    def screenshot(
        self,
        filename: Optional[str] = None,
        directory: Optional[Union[str, Path]] = None
    ) -> str:
        """
        Capture a screenshot.

        Parameters
        ----------
        filename : str, optional
            Output filename. Default: "screenshot.png".
        directory : str or Path, optional
            Output directory. Uses config.screenshot_dir if not provided.

        Returns
        -------
        str
            Full path to saved screenshot.

        Examples
        --------
        Save with default name::

            path = viewer.screenshot()

        Save with custom name::

            path = viewer.screenshot("my_view.png")

        Save to custom directory::

            path = viewer.screenshot("view.png", "/my/images")
        """
        if filename is None:
            filename = "screenshot.png"

        if directory is None:
            directory = self.config.screenshot_path
        else:
            directory = Path(directory)

        directory.mkdir(parents=True, exist_ok=True)
        filepath = directory / filename

        self.render_window.Render()

        w2i = vtk.vtkWindowToImageFilter()
        w2i.SetInput(self.render_window)
        w2i.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(str(filepath))
        writer.SetInputConnection(w2i.GetOutputPort())
        writer.Write()

        size = self.render_window.GetSize()
        self._emit_event(ScreenshotSavedEvent(
            filename=filename,
            filepath=str(filepath),
            size=tuple(size)
        ))

        return str(filepath)

    # =========================================================================
    # Public API - Part Movement
    # =========================================================================

    def move_part(self, name: str, dx: float, dy: float, dz: float):
        """
        Move a part by relative offset.

        Parameters
        ----------
        name : str
            Part name.
        dx : float
            X offset.
        dy : float
            Y offset.
        dz : float
            Z offset.

        Raises
        ------
        KeyError
            If part name is not found.
        """
        if name not in self.actors:
            raise KeyError(f"Part not found: {name}")

        if name not in self.actor_transforms:
            self.actor_transforms[name] = vtk.vtkTransform()
            self.actors[name].SetUserTransform(self.actor_transforms[name])

        transform = self.actor_transforms[name]
        transform.Translate(dx, dy, dz)
        self.render_window.Render()

    def move_part_to(self, name: str, x: float, y: float, z: float):
        """
        Move a part to absolute position (centers part at position).

        Parameters
        ----------
        name : str
            Part name.
        x : float
            Target X coordinate.
        y : float
            Target Y coordinate.
        z : float
            Target Z coordinate.

        Raises
        ------
        KeyError
            If part name is not found.
        """
        if name not in self.actors:
            raise KeyError(f"Part not found: {name}")

        actor = self.actors[name]
        bounds = actor.GetBounds()
        cx = (bounds[0] + bounds[1]) / 2
        cy = (bounds[2] + bounds[3]) / 2
        cz = (bounds[4] + bounds[5]) / 2

        dx, dy, dz = x - cx, y - cy, z - cz
        self.move_part(name, dx, dy, dz)

    # =========================================================================
    # Public API - Event Handling
    # =========================================================================

    def on_event(self, event_type: str, callback: EventCallback):
        """
        Register a callback for viewer events.

        Parameters
        ----------
        event_type : str
            Event type to listen for (e.g., "part_loaded", "xray_changed").
        callback : callable
            Function to call when event occurs. Receives event object.

        Examples
        --------
        Track loaded parts::

            def on_load(event):
                print(f"Loaded: {event.part_name}")

            viewer.on_event("part_loaded", on_load)
        """
        if event_type not in self._event_callbacks:
            self._event_callbacks[event_type] = []
        self._event_callbacks[event_type].append(callback)

    def get_status(self) -> StatusEvent:
        """
        Get current viewer status.

        Returns
        -------
        StatusEvent
            Current viewer status.
        """
        camera_state = self.get_camera_state("ISO")
        return StatusEvent(
            part_count=len(self.actors),
            xray_enabled=self.xray_mode,
            viewport="ISO",
            camera_position=camera_state["position"],
            camera_focal_point=camera_state["focal_point"]
        )

    # =========================================================================
    # Public API - Main Loop
    # =========================================================================

    def render(self):
        """Trigger a render update."""
        self.render_window.Render()

    def start(self):
        """
        Start the viewer main loop (blocking).

        This starts the VTK interactor and blocks until the window is closed.
        Use this for standalone viewer applications.

        Examples
        --------
        Basic standalone usage::

            viewer = VTKViewer(config)
            viewer.load_from_json()
            viewer.start()  # Blocks until window closed
        """
        self.render_window.Render()
        self.interactor.Initialize()

        # Start command file timer if configured
        if self.config.command_path:
            self.interactor.CreateRepeatingTimer(self.config.command_poll_interval_ms)

        self.interactor.Start()

    def start_offscreen(self):
        """
        Initialize for offscreen rendering.

        Use this when running headless or with the API server.
        Does not block.
        """
        self.render_window.SetOffScreenRendering(1)
        self.render_window.Render()

    # =========================================================================
    # Internal - Keyboard Handling
    # =========================================================================

    def _on_key(self, obj, event):
        """Handle keyboard events."""
        key = self.interactor.GetKeySym().lower()

        if key == "x":
            self.toggle_xray()
        elif key == "r":
            self.reset_camera()
        elif key == "q":
            self.interactor.TerminateApp()

    # =========================================================================
    # Internal - Command File Handling (Backward Compatibility)
    # =========================================================================

    def _on_timer(self, obj, event):
        """Timer callback for command file polling."""
        if not self.config.command_path:
            return

        try:
            cmd_file = self.config.command_path
            if cmd_file.exists():
                mtime = cmd_file.stat().st_mtime
                if mtime > self._last_cmd_mtime:
                    self._last_cmd_mtime = mtime
                    cmd = cmd_file.read_text().strip()
                    if cmd:
                        self._process_file_command(cmd)
                        cmd_file.write_text("")
        except Exception as e:
            print(f"Command file error: {e}")

    def _process_file_command(self, cmd: str):
        """
        Process a command from the file-based interface.

        Parameters
        ----------
        cmd : str
            Command string to process.
        """
        parts = cmd.lower().split()
        if not parts:
            return

        action = parts[0]
        print(f"[CMD] {cmd}")

        try:
            if action == "xray":
                self.toggle_xray()
            elif action == "reload":
                self.reload_positions()
            elif action == "reset":
                self.reset_camera()
            elif action == "screenshot":
                filename = parts[1] if len(parts) > 1 else None
                self.screenshot(filename)
            elif action == "focus":
                if len(parts) > 1:
                    self.focus_part(parts[1].upper())
            elif action == "highlight":
                names = [p.upper() for p in parts[1:]]
                self.highlight_parts(names)
            elif action == "camera":
                if len(parts) >= 4:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    fx, fy, fz = (
                        (float(parts[4]), float(parts[5]), float(parts[6]))
                        if len(parts) >= 7 else (0, 0, 100)
                    )
                    self.set_camera((x, y, z), (fx, fy, fz))
            elif action == "orbit":
                if len(parts) >= 4:
                    azimuth = float(parts[1])
                    elevation = float(parts[2])
                    distance = float(parts[3])
                    self.orbit_camera(azimuth, elevation, distance)
            elif action == "move":
                if len(parts) >= 5:
                    name = parts[1].upper()
                    dx, dy, dz = float(parts[2]), float(parts[3]), float(parts[4])
                    self.move_part(name, dx, dy, dz)
            elif action == "moveto":
                if len(parts) >= 5:
                    name = parts[1].upper()
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    self.move_part_to(name, x, y, z)
        except Exception as e:
            print(f"Command error: {e}")
            self._emit_event(ErrorEvent(
                message=str(e),
                code="COMMAND_ERROR",
                details={"command": cmd}
            ))
