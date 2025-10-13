from . import native as native

__all__ = ['native']

try:
    from . import trimesh_engine as trimesh
except Exception:  # optional dependency
    trimesh = None
else:
    __all__.append('trimesh')

ENGINE_REGISTRY = {'native': native}
if trimesh is not None:
    ENGINE_REGISTRY['trimesh'] = trimesh


def get_engine(name: str):
    return ENGINE_REGISTRY.get(name)


__all__.extend(['ENGINE_REGISTRY', 'get_engine'])
