from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .routes.dsl import router as dsl_router
from .routes.geometry import router as geometry_router
from .routes.profile import router as profile_router
from .routes.render import router as render_router
from .routes.system import router as system_router
from .routes.ws import router as ws_router


def create_app() -> FastAPI:
    app = FastAPI(title="yapCAD Service", version="0.1.0")

    # CORS for yapCAD Web Viewer (Vite default)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    app.include_router(system_router)
    app.include_router(dsl_router)
    app.include_router(geometry_router)
    app.include_router(profile_router)
    app.include_router(render_router)
    app.include_router(ws_router)

    return app


app = create_app()
