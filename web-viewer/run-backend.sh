#!/bin/bash
# Run the yapCAD Package Server backend
cd "$(dirname "$0")"
uvicorn src.backend.main:app --reload --host 0.0.0.0 --port 8000
