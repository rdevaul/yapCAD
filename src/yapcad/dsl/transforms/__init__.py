"""
DSL AST Transformation Framework.

Provides a framework for AST transformations that can be applied
before interpretation. This allows for optimizations like:
- Constant folding
- Dead code elimination
- Pattern unrolling
- Common subexpression caching

Usage:
    from yapcad.dsl.transforms import TransformPipeline, ConstantFoldTransform

    pipeline = TransformPipeline()
    pipeline.add(ConstantFoldTransform())

    optimized_module = pipeline.apply(module)
"""

from .base import (
    AstTransform,
    TreeTransform,
    TransformPipeline,
    IdentityTransform,
)

__all__ = [
    'AstTransform',
    'TreeTransform',
    'TransformPipeline',
    'IdentityTransform',
]
