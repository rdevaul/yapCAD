"""
AST transformation framework for DSL optimization.

Provides a base class for AST transformations that can be applied
before interpretation. This allows for future optimizations like:
- Constant folding
- Dead code elimination
- Pattern unrolling
- Common subexpression caching
"""

from abc import ABC, abstractmethod
from typing import List, Optional, TypeVar, Generic
from copy import deepcopy

from ..ast import (
    AstNode, Module, Command, Parameter,
    Statement, LetStatement, AssignmentStatement, RequireStatement,
    EmitStatement, ForStatement, ExpressionStatement, ReturnStatement,
    Block, PythonBlock,
    Expression, Literal, Identifier, BinaryOp, UnaryOp,
    FunctionCall, MethodCall, MemberAccess, IndexAccess,
    ListLiteral, ListComprehension, RangeExpr, DictLiteral,
    IfExpr, MatchExpr, LambdaExpr, PythonExpr,
)


class AstTransform(ABC):
    """
    Base class for AST transformations.

    Transforms are applied to a Module and return a (potentially modified)
    Module. Transforms can be composed in a pipeline.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Name of this transform for debugging/logging."""
        pass

    @abstractmethod
    def transform(self, module: Module) -> Module:
        """
        Apply this transform to a module.

        Args:
            module: The input module AST

        Returns:
            The transformed module (may be the same object or a copy)
        """
        pass


class TreeTransform(AstTransform):
    """
    A transform that walks the tree and can modify nodes.

    Subclasses override visit_* methods to transform specific node types.
    By default, nodes are copied unchanged.
    """

    def transform(self, module: Module) -> Module:
        """Transform a module by visiting all nodes."""
        return self.visit_module(module)

    def visit_module(self, node: Module) -> Module:
        """Visit a module node."""
        new_commands = [self.visit_command(cmd) for cmd in node.commands]
        return Module(
            name=node.name,
            uses=node.uses,  # Don't transform use statements
            commands=new_commands,
            span=node.span,
        )

    def visit_command(self, node: Command) -> Command:
        """Visit a command node."""
        new_body = [self.visit_statement(stmt) for stmt in node.body]
        return Command(
            name=node.name,
            params=node.params,
            return_type=node.return_type,
            body=new_body,
            span=node.span,
        )

    def visit_statement(self, node: Statement) -> Statement:
        """Visit a statement node."""
        if isinstance(node, LetStatement):
            return self.visit_let(node)
        elif isinstance(node, AssignmentStatement):
            return self.visit_assignment(node)
        elif isinstance(node, RequireStatement):
            return self.visit_require(node)
        elif isinstance(node, EmitStatement):
            return self.visit_emit(node)
        elif isinstance(node, ForStatement):
            return self.visit_for(node)
        elif isinstance(node, ExpressionStatement):
            return self.visit_expr_statement(node)
        elif isinstance(node, ReturnStatement):
            return self.visit_return(node)
        elif isinstance(node, Block):
            return self.visit_block(node)
        elif isinstance(node, PythonBlock):
            return self.visit_python_block(node)
        else:
            return node

    def visit_let(self, node: LetStatement) -> Statement:
        """Visit a let statement."""
        new_value = self.visit_expression(node.value)
        return LetStatement(
            name=node.name,
            type_annotation=node.type_annotation,
            value=new_value,
            span=node.span,
        )

    def visit_assignment(self, node: AssignmentStatement) -> Statement:
        """Visit an assignment statement."""
        new_value = self.visit_expression(node.value)
        return AssignmentStatement(
            target=node.target,
            value=new_value,
            span=node.span,
        )

    def visit_require(self, node: RequireStatement) -> Statement:
        """Visit a require statement."""
        new_condition = self.visit_expression(node.condition)
        return RequireStatement(
            condition=new_condition,
            message=node.message,
            span=node.span,
        )

    def visit_emit(self, node: EmitStatement) -> Statement:
        """Visit an emit statement."""
        new_value = self.visit_expression(node.value)
        new_metadata = None
        if node.metadata:
            new_metadata = {k: self.visit_expression(v) for k, v in node.metadata.items()}
        return EmitStatement(
            value=new_value,
            metadata=new_metadata,
            span=node.span,
        )

    def visit_for(self, node: ForStatement) -> Statement:
        """Visit a for statement."""
        new_iterable = self.visit_expression(node.iterable)
        new_body = [self.visit_statement(stmt) for stmt in node.body]
        return ForStatement(
            variable=node.variable,
            iterable=new_iterable,
            body=new_body,
            span=node.span,
        )

    def visit_expr_statement(self, node: ExpressionStatement) -> Statement:
        """Visit an expression statement."""
        new_expr = self.visit_expression(node.expression)
        return ExpressionStatement(expression=new_expr, span=node.span)

    def visit_return(self, node: ReturnStatement) -> Statement:
        """Visit a return statement."""
        new_value = self.visit_expression(node.value) if node.value else None
        return ReturnStatement(value=new_value, span=node.span)

    def visit_block(self, node: Block) -> Statement:
        """Visit a block."""
        new_stmts = [self.visit_statement(stmt) for stmt in node.statements]
        return Block(statements=new_stmts, span=node.span)

    def visit_python_block(self, node: PythonBlock) -> Statement:
        """Visit a Python block (no transformation by default)."""
        return node

    def visit_expression(self, node: Expression) -> Expression:
        """Visit an expression node."""
        if isinstance(node, Literal):
            return self.visit_literal(node)
        elif isinstance(node, Identifier):
            return self.visit_identifier(node)
        elif isinstance(node, BinaryOp):
            return self.visit_binary_op(node)
        elif isinstance(node, UnaryOp):
            return self.visit_unary_op(node)
        elif isinstance(node, FunctionCall):
            return self.visit_function_call(node)
        elif isinstance(node, MethodCall):
            return self.visit_method_call(node)
        elif isinstance(node, MemberAccess):
            return self.visit_member_access(node)
        elif isinstance(node, IndexAccess):
            return self.visit_index_access(node)
        elif isinstance(node, ListLiteral):
            return self.visit_list_literal(node)
        elif isinstance(node, ListComprehension):
            return self.visit_list_comprehension(node)
        elif isinstance(node, RangeExpr):
            return self.visit_range(node)
        elif isinstance(node, DictLiteral):
            return self.visit_dict_literal(node)
        elif isinstance(node, IfExpr):
            return self.visit_if_expr(node)
        elif isinstance(node, MatchExpr):
            return self.visit_match_expr(node)
        elif isinstance(node, LambdaExpr):
            return self.visit_lambda(node)
        elif isinstance(node, PythonExpr):
            return self.visit_python_expr(node)
        else:
            return node

    def visit_literal(self, node: Literal) -> Expression:
        return node

    def visit_identifier(self, node: Identifier) -> Expression:
        return node

    def visit_binary_op(self, node: BinaryOp) -> Expression:
        new_left = self.visit_expression(node.left)
        new_right = self.visit_expression(node.right)
        return BinaryOp(
            operator=node.operator,
            left=new_left,
            right=new_right,
            span=node.span,
        )

    def visit_unary_op(self, node: UnaryOp) -> Expression:
        new_operand = self.visit_expression(node.operand)
        return UnaryOp(
            operator=node.operator,
            operand=new_operand,
            span=node.span,
        )

    def visit_function_call(self, node: FunctionCall) -> Expression:
        new_args = [self.visit_expression(arg) for arg in node.arguments]
        return FunctionCall(
            name=node.name,
            arguments=new_args,
            span=node.span,
        )

    def visit_method_call(self, node: MethodCall) -> Expression:
        new_obj = self.visit_expression(node.object)
        new_args = [self.visit_expression(arg) for arg in node.arguments]
        return MethodCall(
            object=new_obj,
            method=node.method,
            arguments=new_args,
            span=node.span,
        )

    def visit_member_access(self, node: MemberAccess) -> Expression:
        new_obj = self.visit_expression(node.object)
        return MemberAccess(object=new_obj, member=node.member, span=node.span)

    def visit_index_access(self, node: IndexAccess) -> Expression:
        new_obj = self.visit_expression(node.object)
        new_index = self.visit_expression(node.index)
        return IndexAccess(object=new_obj, index=new_index, span=node.span)

    def visit_list_literal(self, node: ListLiteral) -> Expression:
        new_elements = [self.visit_expression(elem) for elem in node.elements]
        return ListLiteral(elements=new_elements, span=node.span)

    def visit_list_comprehension(self, node: ListComprehension) -> Expression:
        new_expr = self.visit_expression(node.expression)
        new_iterable = self.visit_expression(node.iterable)
        new_condition = self.visit_expression(node.condition) if node.condition else None
        return ListComprehension(
            expression=new_expr,
            variable=node.variable,
            iterable=new_iterable,
            condition=new_condition,
            span=node.span,
        )

    def visit_range(self, node: RangeExpr) -> Expression:
        new_start = self.visit_expression(node.start)
        new_end = self.visit_expression(node.end)
        return RangeExpr(start=new_start, end=new_end, span=node.span)

    def visit_dict_literal(self, node: DictLiteral) -> Expression:
        new_entries = {k: self.visit_expression(v) for k, v in node.entries.items()}
        return DictLiteral(entries=new_entries, span=node.span)

    def visit_if_expr(self, node: IfExpr) -> Expression:
        new_cond = self.visit_expression(node.condition)
        new_then = self.visit_expression(node.then_branch)
        new_else = self.visit_expression(node.else_branch) if node.else_branch else None
        return IfExpr(
            condition=new_cond,
            then_branch=new_then,
            else_branch=new_else,
            span=node.span,
        )

    def visit_match_expr(self, node: MatchExpr) -> Expression:
        new_subject = self.visit_expression(node.subject)
        # Arms would need pattern visiting too - simplified for now
        return MatchExpr(subject=new_subject, arms=node.arms, span=node.span)

    def visit_lambda(self, node: LambdaExpr) -> Expression:
        new_body = self.visit_expression(node.body)
        return LambdaExpr(parameters=node.parameters, body=new_body, span=node.span)

    def visit_python_expr(self, node: PythonExpr) -> Expression:
        return node


class TransformPipeline:
    """
    A pipeline of AST transforms to apply in sequence.
    """

    def __init__(self, transforms: List[AstTransform] = None):
        self.transforms = transforms or []

    def add(self, transform: AstTransform) -> "TransformPipeline":
        """Add a transform to the pipeline."""
        self.transforms.append(transform)
        return self

    def apply(self, module: Module) -> Module:
        """Apply all transforms in sequence."""
        result = module
        for transform in self.transforms:
            result = transform.transform(result)
        return result


# Example transform: Identity (does nothing, useful for testing)
class IdentityTransform(AstTransform):
    """Identity transform - returns module unchanged."""

    @property
    def name(self) -> str:
        return "identity"

    def transform(self, module: Module) -> Module:
        return module
