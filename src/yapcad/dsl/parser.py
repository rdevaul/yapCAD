"""
Recursive descent parser for the yapCAD DSL v2 (Pythonic Syntax).

Converts a token stream into an Abstract Syntax Tree (AST).
Supports Python-style indentation-based blocks.
"""

from typing import List, Optional, Callable, Any, Union
from .tokens import Token, TokenType, SourceSpan, SourceLocation, is_type_token
from .ast import (
    # Types
    TypeNode, SimpleType, GenericType, OptionalType,
    # Expressions
    Expression, Literal, Identifier, BinaryOp, UnaryOp,
    FunctionCall, MethodCall, MemberAccess, IndexAccess,
    ListLiteral, ListComprehension, ComprehensionClause, RangeExpr, ConditionalExpr, IfExpr, MatchExpr,
    MatchArm, Pattern, LiteralPattern, IdentifierPattern, WildcardPattern,
    LambdaExpr, PythonExpr, DictLiteral, ElifBranch,
    # Statements
    Statement, VarDecl, LetStatement, AssignmentStatement,
    AssertStatement, RequireStatement,
    EmitStatement, ForStatement, WhileStatement, IfStatement,
    ExpressionStatement, ReturnStatement, PassStatement,
    Block, PythonBlock,
    # Native blocks
    NativeBlock, NativeFunctionDecl, NativeFunction, Decorator,
    # Declarations
    Parameter, FunctionDef, Command, UseStatement, ExportUseStatement, Module,
)
from .errors import (
    ParserError,
    error_unexpected_token,
    error_unexpected_eof,
    error_invalid_expression,
    DiagnosticCollector,
    Diagnostic,
    ErrorSeverity,
)


class Parser:
    """
    Recursive descent parser for the yapCAD DSL v2 with Python-style indentation.

    Usage:
        parser = Parser(tokens)
        module = parser.parse_module()

    The parser implements standard precedence climbing for expressions:
        Lowest:  or
                 and
                 == !=
                 < > <= >=
                 + -
                 * / // %
        Highest: ** (power, right-associative)
                 unary (not -)
    """

    # Operator precedence levels (higher = tighter binding)
    PRECEDENCE = {
        TokenType.OR: 1,
        TokenType.AND: 2,
        TokenType.EQ: 3,
        TokenType.NE: 3,
        TokenType.LT: 4,
        TokenType.GT: 4,
        TokenType.LE: 4,
        TokenType.GE: 4,
        TokenType.PLUS: 5,
        TokenType.MINUS: 5,
        TokenType.STAR: 6,
        TokenType.SLASH: 6,
        TokenType.DOUBLE_SLASH: 6,
        TokenType.PERCENT: 6,
        TokenType.DOUBLE_STAR: 7,  # Power (right-associative)
    }

    # Right-associative operators
    RIGHT_ASSOCIATIVE = {TokenType.DOUBLE_STAR}

    def __init__(self, tokens: List[Token], filename: Optional[str] = None, source: Optional[str] = None):
        self.tokens = tokens
        self.filename = filename
        self.source = source  # Original source code for extracting raw text
        self.pos = 0
        self.diagnostics = DiagnosticCollector()

    # =========================================================================
    # Token Navigation
    # =========================================================================

    def _current(self) -> Token:
        """Get current token."""
        if self.pos >= len(self.tokens):
            return self.tokens[-1]  # EOF
        return self.tokens[self.pos]

    def _peek(self, offset: int = 0) -> Token:
        """Peek at token at current position + offset."""
        idx = self.pos + offset
        if idx >= len(self.tokens):
            return self.tokens[-1]  # EOF
        return self.tokens[idx]

    def _is_at_end(self) -> bool:
        """Check if at end of tokens."""
        return self._current().type == TokenType.EOF

    def _check(self, token_type: TokenType) -> bool:
        """Check if current token is of given type."""
        return self._current().type == token_type

    def _check_any(self, *token_types: TokenType) -> bool:
        """Check if current token is any of given types."""
        return self._current().type in token_types

    def _check_ahead(self, token_type: TokenType, offset: int = 1) -> bool:
        """Check if token at current position + offset is of given type."""
        return self._peek(offset).type == token_type

    def _advance(self) -> Token:
        """Consume and return current token."""
        token = self._current()
        if not self._is_at_end():
            self.pos += 1
        return token

    def _consume(self, token_type: TokenType, expected: str) -> Token:
        """Consume token of expected type, or raise error."""
        if self._check(token_type):
            return self._advance()
        self._error(expected)

    def _match(self, *token_types: TokenType) -> Optional[Token]:
        """Consume token if it matches any of the given types."""
        if self._current().type in token_types:
            return self._advance()
        return None

    def _skip_newlines(self) -> None:
        """Skip any NEWLINE tokens."""
        while self._check(TokenType.NEWLINE):
            self._advance()

    def _expect_newline_or_eof(self) -> None:
        """Expect a NEWLINE or EOF (end of statement in Pythonic syntax)."""
        if self._check(TokenType.EOF):
            return
        if self._check(TokenType.NEWLINE):
            self._advance()
            return
        # Also accept DEDENT (can follow statement at block end)
        if self._check(TokenType.DEDENT):
            return
        # Also accept RBRACE for brace blocks (statement before closing brace)
        if self._check(TokenType.RBRACE):
            return
        # Also accept SEMICOLON for legacy brace-block syntax
        if self._match(TokenType.SEMICOLON):
            return
        self._error("newline or end of file")

    def _error(self, expected: str) -> None:
        """Raise a parser error."""
        token = self._current()
        if token.type == TokenType.EOF:
            raise error_unexpected_eof(expected, token.span)
        raise error_unexpected_token(expected, token.type.name, token.span)

    def _span_from(self, start: Token) -> SourceSpan:
        """Create a span from start token to current position."""
        # Get the previous token's end position
        prev_pos = max(0, self.pos - 1)
        end_token = self.tokens[prev_pos]
        return SourceSpan(start.span.start, end_token.span.end)

    # =========================================================================
    # Type Parsing
    # =========================================================================

    def _parse_type(self) -> TypeNode:
        """Parse a type annotation."""
        start = self._current()

        # Check for type keywords
        if is_type_token(self._current().type):
            type_name = self._advance().value
        elif self._check(TokenType.IDENTIFIER):
            type_name = self._advance().value
        else:
            self._error("type")

        # Check for generic parameters (e.g., list[T] or list<T>)
        if self._match(TokenType.LBRACKET):
            # Python-style: list[T]
            type_args = [self._parse_type()]
            while self._match(TokenType.COMMA):
                type_args.append(self._parse_type())
            self._consume(TokenType.RBRACKET, "']'")
            result = GenericType(
                span=self._span_from(start),
                name=type_name,
                type_args=type_args
            )
        elif self._match(TokenType.LT):
            # Legacy style: list<T>
            type_args = [self._parse_type()]
            while self._match(TokenType.COMMA):
                type_args.append(self._parse_type())
            self._consume(TokenType.GT, "'>'")
            result = GenericType(
                span=self._span_from(start),
                name=type_name,
                type_args=type_args
            )
        else:
            result = SimpleType(span=self._span_from(start), name=type_name)

        # Check for optional marker
        if self._match(TokenType.QUESTION):
            result = OptionalType(span=self._span_from(start), inner=result)

        return result

    # =========================================================================
    # Expression Parsing (Precedence Climbing)
    # =========================================================================

    def _parse_expression(self) -> Expression:
        """Parse an expression.

        Handles ternary conditional expressions at the lowest precedence:
            value if condition else other_value
        """
        return self._parse_conditional_expr()

    def _parse_conditional_expr(self) -> Expression:
        """Parse a conditional expression (ternary if-else).

        Syntax: expr if condition else expr

        This has the lowest precedence of all expression forms.
        """
        # Parse the true branch (what comes before 'if')
        true_branch = self._parse_binary_expr(0)

        # Check for 'if' keyword indicating ternary conditional
        if not self._check(TokenType.IF):
            return true_branch

        # Don't consume 'if' if followed by ':' (that's a block-level if statement)
        # In expression context, we're looking for: expr if cond else expr
        # Save position to backtrack if this is actually a block if
        saved_pos = self.pos
        self._advance()  # consume 'if'

        # Parse the condition
        condition = self._parse_binary_expr(0)

        # Must have 'else' for ternary conditional
        if not self._check(TokenType.ELSE):
            # This might be a list comprehension 'if' or something else
            # Backtrack and return just the first expression
            self.pos = saved_pos
            return true_branch

        self._advance()  # consume 'else'

        # Parse the false branch (recursively to allow chaining)
        false_branch = self._parse_conditional_expr()

        return ConditionalExpr(
            span=SourceSpan(true_branch.span.start, false_branch.span.end),
            condition=condition,
            true_branch=true_branch,
            false_branch=false_branch
        )

    def _parse_binary_expr(self, min_precedence: int) -> Expression:
        """Parse binary expressions with precedence climbing."""
        left = self._parse_unary_expr()

        while True:
            op_token = self._current()
            precedence = self.PRECEDENCE.get(op_token.type)

            if precedence is None or precedence < min_precedence:
                break

            self._advance()  # consume operator

            # Right-associative operators use same precedence, others use precedence + 1
            next_precedence = precedence if op_token.type in self.RIGHT_ASSOCIATIVE else precedence + 1
            right = self._parse_binary_expr(next_precedence)

            left = BinaryOp(
                span=SourceSpan(left.span.start, right.span.end),
                left=left,
                operator=op_token.type,
                right=right
            )

        return left

    def _parse_unary_expr(self) -> Expression:
        """Parse unary expressions (not, -)."""
        if self._check_any(TokenType.NOT, TokenType.MINUS):
            op = self._advance()
            operand = self._parse_unary_expr()
            return UnaryOp(
                span=SourceSpan(op.span.start, operand.span.end),
                operator=op.type,
                operand=operand
            )

        return self._parse_postfix_expr()

    def _parse_postfix_expr(self) -> Expression:
        """Parse postfix expressions (calls, member access, indexing)."""
        expr = self._parse_primary_expr()

        while True:
            if self._check(TokenType.LPAREN):
                # Function call
                expr = self._parse_call(expr)
            elif self._check(TokenType.DOT):
                # Member access or method call
                self._advance()  # consume '.'
                member = self._consume(TokenType.IDENTIFIER, "identifier").value

                if self._check(TokenType.LPAREN):
                    # Method call
                    args, named_args = self._parse_arguments()
                    expr = MethodCall(
                        span=self._span_from(self.tokens[self.pos - 1]),
                        object=expr,
                        method=member,
                        arguments=args,
                        named_arguments=named_args
                    )
                else:
                    # Member access
                    expr = MemberAccess(
                        span=self._span_from(self.tokens[self.pos - 1]),
                        object=expr,
                        member=member
                    )
            elif self._check(TokenType.LBRACKET):
                # Index access
                self._advance()  # consume '['
                index = self._parse_expression()
                self._consume(TokenType.RBRACKET, "']'")
                expr = IndexAccess(
                    span=self._span_from(self.tokens[self.pos - 1]),
                    object=expr,
                    index=index
                )
            else:
                break

        return expr

    def _parse_call(self, callee: Expression) -> FunctionCall:
        """Parse function call arguments."""
        args, named_args = self._parse_arguments()
        return FunctionCall(
            span=SourceSpan(callee.span.start, self.tokens[self.pos - 1].span.end),
            callee=callee,
            arguments=args,
            named_arguments=named_args
        )

    def _parse_arguments(self) -> tuple[List[Expression], dict[str, Expression]]:
        """Parse argument list, returns (positional, named)."""
        self._consume(TokenType.LPAREN, "'('")

        args = []
        named_args = {}

        if not self._check(TokenType.RPAREN):
            # First argument
            self._parse_argument(args, named_args)

            while self._match(TokenType.COMMA):
                if self._check(TokenType.RPAREN):
                    break  # Allow trailing comma
                self._parse_argument(args, named_args)

        self._consume(TokenType.RPAREN, "')'")
        return args, named_args

    def _parse_argument(self, args: List[Expression],
                        named_args: dict[str, Expression]) -> None:
        """Parse a single argument (positional or named)."""
        # Check for named argument: name=value
        if (self._check(TokenType.IDENTIFIER) and
            self._peek(1).type == TokenType.ASSIGN):
            name = self._advance().value
            self._advance()  # consume '='
            value = self._parse_expression()
            named_args[name] = value
        else:
            args.append(self._parse_expression())

    def _parse_primary_expr(self) -> Expression:
        """Parse primary expressions (literals, identifiers, grouped, etc.)."""
        token = self._current()

        # Literals
        if token.type == TokenType.INT_LITERAL:
            self._advance()
            return Literal(
                span=token.span,
                value=token.value,
                literal_type=TokenType.INT_LITERAL
            )

        if token.type == TokenType.FLOAT_LITERAL:
            self._advance()
            return Literal(
                span=token.span,
                value=token.value,
                literal_type=TokenType.FLOAT_LITERAL
            )

        if token.type == TokenType.STRING_LITERAL:
            self._advance()
            return Literal(
                span=token.span,
                value=token.value,
                literal_type=TokenType.STRING_LITERAL
            )

        if token.type == TokenType.BOOL_LITERAL:
            self._advance()
            return Literal(
                span=token.span,
                value=token.value,
                literal_type=TokenType.BOOL_LITERAL
            )

        # Identifiers
        if token.type == TokenType.IDENTIFIER:
            self._advance()
            return Identifier(span=token.span, name=token.value)

        # Type names can also be used as constructors (e.g., point(...))
        if is_type_token(token.type):
            self._advance()
            return Identifier(span=token.span, name=token.value)

        # Grouped expression or lambda
        if token.type == TokenType.LPAREN:
            return self._parse_grouped_or_lambda()

        # List literal or comprehension
        if token.type == TokenType.LBRACKET:
            return self._parse_list_literal()

        # Dict literal
        if token.type == TokenType.LBRACE:
            return self._parse_dict_literal()

        # If expression
        if token.type == TokenType.IF:
            return self._parse_if_expr()

        # Match expression
        if token.type == TokenType.MATCH:
            return self._parse_match_expr()

        # Range expression handled in binary (but check for standalone)
        if token.type == TokenType.RANGE:
            self._error("expression before '..'")

        self._error("expression")

    def _parse_grouped_or_lambda(self) -> Expression:
        """Parse grouped expression or lambda."""
        start = self._advance()  # consume '('

        # Check for lambda: (params) => expr
        # This is a simplified check - full parsing would need lookahead
        if self._check(TokenType.RPAREN):
            self._advance()
            if self._match(TokenType.DOUBLE_ARROW):
                body = self._parse_expression()
                return LambdaExpr(span=self._span_from(start), parameters=[], body=body)
            self._error("expression")

        # Check if it looks like lambda parameters
        if self._check(TokenType.IDENTIFIER):
            # Could be lambda or grouped expression - need to look ahead
            saved_pos = self.pos
            param_names = [self._advance().value]

            while self._match(TokenType.COMMA):
                if self._check(TokenType.IDENTIFIER):
                    param_names.append(self._advance().value)
                else:
                    break

            if self._match(TokenType.RPAREN) and self._match(TokenType.DOUBLE_ARROW):
                body = self._parse_expression()
                return LambdaExpr(span=self._span_from(start), parameters=param_names, body=body)

            # Not a lambda, restore position and parse as grouped expression
            self.pos = saved_pos

        expr = self._parse_expression()
        self._consume(TokenType.RPAREN, "')'")
        return expr

    def _parse_list_literal(self) -> Expression:
        """Parse list literal or list comprehension.

        Supports nested comprehensions:
            [expr for x in xs]                    # single clause
            [expr for x in xs if cond]            # with condition
            [expr for x in xs for y in ys]        # nested (multiple clauses)
            [expr for x in xs if c1 for y in ys if c2]  # with conditions
        """
        start = self._advance()  # consume '['

        if self._check(TokenType.RBRACKET):
            self._advance()
            return ListLiteral(span=self._span_from(start), elements=[])

        first = self._parse_expression()

        # Check for list comprehension: [expr for x in iterable ...]
        if self._check(TokenType.FOR):
            clauses: List[ComprehensionClause] = []

            # Parse one or more for clauses
            while self._check(TokenType.FOR):
                clause_start = self._current()
                self._advance()  # consume 'for'
                var = self._consume(TokenType.IDENTIFIER, "identifier").value
                self._consume(TokenType.IN, "'in'")
                iterable = self._parse_expression()

                # Parse zero or more 'if' conditions for this clause
                # Each 'if' is followed by a condition expression. The condition
                # expression will not consume 'for' (it's a keyword, not identifier),
                # so the outer while loop correctly stops at the next 'for'.
                conditions: List[Expression] = []
                while self._check(TokenType.IF):
                    self._advance()  # consume 'if'
                    conditions.append(self._parse_expression())

                clauses.append(ComprehensionClause(
                    span=self._span_from(clause_start),
                    variable=var,
                    iterable=iterable,
                    conditions=conditions
                ))

            self._consume(TokenType.RBRACKET, "']'")
            return ListComprehension(
                span=self._span_from(start),
                element_expr=first,
                clauses=clauses
            )

        # Check for range: [start..end]
        if self._check(TokenType.RANGE):
            self._advance()  # consume '..'
            end = self._parse_expression()
            self._consume(TokenType.RBRACKET, "']'")
            return RangeExpr(span=self._span_from(start), start=first, end=end)

        # Regular list literal
        elements = [first]
        while self._match(TokenType.COMMA):
            if self._check(TokenType.RBRACKET):
                break  # Allow trailing comma
            elements.append(self._parse_expression())

        self._consume(TokenType.RBRACKET, "']'")
        return ListLiteral(span=self._span_from(start), elements=elements)

    def _parse_dict_literal(self) -> DictLiteral:
        """Parse a dictionary literal { key: value, ... } or {"key": value}."""
        start = self._advance()  # consume '{'
        entries = {}

        if not self._check(TokenType.RBRACE):
            # First entry - key can be identifier or string
            if self._check(TokenType.STRING_LITERAL):
                key = self._advance().value
            else:
                key = self._consume(TokenType.IDENTIFIER, "identifier or string").value
            self._consume(TokenType.COLON, "':'")
            value = self._parse_expression()
            entries[key] = value

            while self._match(TokenType.COMMA):
                if self._check(TokenType.RBRACE):
                    break
                if self._check(TokenType.STRING_LITERAL):
                    key = self._advance().value
                else:
                    key = self._consume(TokenType.IDENTIFIER, "identifier or string").value
                self._consume(TokenType.COLON, "':'")
                value = self._parse_expression()
                entries[key] = value

        self._consume(TokenType.RBRACE, "'}'")
        return DictLiteral(span=self._span_from(start), entries=entries)

    def _parse_if_expr(self) -> IfExpr:
        """Parse an if expression (for use in expressions)."""
        start = self._advance()  # consume 'if'
        condition = self._parse_expression()
        self._consume(TokenType.COLON, "':'")
        then_branch = self._parse_block()

        elif_branches = []
        while self._check(TokenType.ELIF):
            self._advance()  # consume 'elif'
            elif_cond = self._parse_expression()
            self._consume(TokenType.COLON, "':'")
            elif_body = self._parse_block()
            elif_branches.append(ElifBranch(
                span=self._span_from(start),
                condition=elif_cond,
                body=elif_body
            ))

        else_branch = None
        if self._match(TokenType.ELSE):
            self._consume(TokenType.COLON, "':'")
            else_branch = self._parse_block()

        return IfExpr(
            span=self._span_from(start),
            condition=condition,
            then_branch=then_branch,
            elif_branches=elif_branches,
            else_branch=else_branch
        )

    def _parse_match_expr(self) -> MatchExpr:
        """Parse a match expression."""
        start = self._advance()  # consume 'match'
        subject = self._parse_expression()
        self._consume(TokenType.COLON, "':'")
        self._consume(TokenType.NEWLINE, "newline")
        self._consume(TokenType.INDENT, "indented block")

        arms = []
        while not self._check(TokenType.DEDENT) and not self._is_at_end():
            self._skip_newlines()
            if self._check(TokenType.DEDENT):
                break
            arms.append(self._parse_match_arm())

        self._consume(TokenType.DEDENT, "end of match block")
        return MatchExpr(span=self._span_from(start), subject=subject, arms=arms)

    def _parse_match_arm(self) -> MatchArm:
        """Parse a single match arm."""
        start = self._current()
        # Optional 'case' keyword for Python-style match
        self._match(TokenType.IDENTIFIER)  # Skip 'case' if present
        pattern = self._parse_pattern()
        self._consume(TokenType.COLON, "':'")
        body = self._parse_expression()
        self._skip_newlines()

        return MatchArm(span=self._span_from(start), pattern=pattern, body=body)

    def _parse_pattern(self) -> Pattern:
        """Parse a match pattern."""
        token = self._current()

        # Wildcard
        if token.type == TokenType.UNDERSCORE:
            self._advance()
            return WildcardPattern(span=token.span)

        # String literal
        if token.type == TokenType.STRING_LITERAL:
            self._advance()
            lit = Literal(span=token.span, value=token.value,
                          literal_type=TokenType.STRING_LITERAL)
            return LiteralPattern(span=token.span, value=lit)

        # Numeric literal
        if token.type in (TokenType.INT_LITERAL, TokenType.FLOAT_LITERAL):
            self._advance()
            lit = Literal(span=token.span, value=token.value,
                          literal_type=token.type)
            return LiteralPattern(span=token.span, value=lit)

        # Boolean literal
        if token.type == TokenType.BOOL_LITERAL:
            self._advance()
            lit = Literal(span=token.span, value=token.value,
                          literal_type=TokenType.BOOL_LITERAL)
            return LiteralPattern(span=token.span, value=lit)

        # Identifier (binding)
        if token.type == TokenType.IDENTIFIER:
            self._advance()
            return IdentifierPattern(span=token.span, name=token.value)

        self._error("pattern")

    # =========================================================================
    # Statement Parsing
    # =========================================================================

    def _parse_statement(self) -> Statement:
        """Parse a statement."""
        self._skip_newlines()
        token = self._current()

        # Legacy 'let' keyword (deprecated but supported)
        if token.type == TokenType.LET:
            return self._parse_let_statement()

        # Legacy 'require' keyword (deprecated, use assert)
        if token.type == TokenType.REQUIRE:
            return self._parse_require_statement()

        # Assert statement
        if token.type == TokenType.ASSERT:
            return self._parse_assert_statement()

        # Emit statement
        if token.type == TokenType.EMIT:
            return self._parse_emit_statement()

        # For loop
        if token.type == TokenType.FOR:
            return self._parse_for_statement()

        # While loop - DEPRECATED, give helpful error
        if token.type == TokenType.WHILE:
            diag = Diagnostic(
                code="E104",
                message="'while' loops are not supported (removed for static verifiability). "
                        "Use 'for i in range(max_iterations)' with early return instead.",
                severity=ErrorSeverity.ERROR,
                span=token.span,
            )
            raise ParserError(diag)

        # If statement
        if token.type == TokenType.IF:
            return self._parse_if_statement()

        # Pass statement
        if token.type == TokenType.PASS:
            return self._parse_pass_statement()

        # Return statement
        if token.type == TokenType.RETURN:
            return self._parse_return_statement()

        # Legacy python block
        if token.type == TokenType.PYTHON:
            return self._parse_python_block()

        # Variable declaration or expression/assignment
        # Check for: name: type = value or name = value
        if self._check(TokenType.IDENTIFIER):
            # Look ahead to determine if this is a var decl or expression
            if self._peek(1).type == TokenType.COLON:
                return self._parse_var_decl()
            # Check for simple assignment that creates new variable
            if self._peek(1).type == TokenType.ASSIGN:
                # This could be assignment or new variable - we'll treat as new var
                return self._parse_simple_var_decl()

        # Expression statement or assignment
        expr = self._parse_expression()

        # Check for assignment
        if self._match(TokenType.ASSIGN):
            value = self._parse_expression()
            self._expect_newline_or_eof()
            return AssignmentStatement(
                span=SourceSpan(expr.span.start, value.span.end),
                target=expr,
                value=value
            )

        self._expect_newline_or_eof()
        return ExpressionStatement(span=expr.span, expression=expr)

    def _parse_var_decl(self) -> VarDecl:
        """Parse a variable declaration with type annotation: name: type = value."""
        start = self._current()
        name = self._consume(TokenType.IDENTIFIER, "identifier").value
        self._consume(TokenType.COLON, "':'")
        type_annotation = self._parse_type()

        initializer = None
        if self._match(TokenType.ASSIGN):
            initializer = self._parse_expression()

        self._expect_newline_or_eof()
        return VarDecl(
            span=self._span_from(start),
            name=name,
            type_annotation=type_annotation,
            initializer=initializer
        )

    def _parse_simple_var_decl(self) -> VarDecl:
        """Parse a simple variable declaration: name = value (type inferred)."""
        start = self._current()
        name = self._consume(TokenType.IDENTIFIER, "identifier").value
        self._consume(TokenType.ASSIGN, "'='")
        initializer = self._parse_expression()
        self._expect_newline_or_eof()

        return VarDecl(
            span=self._span_from(start),
            name=name,
            type_annotation=None,
            initializer=initializer
        )

    def _parse_let_statement(self) -> LetStatement:
        """Parse a legacy let statement (deprecated)."""
        start = self._advance()  # consume 'let'
        name = self._consume(TokenType.IDENTIFIER, "identifier").value

        type_annotation = None
        if self._match(TokenType.COLON):
            type_annotation = self._parse_type()

        self._consume(TokenType.ASSIGN, "'='")
        initializer = self._parse_expression()

        # Accept either semicolon (legacy) or newline
        if not self._match(TokenType.SEMICOLON):
            self._expect_newline_or_eof()

        return LetStatement(
            span=self._span_from(start),
            name=name,
            type_annotation=type_annotation,
            initializer=initializer
        )

    def _parse_require_statement(self) -> RequireStatement:
        """Parse a legacy require statement (deprecated, use assert)."""
        start = self._advance()  # consume 'require'
        condition = self._parse_expression()

        message = None
        if self._match(TokenType.COMMA):
            message = self._parse_expression()

        # Accept either semicolon (legacy) or newline
        if not self._match(TokenType.SEMICOLON):
            self._expect_newline_or_eof()

        return RequireStatement(
            span=self._span_from(start),
            condition=condition,
            message=message
        )

    def _parse_assert_statement(self) -> AssertStatement:
        """Parse an assert statement."""
        start = self._advance()  # consume 'assert'
        condition = self._parse_expression()

        message = None
        if self._match(TokenType.COMMA):
            message = self._parse_expression()

        self._expect_newline_or_eof()
        return AssertStatement(
            span=self._span_from(start),
            condition=condition,
            message=message
        )

    def _parse_emit_statement(self) -> EmitStatement:
        """Parse an emit statement with optional kwargs metadata."""
        start = self._advance()  # consume 'emit'
        value = self._parse_expression()

        # Parse optional metadata as keyword arguments
        metadata = {}
        while self._match(TokenType.COMMA):
            if self._check(TokenType.IDENTIFIER) and self._peek(1).type == TokenType.ASSIGN:
                key = self._advance().value
                self._advance()  # consume '='
                metadata[key] = self._parse_expression()
            else:
                break

        # Legacy: accept 'with { ... }' syntax
        if self._match(TokenType.WITH):
            dict_lit = self._parse_dict_literal()
            metadata.update(dict_lit.entries)

        # Accept either semicolon (legacy) or newline
        if not self._match(TokenType.SEMICOLON):
            self._expect_newline_or_eof()

        return EmitStatement(
            span=self._span_from(start),
            value=value,
            metadata=metadata
        )

    def _parse_for_statement(self) -> ForStatement:
        """Parse a for statement."""
        start = self._advance()  # consume 'for'
        variable = self._consume(TokenType.IDENTIFIER, "identifier").value
        self._consume(TokenType.IN, "'in'")
        iterable = self._parse_expression()
        self._consume(TokenType.COLON, "':'")
        body = self._parse_block()

        return ForStatement(
            span=self._span_from(start),
            variable=variable,
            iterable=iterable,
            body=body
        )

    # NOTE: _parse_while_statement removed - while loops not supported for static verifiability

    def _parse_if_statement(self) -> IfStatement:
        """Parse an if statement."""
        start = self._advance()  # consume 'if'
        condition = self._parse_expression()
        self._consume(TokenType.COLON, "':'")
        then_branch = self._parse_block()

        elif_branches = []
        while self._check(TokenType.ELIF):
            self._advance()  # consume 'elif'
            elif_cond = self._parse_expression()
            self._consume(TokenType.COLON, "':'")
            elif_body = self._parse_block()
            elif_branches.append(ElifBranch(
                span=self._span_from(start),
                condition=elif_cond,
                body=elif_body
            ))

        else_branch = None
        if self._match(TokenType.ELSE):
            self._consume(TokenType.COLON, "':'")
            else_branch = self._parse_block()

        return IfStatement(
            span=self._span_from(start),
            condition=condition,
            then_branch=then_branch,
            elif_branches=elif_branches,
            else_branch=else_branch
        )

    def _parse_pass_statement(self) -> PassStatement:
        """Parse a pass statement."""
        start = self._advance()  # consume 'pass'
        self._expect_newline_or_eof()
        return PassStatement(span=self._span_from(start))

    def _parse_python_block(self) -> PythonBlock:
        """Parse a legacy python block."""
        start = self._advance()  # consume 'python'
        self._consume(TokenType.LBRACE, "'{'")

        # Capture everything until matching '}'
        brace_depth = 1
        code_start = self._current().span.start.offset
        code_tokens = []

        while not self._is_at_end() and brace_depth > 0:
            token = self._advance()
            if token.type == TokenType.LBRACE:
                brace_depth += 1
            elif token.type == TokenType.RBRACE:
                brace_depth -= 1
            if brace_depth > 0:
                code_tokens.append(token)

        code = " ".join(t.lexeme for t in code_tokens)
        return PythonBlock(span=self._span_from(start), code=code)

    def _parse_return_statement(self) -> ReturnStatement:
        """Parse a return statement."""
        start = self._advance()  # consume 'return'

        value = None
        # Check if there's a return value (not just 'return' on its own)
        if not self._check_any(TokenType.NEWLINE, TokenType.DEDENT, TokenType.EOF):
            value = self._parse_expression()

        # Legacy: 'return value as type' syntax
        if self._match(TokenType.AS):
            # Skip the type annotation for legacy support
            self._parse_type()

        # Accept either semicolon (legacy) or newline
        if not self._match(TokenType.SEMICOLON):
            self._expect_newline_or_eof()

        return ReturnStatement(span=self._span_from(start), value=value)

    def _parse_block(self) -> Block:
        """Parse a block of statements (indentation-based)."""
        # Skip newlines (including blank lines) after the colon
        self._skip_newlines()

        # For legacy brace-based blocks
        if self._check(TokenType.LBRACE):
            return self._parse_brace_block()

        # Python-style indentation block
        start = self._consume(TokenType.INDENT, "indented block")
        statements = []

        while not self._check(TokenType.DEDENT) and not self._is_at_end():
            self._skip_newlines()
            if self._check(TokenType.DEDENT):
                break
            statements.append(self._parse_statement())

        self._consume(TokenType.DEDENT, "end of block")

        # Check if last statement is an expression (for expression-valued blocks)
        final_expr = None
        if statements and isinstance(statements[-1], ExpressionStatement):
            final_expr = statements[-1].expression
            statements = statements[:-1]

        return Block(
            span=self._span_from(start),
            statements=statements,
            final_expression=final_expr
        )

    def _parse_brace_block(self) -> Block:
        """Parse a legacy brace-delimited block."""
        start = self._consume(TokenType.LBRACE, "'{'")
        statements = []

        while not self._check(TokenType.RBRACE) and not self._is_at_end():
            self._skip_newlines()
            if self._check(TokenType.RBRACE):
                break
            statements.append(self._parse_statement())

        self._consume(TokenType.RBRACE, "'}'")

        final_expr = None
        if statements and isinstance(statements[-1], ExpressionStatement):
            final_expr = statements[-1].expression
            statements = statements[:-1]

        return Block(
            span=self._span_from(start),
            statements=statements,
            final_expression=final_expr
        )

    # =========================================================================
    # Declaration Parsing
    # =========================================================================

    def _parse_parameter(self) -> Parameter:
        """Parse a function parameter."""
        start = self._current()
        name = self._consume(TokenType.IDENTIFIER, "parameter name").value

        type_annotation = None
        if self._match(TokenType.COLON):
            type_annotation = self._parse_type()

        default_value = None
        if self._match(TokenType.ASSIGN):
            default_value = self._parse_expression()

        return Parameter(
            span=self._span_from(start),
            name=name,
            type_annotation=type_annotation,
            default_value=default_value
        )

    def _parse_decorator(self) -> Decorator:
        """Parse a decorator."""
        start = self._advance()  # consume '@'
        name = self._consume(TokenType.IDENTIFIER, "decorator name").value

        arguments = []
        if self._match(TokenType.LPAREN):
            if not self._check(TokenType.RPAREN):
                arguments.append(self._parse_expression())
                while self._match(TokenType.COMMA):
                    arguments.append(self._parse_expression())
            self._consume(TokenType.RPAREN, "')'")

        self._expect_newline_or_eof()
        return Decorator(span=self._span_from(start), name=name, arguments=arguments)

    def _parse_function_def(self, decorators: List[Decorator] = None) -> FunctionDef:
        """Parse a function definition (def keyword)."""
        start = self._advance()  # consume 'def'
        name = self._consume(TokenType.IDENTIFIER, "function name").value

        self._consume(TokenType.LPAREN, "'('")
        parameters = []
        if not self._check(TokenType.RPAREN):
            parameters.append(self._parse_parameter())
            while self._match(TokenType.COMMA):
                if self._check(TokenType.RPAREN):
                    break
                parameters.append(self._parse_parameter())
        self._consume(TokenType.RPAREN, "')'")

        return_type = None
        if self._match(TokenType.ARROW):
            return_type = self._parse_type()

        self._consume(TokenType.COLON, "':'")
        body = self._parse_block()

        return FunctionDef(
            span=self._span_from(start),
            name=name,
            parameters=parameters,
            return_type=return_type,
            body=body,
            decorators=decorators or []
        )

    def _parse_command(self) -> Command:
        """Parse a legacy command definition."""
        start = self._advance()  # consume 'command'
        name = self._consume(TokenType.IDENTIFIER, "command name").value

        self._consume(TokenType.LPAREN, "'('")
        parameters = []
        if not self._check(TokenType.RPAREN):
            parameters.append(self._parse_parameter())
            while self._match(TokenType.COMMA):
                if self._check(TokenType.RPAREN):
                    break
                parameters.append(self._parse_parameter())
        self._consume(TokenType.RPAREN, "')'")

        self._consume(TokenType.ARROW, "'->'")
        return_type = self._parse_type()

        # Consume colon for Python-style syntax, or fall through to brace block
        if self._check(TokenType.COLON):
            self._advance()

        body = self._parse_block()

        return Command(
            span=self._span_from(start),
            name=name,
            parameters=parameters,
            return_type=return_type,
            body=body,
            decorators=[]
        )

    def _parse_module_path_component(self) -> str:
        """Parse a module path component."""
        token = self._current()
        if token.type == TokenType.IDENTIFIER:
            self._advance()
            return token.value
        if token.value is not None and isinstance(token.value, str):
            self._advance()
            return token.value
        self._error("module name")

    def _parse_use_statement(self) -> Union[UseStatement, ExportUseStatement]:
        """Parse a use statement."""
        start = self._current()
        is_export = self._match(TokenType.EXPORT)
        self._consume(TokenType.USE, "'use'")

        path = [self._parse_module_path_component()]
        while self._match(TokenType.DOT):
            path.append(self._parse_module_path_component())

        alias = None
        if not is_export and self._match(TokenType.AS):
            alias = self._consume(TokenType.IDENTIFIER, "alias").value

        # Accept either semicolon (legacy) or newline
        if not self._match(TokenType.SEMICOLON):
            self._expect_newline_or_eof()

        if is_export:
            return ExportUseStatement(span=self._span_from(start), module_path=path)
        return UseStatement(span=self._span_from(start), module_path=path, alias=alias)

    def _parse_native_function_decl(self) -> NativeFunctionDecl:
        """Parse a native function declaration (legacy exports block)."""
        start = self._advance()  # consume 'fn'
        name = self._consume(TokenType.IDENTIFIER, "function name").value

        self._consume(TokenType.LPAREN, "'('")
        parameters = []
        if not self._check(TokenType.RPAREN):
            parameters.append(self._parse_parameter())
            while self._match(TokenType.COMMA):
                if self._check(TokenType.RPAREN):
                    break
                parameters.append(self._parse_parameter())
        self._consume(TokenType.RPAREN, "')'")

        self._consume(TokenType.ARROW, "'->'")
        return_type = self._parse_type()

        self._consume(TokenType.SEMICOLON, "';'")

        return NativeFunctionDecl(
            span=self._span_from(start),
            name=name,
            parameters=parameters,
            return_type=return_type
        )

    def _parse_native_block(self) -> NativeBlock:
        """Parse a legacy native python block."""
        start = self._advance()  # consume 'native'
        self._consume(TokenType.PYTHON, "'python'")
        self._consume(TokenType.LBRACE, "'{'")

        brace_depth = 1
        code_start_offset = self._current().span.start.offset

        while not self._is_at_end() and brace_depth > 0:
            token = self._advance()
            if token.type == TokenType.LBRACE:
                brace_depth += 1
            elif token.type == TokenType.RBRACE:
                brace_depth -= 1

        closing_brace = self.tokens[self.pos - 1]
        code_end_offset = closing_brace.span.start.offset

        if self.source is not None:
            code = self.source[code_start_offset:code_end_offset].strip()
        else:
            code = ""

        self._consume(TokenType.EXPORTS, "'exports'")
        self._consume(TokenType.LBRACE, "'{'")

        exports = []
        while self._check(TokenType.FN) and not self._is_at_end():
            exports.append(self._parse_native_function_decl())

        self._consume(TokenType.RBRACE, "'}'")

        return NativeBlock(
            span=self._span_from(start),
            code=code,
            exports=exports
        )

    def _parse_native_function(self, decorators: List[Decorator]) -> NativeFunction:
        """Parse a @native decorated function."""
        start = self._current()
        self._advance()  # consume 'def'
        name = self._consume(TokenType.IDENTIFIER, "function name").value

        self._consume(TokenType.LPAREN, "'('")
        parameters = []
        if not self._check(TokenType.RPAREN):
            parameters.append(self._parse_parameter())
            while self._match(TokenType.COMMA):
                if self._check(TokenType.RPAREN):
                    break
                parameters.append(self._parse_parameter())
        self._consume(TokenType.RPAREN, "')'")

        self._consume(TokenType.ARROW, "'->'")
        return_type = self._parse_type()

        self._consume(TokenType.COLON, "':'")

        # Parse the function body as Python code
        # For now, extract from source using indentation
        if self._check(TokenType.NEWLINE):
            self._advance()

        # Capture the indented Python code
        python_code = ""
        if self._check(TokenType.INDENT):
            self._advance()
            code_start = self._current().span.start.offset

            # Read until DEDENT
            while not self._check(TokenType.DEDENT) and not self._is_at_end():
                self._advance()

            code_end = self._current().span.start.offset
            if self.source:
                python_code = self.source[code_start:code_end].strip()

            self._consume(TokenType.DEDENT, "end of native function")

        return NativeFunction(
            span=self._span_from(start),
            name=name,
            parameters=parameters,
            return_type=return_type,
            python_code=python_code
        )

    def parse_module(self) -> Module:
        """Parse a complete module."""
        start = self._current()
        self._skip_newlines()

        # Optional module declaration
        name = None
        if self._match(TokenType.MODULE):
            name = self._consume(TokenType.IDENTIFIER, "module name").value
            # Accept either semicolon (legacy) or newline
            if not self._match(TokenType.SEMICOLON):
                self._expect_newline_or_eof()

        self._skip_newlines()

        # Use statements
        uses = []
        while self._check_any(TokenType.USE, TokenType.EXPORT):
            if self._check(TokenType.EXPORT) and self._peek(1).type != TokenType.USE:
                break  # Not an export use statement
            uses.append(self._parse_use_statement())
            self._skip_newlines()

        # Native blocks and function definitions
        native_blocks = []
        native_functions = []
        functions = []

        while not self._is_at_end():
            self._skip_newlines()
            if self._is_at_end():
                break

            # Decorators
            decorators = []
            while self._check(TokenType.AT):
                decorators.append(self._parse_decorator())
                self._skip_newlines()

            if self._is_at_end():
                break

            # Check what follows
            if self._check(TokenType.NATIVE):
                # Legacy native block
                native_blocks.append(self._parse_native_block())
            elif self._check(TokenType.DEF):
                # Check if this is a @native decorated function
                is_native = any(d.name == "native" for d in decorators)
                if is_native:
                    native_functions.append(self._parse_native_function(decorators))
                else:
                    functions.append(self._parse_function_def(decorators))
            elif self._check(TokenType.COMMAND):
                # Legacy command
                functions.append(self._parse_command())
            else:
                self._error("function definition or end of file")

            self._skip_newlines()

        return Module(
            span=self._span_from(start),
            name=name,
            uses=uses,
            native_blocks=native_blocks,
            native_functions=native_functions,
            functions=functions
        )


def parse(tokens: List[Token], filename: Optional[str] = None, source: Optional[str] = None) -> Module:
    """
    Convenience function to parse tokens into a module.

    Args:
        tokens: List of tokens from the lexer
        filename: Optional filename for error messages
        source: Optional original source code for extracting raw text

    Returns:
        Parsed Module AST

    Raises:
        ParserError: If parsing fails
    """
    parser = Parser(tokens, filename, source)
    return parser.parse_module()
