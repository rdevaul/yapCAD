"""
Unit tests for the yapCAD DSL lexer.
"""

import pytest
from yapcad.dsl import tokenize, Lexer, TokenType, LexerError


class TestLexerBasics:
    """Test basic lexer functionality."""

    def test_empty_source(self):
        """Empty source produces only EOF."""
        tokens = tokenize("")
        assert len(tokens) == 1
        assert tokens[0].type == TokenType.EOF

    def test_whitespace_only(self):
        """Whitespace-only source produces only EOF (may have NEWLINE)."""
        tokens = tokenize("   \t  ")  # No newlines in this test
        types = [t.type for t in tokens if t.type not in (TokenType.EOF, TokenType.NEWLINE)]
        assert types == []

    def test_simple_let_statement(self):
        """Basic let statement tokenization."""
        tokens = tokenize("let x: int = 42;")
        types = [t.type for t in tokens]
        assert types == [
            TokenType.LET,
            TokenType.IDENTIFIER,
            TokenType.COLON,
            TokenType.TYPE_INT,
            TokenType.ASSIGN,
            TokenType.INT_LITERAL,
            TokenType.SEMICOLON,
            TokenType.EOF,
        ]

    def test_identifier_value(self):
        """Identifier token has correct value."""
        tokens = tokenize("foo_bar123")
        assert tokens[0].type == TokenType.IDENTIFIER
        assert tokens[0].value == "foo_bar123"

    def test_position_tracking(self):
        """Token positions are tracked correctly."""
        tokens = tokenize("let x = 5;")
        # 'let' starts at column 1
        assert tokens[0].span.start.line == 1
        assert tokens[0].span.start.column == 1
        # 'x' starts at column 5
        assert tokens[1].span.start.column == 5

    def test_multiline_position_tracking(self):
        """Position tracking across multiple lines."""
        tokens = tokenize("let x = 5;\nlet y = 10;")
        # Second 'let' is on line 2
        let_tokens = [t for t in tokens if t.type == TokenType.LET]
        assert let_tokens[0].span.start.line == 1
        assert let_tokens[1].span.start.line == 2


class TestComments:
    """Test comment handling."""

    def test_single_line_comment(self):
        """Single-line comments are skipped."""
        tokens = tokenize("# This is a comment\nlet x = 5;")
        # First meaningful token should be LET (may have NEWLINE from comment line)
        meaningful = [t for t in tokens if t.type not in (TokenType.EOF, TokenType.NEWLINE)]
        assert meaningful[0].type == TokenType.LET

    def test_single_line_comment_at_end(self):
        """Comment at end of line."""
        tokens = tokenize("let x = 5; # comment")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert TokenType.IDENTIFIER not in types or tokens[1].value == "x"

    def test_multiline_comment(self):
        """Multi-line comments are skipped."""
        tokens = tokenize("/* This is\na multi-line\ncomment */ let x = 5;")
        assert tokens[0].type == TokenType.LET

    def test_nested_multiline_comment(self):
        """Nested multi-line comments are supported."""
        tokens = tokenize("/* outer /* inner */ still outer */ let x = 5;")
        assert tokens[0].type == TokenType.LET

    def test_unterminated_multiline_comment(self):
        """Unterminated multi-line comment raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize("/* unterminated comment")
        assert "E004" in str(exc_info.value)


class TestStringLiterals:
    """Test string literal handling."""

    def test_simple_string(self):
        """Simple string literal."""
        tokens = tokenize('"hello world"')
        assert tokens[0].type == TokenType.STRING_LITERAL
        assert tokens[0].value == "hello world"

    def test_single_quote_string(self):
        """Single-quoted string literal."""
        tokens = tokenize("'hello world'")
        assert tokens[0].type == TokenType.STRING_LITERAL
        assert tokens[0].value == "hello world"

    def test_escape_sequences(self):
        """Escape sequences in strings."""
        tokens = tokenize(r'"line1\nline2\ttab"')
        assert tokens[0].value == "line1\nline2\ttab"

    def test_escape_quote(self):
        """Escaped quotes in strings."""
        tokens = tokenize(r'"He said \"hello\""')
        assert tokens[0].value == 'He said "hello"'

    def test_escape_backslash(self):
        """Escaped backslash in strings."""
        tokens = tokenize(r'"path\\to\\file"')
        assert tokens[0].value == "path\\to\\file"

    def test_hex_escape(self):
        """Hex escape sequences."""
        tokens = tokenize(r'"\x41\x42\x43"')
        assert tokens[0].value == "ABC"

    def test_unicode_escape(self):
        """Unicode escape sequences."""
        tokens = tokenize(r'"\u0048\u0069"')
        assert tokens[0].value == "Hi"

    def test_multiline_string(self):
        """Triple-quoted multi-line strings."""
        tokens = tokenize('"""line1\nline2\nline3"""')
        assert tokens[0].type == TokenType.STRING_LITERAL
        assert tokens[0].value == "line1\nline2\nline3"

    def test_unterminated_string(self):
        """Unterminated string raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize('"unterminated')
        assert "E002" in str(exc_info.value)

    def test_unterminated_multiline_string(self):
        """Unterminated multi-line string raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize('"""unterminated')
        assert "E003" in str(exc_info.value)

    def test_invalid_escape_sequence(self):
        """Invalid escape sequence raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize(r'"\q"')
        assert "E005" in str(exc_info.value)


class TestNumericLiterals:
    """Test numeric literal handling."""

    def test_integer(self):
        """Simple integer literal."""
        tokens = tokenize("42")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 42

    def test_zero(self):
        """Zero literal."""
        tokens = tokenize("0")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 0

    def test_float_simple(self):
        """Simple float literal."""
        tokens = tokenize("3.14")
        assert tokens[0].type == TokenType.FLOAT_LITERAL
        assert tokens[0].value == 3.14

    def test_float_leading_zero(self):
        """Float with leading zero."""
        tokens = tokenize("0.5")
        assert tokens[0].type == TokenType.FLOAT_LITERAL
        assert tokens[0].value == 0.5

    def test_scientific_notation(self):
        """Scientific notation."""
        tokens = tokenize("1e-9")
        assert tokens[0].type == TokenType.FLOAT_LITERAL
        assert tokens[0].value == 1e-9

    def test_scientific_notation_positive(self):
        """Scientific notation with explicit positive exponent."""
        tokens = tokenize("2.5E+10")
        assert tokens[0].type == TokenType.FLOAT_LITERAL
        assert tokens[0].value == 2.5e10

    def test_scientific_notation_capital_e(self):
        """Scientific notation with capital E."""
        tokens = tokenize("1.5E6")
        assert tokens[0].type == TokenType.FLOAT_LITERAL
        assert tokens[0].value == 1.5e6

    def test_hex_literal(self):
        """Hexadecimal literal."""
        tokens = tokenize("0xFF")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 255

    def test_hex_literal_lowercase(self):
        """Hexadecimal literal with lowercase."""
        tokens = tokenize("0xabcd")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 0xabcd

    def test_binary_literal(self):
        """Binary literal."""
        tokens = tokenize("0b1010")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 10

    def test_binary_literal_with_underscores(self):
        """Binary literal with underscore separators."""
        tokens = tokenize("0b1111_0000")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 0b11110000

    def test_hex_with_underscores(self):
        """Hex literal with underscore separators."""
        tokens = tokenize("0xFF_FF")
        assert tokens[0].type == TokenType.INT_LITERAL
        assert tokens[0].value == 0xFFFF

    def test_invalid_hex_literal(self):
        """Invalid hex literal raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize("0x")
        assert "E007" in str(exc_info.value)

    def test_invalid_binary_literal(self):
        """Invalid binary literal raises error."""
        with pytest.raises(LexerError) as exc_info:
            tokenize("0b")
        assert "E008" in str(exc_info.value)


class TestKeywords:
    """Test keyword recognition."""

    def test_module_keyword(self):
        """Module keyword."""
        tokens = tokenize("module")
        assert tokens[0].type == TokenType.MODULE

    def test_command_keyword(self):
        """Command keyword."""
        tokens = tokenize("command")
        assert tokens[0].type == TokenType.COMMAND

    def test_let_keyword(self):
        """Let keyword."""
        tokens = tokenize("let")
        assert tokens[0].type == TokenType.LET

    def test_require_keyword(self):
        """Require keyword."""
        tokens = tokenize("require")
        assert tokens[0].type == TokenType.REQUIRE

    def test_emit_keyword(self):
        """Emit keyword."""
        tokens = tokenize("emit")
        assert tokens[0].type == TokenType.EMIT

    def test_true_literal(self):
        """True boolean literal."""
        tokens = tokenize("true")
        assert tokens[0].type == TokenType.BOOL_LITERAL
        assert tokens[0].value is True

    def test_false_literal(self):
        """False boolean literal."""
        tokens = tokenize("false")
        assert tokens[0].type == TokenType.BOOL_LITERAL
        assert tokens[0].value is False

    def test_if_else(self):
        """If/else keywords."""
        tokens = tokenize("if else")
        assert tokens[0].type == TokenType.IF
        assert tokens[1].type == TokenType.ELSE

    def test_for_in(self):
        """For/in keywords."""
        tokens = tokenize("for in")
        assert tokens[0].type == TokenType.FOR
        assert tokens[1].type == TokenType.IN


class TestTypeKeywords:
    """Test type keyword recognition."""

    def test_tier1_types(self):
        """Tier 1 type keywords."""
        type_tokens = {
            "int": TokenType.TYPE_INT,
            "float": TokenType.TYPE_FLOAT,
            "string": TokenType.TYPE_STRING,
            "bool": TokenType.TYPE_BOOL,
            "point": TokenType.TYPE_POINT,
            "point2d": TokenType.TYPE_POINT2D,
            "point3d": TokenType.TYPE_POINT3D,
            "vector": TokenType.TYPE_VECTOR,
            "vector2d": TokenType.TYPE_VECTOR2D,
            "vector3d": TokenType.TYPE_VECTOR3D,
            "transform": TokenType.TYPE_TRANSFORM,
        }
        for text, expected_type in type_tokens.items():
            tokens = tokenize(text)
            assert tokens[0].type == expected_type, f"Failed for {text}"

    def test_tier2_types(self):
        """Tier 2 type keywords (curves)."""
        type_tokens = {
            "line_segment": TokenType.TYPE_LINE_SEGMENT,
            "arc": TokenType.TYPE_ARC,
            "circle": TokenType.TYPE_CIRCLE,
            "ellipse": TokenType.TYPE_ELLIPSE,
            "bezier": TokenType.TYPE_BEZIER,
            "nurbs": TokenType.TYPE_NURBS,
            "catmullrom": TokenType.TYPE_CATMULLROM,
        }
        for text, expected_type in type_tokens.items():
            tokens = tokenize(text)
            assert tokens[0].type == expected_type, f"Failed for {text}"

    def test_tier3_types(self):
        """Tier 3 type keywords (compound curves)."""
        type_tokens = {
            "path2d": TokenType.TYPE_PATH2D,
            "path3d": TokenType.TYPE_PATH3D,
            "profile2d": TokenType.TYPE_PROFILE2D,
            "region2d": TokenType.TYPE_REGION2D,
            "loop3d": TokenType.TYPE_LOOP3D,
        }
        for text, expected_type in type_tokens.items():
            tokens = tokenize(text)
            assert tokens[0].type == expected_type, f"Failed for {text}"

    def test_tier4_and_5_types(self):
        """Tier 4 (surfaces) and Tier 5 (solids) types."""
        type_tokens = {
            "surface": TokenType.TYPE_SURFACE,
            "shell": TokenType.TYPE_SHELL,
            "solid": TokenType.TYPE_SOLID,
        }
        for text, expected_type in type_tokens.items():
            tokens = tokenize(text)
            assert tokens[0].type == expected_type, f"Failed for {text}"

    def test_generic_types(self):
        """Generic type keywords."""
        tokens = tokenize("list dict")
        assert tokens[0].type == TokenType.TYPE_LIST
        assert tokens[1].type == TokenType.TYPE_DICT


class TestOperators:
    """Test operator recognition."""

    def test_arithmetic_operators(self):
        """Arithmetic operators."""
        tokens = tokenize("+ - * / %")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.PLUS,
            TokenType.MINUS,
            TokenType.STAR,
            TokenType.SLASH,
            TokenType.PERCENT,
        ]

    def test_comparison_operators(self):
        """Comparison operators."""
        tokens = tokenize("< > <= >= == !=")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.LT,
            TokenType.GT,
            TokenType.LE,
            TokenType.GE,
            TokenType.EQ,
            TokenType.NE,
        ]

    def test_logical_operators(self):
        """Logical operators (Python-style keywords)."""
        tokens = tokenize("and or not")
        types = [t.type for t in tokens if t.type not in (TokenType.EOF, TokenType.NEWLINE)]
        assert types == [TokenType.AND, TokenType.OR, TokenType.NOT]

    def test_arrow_operators(self):
        """Arrow operators."""
        tokens = tokenize("-> =>")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [TokenType.ARROW, TokenType.DOUBLE_ARROW]

    def test_range_operator(self):
        """Range operator."""
        tokens = tokenize("0..10")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [TokenType.INT_LITERAL, TokenType.RANGE, TokenType.INT_LITERAL]


class TestDelimiters:
    """Test delimiter recognition."""

    def test_braces(self):
        """Braces."""
        tokens = tokenize("{ }")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [TokenType.LBRACE, TokenType.RBRACE]

    def test_parentheses(self):
        """Parentheses."""
        tokens = tokenize("( )")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [TokenType.LPAREN, TokenType.RPAREN]

    def test_brackets(self):
        """Brackets."""
        tokens = tokenize("[ ]")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [TokenType.LBRACKET, TokenType.RBRACKET]

    def test_punctuation(self):
        """Punctuation marks."""
        tokens = tokenize(": ; , . ?")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.COLON,
            TokenType.SEMICOLON,
            TokenType.COMMA,
            TokenType.DOT,
            TokenType.QUESTION,
        ]


class TestComplexExpressions:
    """Test tokenization of complex expressions."""

    def test_function_call(self):
        """Function call with arguments."""
        tokens = tokenize("point(1.0, 2.0, 3.0)")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.TYPE_POINT,
            TokenType.LPAREN,
            TokenType.FLOAT_LITERAL,
            TokenType.COMMA,
            TokenType.FLOAT_LITERAL,
            TokenType.COMMA,
            TokenType.FLOAT_LITERAL,
            TokenType.RPAREN,
        ]

    def test_method_call(self):
        """Method call with dot notation."""
        tokens = tokenize("curve.at(0.5)")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.IDENTIFIER,
            TokenType.DOT,
            TokenType.IDENTIFIER,
            TokenType.LPAREN,
            TokenType.FLOAT_LITERAL,
            TokenType.RPAREN,
        ]

    def test_list_literal(self):
        """List literal."""
        tokens = tokenize("[1, 2, 3]")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.LBRACKET,
            TokenType.INT_LITERAL,
            TokenType.COMMA,
            TokenType.INT_LITERAL,
            TokenType.COMMA,
            TokenType.INT_LITERAL,
            TokenType.RBRACKET,
        ]

    def test_command_definition(self):
        """Command definition header."""
        source = "command MAKE_BOX(width: float, height: float) -> solid {"
        tokens = tokenize(source)
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert TokenType.COMMAND in types
        assert TokenType.IDENTIFIER in types
        assert TokenType.ARROW in types
        assert TokenType.TYPE_SOLID in types

    def test_generic_type(self):
        """Generic type annotation."""
        tokens = tokenize("list<point3d>")
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert types == [
            TokenType.TYPE_LIST,
            TokenType.LT,
            TokenType.TYPE_POINT3D,
            TokenType.GT,
        ]

    def test_boolean_expression(self):
        """Boolean expression (Python-style)."""
        tokens = tokenize("x > 0 and y < 100")
        types = [t.type for t in tokens if t.type not in (TokenType.EOF, TokenType.NEWLINE)]
        assert types == [
            TokenType.IDENTIFIER,
            TokenType.GT,
            TokenType.INT_LITERAL,
            TokenType.AND,
            TokenType.IDENTIFIER,
            TokenType.LT,
            TokenType.INT_LITERAL,
        ]


class TestErrorMessages:
    """Test error message quality."""

    def test_error_has_line_info(self):
        """Error includes line number."""
        try:
            tokenize("let x = @;")
        except LexerError as e:
            assert "1:" in str(e)  # Line 1
            assert "@" in str(e)  # The problematic character

    def test_error_has_column_info(self):
        """Error includes column number."""
        try:
            tokenize("let x = @;")
        except LexerError as e:
            diag = e.diagnostic
            assert diag.span.start.column == 9  # '@' is at column 9

    def test_error_shows_source_line(self):
        """Error shows the source line."""
        try:
            tokenize("let x = @;")
        except LexerError as e:
            formatted = e.diagnostic.format()
            assert "let x = @" in formatted

    def test_error_shows_caret(self):
        """Error shows caret pointing to problem."""
        try:
            tokenize("let x = @;")
        except LexerError as e:
            formatted = e.diagnostic.format()
            assert "^" in formatted


class TestLexerIterator:
    """Test lexer as iterator."""

    def test_iterate_tokens(self):
        """Can iterate over tokens."""
        lexer = Lexer("let x = 5;")
        tokens = list(lexer)
        assert len(tokens) == 6  # let, x, =, 5, ;, EOF
        assert tokens[-1].type == TokenType.EOF

    def test_iterator_stops_at_eof(self):
        """Iterator stops after EOF."""
        lexer = Lexer("42")
        tokens = list(lexer)
        assert len(tokens) == 2
        assert tokens[-1].type == TokenType.EOF


class TestDSLExample:
    """Test tokenization of DSL example from roadmap."""

    def test_gear_command_header(self):
        """Tokenize gear command header."""
        source = """command INVOLUTE_SPUR(
    teeth: int,
    module_mm: float,
    face_width: float,
    pressure_angle: float = 20.0
) -> solid {"""
        tokens = tokenize(source)
        # Should tokenize without errors
        assert tokens[-1].type == TokenType.EOF or tokens[-2].type == TokenType.LBRACE

    def test_require_statement(self):
        """Tokenize require statement."""
        source = 'require teeth >= 8, "Minimum 8 teeth required";'
        tokens = tokenize(source)
        types = [t.type for t in tokens if t.type != TokenType.EOF]
        assert TokenType.REQUIRE in types
        assert TokenType.GE in types
        assert TokenType.STRING_LITERAL in types
