"""
Lexer for the yapCAD DSL v2 (Pythonic Syntax).

Converts source text into a stream of tokens for the parser.
Supports:
- Python-style indentation (INDENT/DEDENT tokens)
- Significant newlines (NEWLINE tokens)
- Implicit line continuation inside brackets
- Single-line comments (#)
- Multi-line comments (/* */)
- String literals with escape sequences
- Multi-line strings (triple quotes)
- Integer literals (decimal, hex, binary)
- Float literals (including scientific notation)
- All DSL keywords and operators
"""

from typing import List, Optional, Iterator
from .tokens import (
    Token, TokenType, SourceLocation, SourceSpan, KEYWORDS,
    is_deprecated_keyword, get_deprecation_message
)
from .errors import (
    LexerError,
    error_unexpected_character,
    error_unterminated_string,
    error_unterminated_multiline_string,
    error_unterminated_comment,
    error_invalid_escape_sequence,
    error_invalid_number_literal,
    error_invalid_hex_literal,
    error_invalid_binary_literal,
)


class Lexer:
    """
    Tokenizer for the yapCAD DSL v2 with Python-style indentation.

    This lexer generates INDENT and DEDENT tokens based on changes in
    leading whitespace, similar to Python's tokenizer.

    Usage:
        lexer = Lexer(source_code)
        tokens = lexer.tokenize()

    Or for streaming:
        lexer = Lexer(source_code)
        for token in lexer:
            process(token)
    """

    def __init__(self, source: str, filename: Optional[str] = None):
        self.source = source
        self.filename = filename
        self.pos = 0            # Current position in source
        self.line = 1           # Current line (1-indexed)
        self.column = 1         # Current column (1-indexed)
        self.line_start = 0     # Position of current line start
        self._lines: Optional[List[str]] = None  # Cached line list

        # Indentation tracking
        self.indent_stack = [0]  # Stack of indentation levels (starts at 0)
        self.at_line_start = True  # Are we at the start of a line?
        self.pending_tokens: List[Token] = []  # Tokens to emit before next scan

        # Bracket nesting for implicit line continuation
        self.bracket_depth = 0  # Count of (, [, { nesting

    @property
    def lines(self) -> List[str]:
        """Lazy-load line list for error reporting."""
        if self._lines is None:
            self._lines = self.source.splitlines()
        return self._lines

    def get_source_line(self, line_num: int) -> Optional[str]:
        """Get a specific line of source (1-indexed)."""
        if 1 <= line_num <= len(self.lines):
            return self.lines[line_num - 1]
        return None

    def _location(self) -> SourceLocation:
        """Get current source location."""
        return SourceLocation(self.line, self.column, self.pos, self.filename)

    def _span(self, start: SourceLocation) -> SourceSpan:
        """Create a span from start to current position."""
        return SourceSpan(start, self._location())

    def _peek(self, offset: int = 0) -> str:
        """Look at character at current position + offset without consuming."""
        idx = self.pos + offset
        if idx >= len(self.source):
            return '\0'
        return self.source[idx]

    def _advance(self) -> str:
        """Consume and return current character."""
        if self.pos >= len(self.source):
            return '\0'
        ch = self.source[self.pos]
        self.pos += 1
        if ch == '\n':
            self.line += 1
            self.column = 1
            self.line_start = self.pos
            self.at_line_start = True
        else:
            self.column += 1
        return ch

    def _match(self, expected: str) -> bool:
        """Consume character if it matches expected."""
        if self._peek() == expected:
            self._advance()
            return True
        return False

    def _is_at_end(self) -> bool:
        """Check if we've reached end of source."""
        return self.pos >= len(self.source)

    def _skip_comment(self) -> None:
        """Skip a single-line comment (# to end of line)."""
        while self._peek() != '\n' and not self._is_at_end():
            self._advance()

    def _skip_multiline_comment(self) -> None:
        """Skip /* ... */ comment."""
        start = self._location()
        self._advance()  # consume '/'
        self._advance()  # consume '*'
        depth = 1  # Support nested comments

        while not self._is_at_end() and depth > 0:
            if self._peek() == '/' and self._peek(1) == '*':
                self._advance()
                self._advance()
                depth += 1
            elif self._peek() == '*' and self._peek(1) == '/':
                self._advance()
                self._advance()
                depth -= 1
            else:
                self._advance()

        if depth > 0:
            raise error_unterminated_comment(
                self._span(start),
                self.get_source_line(start.line)
            )

    def _handle_line_start(self) -> Optional[Token]:
        """
        Handle indentation at the start of a line.
        Returns a NEWLINE token if we have a significant newline,
        and queues INDENT/DEDENT tokens as needed.
        """
        if not self.at_line_start:
            return None

        # Calculate indentation (spaces and tabs)
        # Note: We treat tabs as single indent units, but warn about mixing
        start = self._location()
        indent = 0
        while self._peek() in ' \t':
            if self._peek() == ' ':
                indent += 1
            else:  # tab
                indent += 1  # Simple approach: tab = 1 indent unit
            self._advance()

        # Skip blank lines and comment-only lines
        if self._peek() == '\n' or self._peek() == '\0':
            self.at_line_start = True
            return None
        if self._peek() == '#':
            self._skip_comment()
            return None
        if self._peek() == '/' and self._peek(1) == '*':
            self._skip_multiline_comment()
            # After multiline comment, check if we're at line start again
            return None

        self.at_line_start = False

        # Inside brackets, indentation doesn't matter
        if self.bracket_depth > 0:
            return None

        current_indent = self.indent_stack[-1]

        if indent > current_indent:
            # Indent
            self.indent_stack.append(indent)
            return self._make_token(TokenType.INDENT, None, start, "")
        elif indent < current_indent:
            # Dedent - may need multiple DEDENT tokens
            dedent_tokens = []
            while self.indent_stack[-1] > indent:
                self.indent_stack.pop()
                dedent_tokens.append(
                    self._make_token(TokenType.DEDENT, None, start, "")
                )
            # Check for inconsistent dedent
            if self.indent_stack[-1] != indent:
                # Dedent to a level that was never indented to
                # This is an error in Python, but we'll be lenient
                self.indent_stack.append(indent)

            if dedent_tokens:
                # Queue all but the first
                self.pending_tokens.extend(dedent_tokens[1:])
                return dedent_tokens[0]

        return None

    def _skip_whitespace_within_line(self) -> None:
        """Skip horizontal whitespace (spaces and tabs), not newlines."""
        while self._peek() in ' \t':
            self._advance()

    def _make_token(self, token_type: TokenType, value, start: SourceLocation,
                    lexeme: Optional[str] = None) -> Token:
        """Create a token."""
        span = self._span(start)
        if lexeme is None:
            lexeme = self.source[start.offset:self.pos]
        return Token(token_type, value, lexeme, span)

    def _scan_string(self) -> Token:
        """Scan a string literal."""
        start = self._location()
        quote = self._advance()  # consume opening quote

        # Check for triple-quote (multi-line string)
        if self._peek() == quote and self._peek(1) == quote:
            self._advance()  # consume second quote
            self._advance()  # consume third quote
            return self._scan_multiline_string(start, quote)

        # Single-line string
        chars = []
        while not self._is_at_end() and self._peek() != quote:
            ch = self._peek()
            if ch == '\n':
                raise error_unterminated_string(
                    self._span(start),
                    self.get_source_line(start.line)
                )
            if ch == '\\':
                self._advance()  # consume backslash
                chars.append(self._scan_escape_sequence(start))
            else:
                chars.append(self._advance())

        if self._is_at_end():
            raise error_unterminated_string(
                self._span(start),
                self.get_source_line(start.line)
            )

        self._advance()  # consume closing quote
        value = ''.join(chars)
        return self._make_token(TokenType.STRING_LITERAL, value, start)

    def _scan_multiline_string(self, start: SourceLocation, quote: str) -> Token:
        """Scan a triple-quoted multi-line string."""
        chars = []
        while not self._is_at_end():
            if self._peek() == quote and self._peek(1) == quote and self._peek(2) == quote:
                self._advance()  # consume quotes
                self._advance()
                self._advance()
                value = ''.join(chars)
                return self._make_token(TokenType.STRING_LITERAL, value, start)

            ch = self._peek()
            if ch == '\\':
                self._advance()
                chars.append(self._scan_escape_sequence(start))
            else:
                chars.append(self._advance())

        raise error_unterminated_multiline_string(
            self._span(start),
            self.get_source_line(start.line)
        )

    def _scan_escape_sequence(self, string_start: SourceLocation) -> str:
        """Parse an escape sequence after backslash."""
        esc_start = self._location()
        if self._is_at_end():
            raise error_invalid_escape_sequence(
                "", self._span(esc_start), self.get_source_line(esc_start.line)
            )

        ch = self._advance()
        escape_chars = {
            'n': '\n',
            't': '\t',
            'r': '\r',
            '\\': '\\',
            '"': '"',
            "'": "'",
            '0': '\0',
        }

        if ch in escape_chars:
            return escape_chars[ch]
        elif ch == 'x':
            # Hex escape: \xHH
            hex_chars = self._advance() + self._advance()
            try:
                return chr(int(hex_chars, 16))
            except ValueError:
                raise error_invalid_escape_sequence(
                    f"x{hex_chars}", self._span(esc_start),
                    self.get_source_line(esc_start.line)
                )
        elif ch == 'u':
            # Unicode escape: \uHHHH
            hex_chars = ''
            for _ in range(4):
                hex_chars += self._advance()
            try:
                return chr(int(hex_chars, 16))
            except ValueError:
                raise error_invalid_escape_sequence(
                    f"u{hex_chars}", self._span(esc_start),
                    self.get_source_line(esc_start.line)
                )
        else:
            raise error_invalid_escape_sequence(
                ch, self._span(esc_start), self.get_source_line(esc_start.line)
            )

    def _scan_number(self) -> Token:
        """Scan a numeric literal (int or float)."""
        start = self._location()
        first_char = self._peek()

        # Check for hex or binary prefix
        if first_char == '0':
            self._advance()
            if self._peek() in 'xX':
                return self._scan_hex_number(start)
            elif self._peek() in 'bB':
                return self._scan_binary_number(start)
            # Otherwise it's a decimal starting with 0

        # Scan decimal digits
        while self._peek().isdigit():
            self._advance()

        # Check for float (decimal point)
        is_float = False
        if self._peek() == '.' and self._peek(1).isdigit():
            is_float = True
            self._advance()  # consume '.'
            while self._peek().isdigit():
                self._advance()

        # Check for scientific notation
        if self._peek() in 'eE':
            is_float = True
            self._advance()  # consume 'e'
            if self._peek() in '+-':
                self._advance()
            if not self._peek().isdigit():
                lexeme = self.source[start.offset:self.pos]
                raise error_invalid_number_literal(
                    lexeme, self._span(start), self.get_source_line(start.line)
                )
            while self._peek().isdigit():
                self._advance()

        lexeme = self.source[start.offset:self.pos]
        try:
            if is_float:
                value = float(lexeme)
                return self._make_token(TokenType.FLOAT_LITERAL, value, start, lexeme)
            else:
                value = int(lexeme)
                return self._make_token(TokenType.INT_LITERAL, value, start, lexeme)
        except ValueError:
            raise error_invalid_number_literal(
                lexeme, self._span(start), self.get_source_line(start.line)
            )

    def _scan_hex_number(self, start: SourceLocation) -> Token:
        """Scan a hexadecimal integer literal (0x...)."""
        self._advance()  # consume 'x' or 'X'

        if not self._peek() in '0123456789abcdefABCDEF':
            lexeme = self.source[start.offset:self.pos]
            raise error_invalid_hex_literal(
                lexeme, self._span(start), self.get_source_line(start.line)
            )

        while self._peek() in '0123456789abcdefABCDEF_':
            if self._peek() != '_':  # Allow underscores as separators
                pass
            self._advance()

        lexeme = self.source[start.offset:self.pos]
        clean_lexeme = lexeme.replace('_', '')
        try:
            value = int(clean_lexeme, 16)
            return self._make_token(TokenType.INT_LITERAL, value, start, lexeme)
        except ValueError:
            raise error_invalid_hex_literal(
                lexeme, self._span(start), self.get_source_line(start.line)
            )

    def _scan_binary_number(self, start: SourceLocation) -> Token:
        """Scan a binary integer literal (0b...)."""
        self._advance()  # consume 'b' or 'B'

        if not self._peek() in '01':
            lexeme = self.source[start.offset:self.pos]
            raise error_invalid_binary_literal(
                lexeme, self._span(start), self.get_source_line(start.line)
            )

        while self._peek() in '01_':
            self._advance()

        lexeme = self.source[start.offset:self.pos]
        clean_lexeme = lexeme.replace('_', '')
        try:
            value = int(clean_lexeme, 2)
            return self._make_token(TokenType.INT_LITERAL, value, start, lexeme)
        except ValueError:
            raise error_invalid_binary_literal(
                lexeme, self._span(start), self.get_source_line(start.line)
            )

    def _scan_identifier_or_keyword(self) -> Token:
        """Scan an identifier or keyword."""
        start = self._location()

        while self._peek().isalnum() or self._peek() == '_':
            self._advance()

        lexeme = self.source[start.offset:self.pos]

        # Check if it's a keyword
        if lexeme in KEYWORDS:
            token_type = KEYWORDS[lexeme]
            # Handle boolean literals specially (True/False and true/false)
            if lexeme in ('True', 'False', 'true', 'false'):
                value = lexeme in ('True', 'true')
            else:
                value = lexeme
            return self._make_token(token_type, value, start, lexeme)

        # It's an identifier
        return self._make_token(TokenType.IDENTIFIER, lexeme, start, lexeme)

    def _scan_token(self) -> Token:
        """Scan the next token."""
        # Return any pending tokens first (from DEDENT sequences)
        if self.pending_tokens:
            return self.pending_tokens.pop(0)

        # Handle line start indentation
        indent_token = self._handle_line_start()
        if indent_token:
            return indent_token

        # Skip horizontal whitespace
        self._skip_whitespace_within_line()

        # Check for comments (but not at line start, handled above)
        while self._peek() == '#' or (self._peek() == '/' and self._peek(1) == '*'):
            if self._peek() == '#':
                self._skip_comment()
            else:
                self._skip_multiline_comment()
            self._skip_whitespace_within_line()

        # Check for newline
        if self._peek() == '\n':
            start = self._location()
            self._advance()
            # Only emit NEWLINE if not inside brackets
            if self.bracket_depth == 0:
                # Handle indentation for next line
                return self._make_token(TokenType.NEWLINE, None, start, "\\n")
            else:
                # Inside brackets - skip newline, handle indentation
                return self._scan_token()

        if self._is_at_end():
            # At EOF, emit any remaining DEDENT tokens
            start = self._location()
            while len(self.indent_stack) > 1:
                self.indent_stack.pop()
                self.pending_tokens.append(
                    self._make_token(TokenType.DEDENT, None, start, "")
                )
            if self.pending_tokens:
                return self.pending_tokens.pop(0)
            return self._make_token(TokenType.EOF, None, self._location(), "")

        start = self._location()
        ch = self._peek()

        # String literals
        if ch in '"\'':
            return self._scan_string()

        # Numbers
        if ch.isdigit():
            return self._scan_number()

        # Identifiers and keywords
        if ch.isalpha() or ch == '_':
            return self._scan_identifier_or_keyword()

        # Single-character tokens first (advance after match)
        self._advance()

        # Two-character operators
        if ch == '-' and self._match('>'):
            return self._make_token(TokenType.ARROW, "->", start)
        if ch == '=' and self._match('>'):
            return self._make_token(TokenType.DOUBLE_ARROW, "=>", start)
        if ch == '=' and self._match('='):
            return self._make_token(TokenType.EQ, "==", start)
        if ch == '!' and self._match('='):
            return self._make_token(TokenType.NE, "!=", start)
        if ch == '<' and self._match('='):
            return self._make_token(TokenType.LE, "<=", start)
        if ch == '>' and self._match('='):
            return self._make_token(TokenType.GE, ">=", start)
        if ch == '.' and self._match('.'):
            return self._make_token(TokenType.RANGE, "..", start)
        if ch == '*' and self._match('*'):
            return self._make_token(TokenType.DOUBLE_STAR, "**", start)
        if ch == '/' and self._match('/'):
            return self._make_token(TokenType.DOUBLE_SLASH, "//", start)
        if ch == '+' and self._match('='):
            return self._make_token(TokenType.PLUS_ASSIGN, "+=", start)
        if ch == '-' and self._match('='):
            return self._make_token(TokenType.MINUS_ASSIGN, "-=", start)

        # Track bracket depth for implicit line continuation
        if ch == '(':
            self.bracket_depth += 1
            return self._make_token(TokenType.LPAREN, ch, start)
        if ch == ')':
            self.bracket_depth = max(0, self.bracket_depth - 1)
            return self._make_token(TokenType.RPAREN, ch, start)
        if ch == '[':
            self.bracket_depth += 1
            return self._make_token(TokenType.LBRACKET, ch, start)
        if ch == ']':
            self.bracket_depth = max(0, self.bracket_depth - 1)
            return self._make_token(TokenType.RBRACKET, ch, start)
        if ch == '{':
            self.bracket_depth += 1
            return self._make_token(TokenType.LBRACE, ch, start)
        if ch == '}':
            self.bracket_depth = max(0, self.bracket_depth - 1)
            return self._make_token(TokenType.RBRACE, ch, start)

        # Single-character tokens
        single_char_tokens = {
            '+': TokenType.PLUS,
            '-': TokenType.MINUS,
            '*': TokenType.STAR,
            '/': TokenType.SLASH,
            '%': TokenType.PERCENT,
            '<': TokenType.LT,
            '>': TokenType.GT,
            '=': TokenType.ASSIGN,
            ':': TokenType.COLON,
            ';': TokenType.SEMICOLON,  # Deprecated but kept for error messages
            ',': TokenType.COMMA,
            '.': TokenType.DOT,
            '?': TokenType.QUESTION,
            '@': TokenType.AT,
            '#': TokenType.HASH,  # Usually not reached (comments)
        }

        if ch in single_char_tokens:
            return self._make_token(single_char_tokens[ch], ch, start)

        # Unknown character
        raise error_unexpected_character(
            ch, self._span(start), self.get_source_line(start.line)
        )

    def tokenize(self) -> List[Token]:
        """Tokenize the entire source, returning a list of tokens."""
        tokens = []
        while True:
            token = self._scan_token()
            tokens.append(token)
            if token.type == TokenType.EOF:
                break
        return tokens

    def __iter__(self) -> Iterator[Token]:
        """Iterate over tokens."""
        while True:
            token = self._scan_token()
            yield token
            if token.type == TokenType.EOF:
                break


def tokenize(source: str, filename: Optional[str] = None) -> List[Token]:
    """
    Convenience function to tokenize source code.

    Args:
        source: The source code to tokenize
        filename: Optional filename for error messages

    Returns:
        List of tokens

    Raises:
        LexerError: If tokenization fails
    """
    lexer = Lexer(source, filename)
    return lexer.tokenize()
