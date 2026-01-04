;;; yapcad-dsl-mode.el --- Major mode for yapCAD DSL files -*- lexical-binding: t; -*-

;; Copyright (C) 2024 yapCAD Authors
;; Author: yapCAD Project
;; Version: 1.0.0
;; Keywords: languages, cad
;; URL: https://github.com/yapCAD/yapCAD
;; Package-Requires: ((emacs "25.1"))

;; This file is not part of GNU Emacs.

;;; Commentary:

;; Major mode for editing yapCAD DSL (.dsl) files.
;;
;; yapCAD DSL is a domain-specific language for parametric CAD design
;; with Python-like syntax, strong typing, and provenance tracking.
;;
;; Features:
;; - Syntax highlighting for keywords, types, operators, and literals
;; - Python-style indentation (4 spaces)
;; - Comment support (# for line comments)
;; - Automatic pairing of brackets and quotes
;;
;; Installation:
;;   Add to your init.el:
;;     (add-to-list 'load-path "/path/to/yapCAD/editors")
;;     (require 'yapcad-dsl-mode)
;;
;;   Or use use-package:
;;     (use-package yapcad-dsl-mode
;;       :load-path "/path/to/yapCAD/editors"
;;       :mode "\\.dsl\\'")

;;; Code:

(require 'rx)

(defgroup yapcad-dsl nil
  "Major mode for yapCAD DSL files."
  :group 'languages
  :prefix "yapcad-dsl-")

(defcustom yapcad-dsl-indent-offset 4
  "Indentation offset for yapCAD DSL."
  :type 'integer
  :group 'yapcad-dsl)

;;; Syntax Table

(defvar yapcad-dsl-mode-syntax-table
  (let ((table (make-syntax-table)))
    ;; Comments: # to end of line
    (modify-syntax-entry ?# "<" table)
    (modify-syntax-entry ?\n ">" table)

    ;; Strings
    (modify-syntax-entry ?\" "\"" table)

    ;; Operators and punctuation
    (modify-syntax-entry ?+ "." table)
    (modify-syntax-entry ?- "." table)
    (modify-syntax-entry ?* "." table)
    (modify-syntax-entry ?/ "." table)
    (modify-syntax-entry ?% "." table)
    (modify-syntax-entry ?< "." table)
    (modify-syntax-entry ?> "." table)
    (modify-syntax-entry ?= "." table)
    (modify-syntax-entry ?! "." table)

    ;; Underscore is part of identifiers
    (modify-syntax-entry ?_ "w" table)

    ;; Brackets
    (modify-syntax-entry ?\( "()" table)
    (modify-syntax-entry ?\) ")(" table)
    (modify-syntax-entry ?\[ "(]" table)
    (modify-syntax-entry ?\] ")[" table)
    (modify-syntax-entry ?\{ "(}" table)
    (modify-syntax-entry ?\} "){" table)

    table)
  "Syntax table for `yapcad-dsl-mode'.")

;;; Font Lock (Syntax Highlighting)

(defconst yapcad-dsl-keywords
  '(;; Core keywords
    "module" "use" "def" "command" "return" "emit"
    "if" "elif" "else" "for" "in" "while"
    "assert" "require" "pass" "as" "match" "export"
    "native" "let" "with" "fn" "exports" "python"
    ;; Logical operators
    "and" "or" "not"
    ;; Closure
    "close" "closeC0" "closeC1")
  "yapCAD DSL keywords.")

(defconst yapcad-dsl-types
  '(;; Primitive types
    "int" "float" "string" "str" "bool"
    ;; Point/vector types
    "point" "point2d" "point3d"
    "vector" "vector2d" "vector3d"
    "transform"
    ;; Curve types
    "line_segment" "arc" "circle" "ellipse"
    "parabola" "hyperbola" "catmullrom" "nurbs" "bezier"
    ;; Path types
    "path2d" "path3d" "profile2d" "region2d" "loop3d"
    ;; Surface/solid types
    "surface" "shell" "solid"
    ;; Generic types
    "list" "dict")
  "yapCAD DSL type names.")

(defconst yapcad-dsl-constants
  '("True" "False" "true" "false")
  "yapCAD DSL constants.")

(defconst yapcad-dsl-builtins
  '(;; Math functions
    "sin" "cos" "tan" "asin" "acos" "atan" "atan2"
    "sqrt" "pow" "exp" "log" "log10"
    "abs" "floor" "ceil" "round" "min" "max"
    "radians" "degrees" "pi"
    ;; Utility functions
    "range" "len" "print" "type"
    ;; Geometry constructors
    "point" "vector" "line" "arc" "circle" "ellipse"
    "bezier" "catmullrom" "polygon" "rectangle"
    ;; Solid primitives
    "box" "sphere" "cylinder" "cone" "torus"
    "prolate_spheroid" "oblate_spheroid"
    "involute_gear"
    ;; Fasteners
    "metric_hex_bolt" "metric_hex_nut"
    "unified_hex_bolt" "unified_hex_nut"
    ;; 2D operations
    "extrude" "revolve" "sweep" "sweep_hollow"
    "sweep_adaptive" "sweep_adaptive_hollow"
    "sweep_adaptive_frenet" "sweep_adaptive_frenet_hollow"
    "loft" "loft_hollow"
    ;; Boolean operations
    "union" "difference" "intersection" "compound"
    ;; Transforms
    "translate" "rotate" "scale" "mirror"
    "translate_xform" "rotate_xform" "scale_xform"
    ;; Regions/paths
    "region" "path2d" "path3d" "close" "closeC0" "closeC1"
    ;; Query functions
    "volume" "area" "centroid" "bbox")
  "yapCAD DSL builtin functions.")

(defvar yapcad-dsl-font-lock-keywords
  `(
    ;; Comments (handled by syntax table, but ensure highlighting)
    (,(rx line-start (* space) "#" (* any)) . font-lock-comment-face)

    ;; Module declaration
    (,(rx word-start "module" word-end (+ space) (group (+ (any word ?_))))
     (1 font-lock-function-name-face))

    ;; Command/function definitions
    (,(rx word-start (or "def" "command") word-end (+ space)
          (group (+ (any word ?_))))
     (1 font-lock-function-name-face))

    ;; Keywords
    (,(regexp-opt yapcad-dsl-keywords 'words) . font-lock-keyword-face)

    ;; Types
    (,(regexp-opt yapcad-dsl-types 'words) . font-lock-type-face)

    ;; Constants (True/False)
    (,(regexp-opt yapcad-dsl-constants 'words) . font-lock-constant-face)

    ;; Builtins
    (,(regexp-opt yapcad-dsl-builtins 'words) . font-lock-builtin-face)

    ;; Decorators (@native, etc.)
    (,(rx "@" (+ (any word ?_))) . font-lock-preprocessor-face)

    ;; Type annotations (variable: type)
    (,(rx ":" (+ space) (group (+ (any word ?_ ?< ?> ?,))))
     (1 font-lock-type-face))

    ;; Return type (-> type)
    (,(rx "->" (+ space) (group (+ (any word ?_))))
     (1 font-lock-type-face))

    ;; Numeric literals
    (,(rx word-start
          (or
           ;; Hex
           (seq "0" (any "xX") (+ (any "0-9a-fA-F")))
           ;; Binary
           (seq "0" (any "bB") (+ (any "01")))
           ;; Float with exponent
           (seq (? (any "+-")) (+ digit) (? "." (* digit))
                (any "eE") (? (any "+-")) (+ digit))
           ;; Float with decimal
           (seq (? (any "+-")) (+ digit) "." (* digit))
           ;; Integer
           (seq (? (any "+-")) (+ digit)))
          word-end)
     . font-lock-constant-face)

    ;; String literals
    (,(rx "\"" (* (or (not (any "\"\\")) (seq "\\" any))) "\"")
     . font-lock-string-face)

    ;; Variable assignment (name = or name: type =)
    (,(rx word-start (group (+ (any word ?_))) (+ space)
          (? ":" (* (not (any "="))) "="))
     (1 font-lock-variable-name-face))

    ;; Parameter names in function definitions
    (,(rx "(" (* space) (group (+ (any word ?_))) (* space) ":")
     (1 font-lock-variable-name-face)))
  "Font lock keywords for `yapcad-dsl-mode'.")

;;; Indentation

(defun yapcad-dsl-indent-line ()
  "Indent current line as yapCAD DSL code."
  (interactive)
  (let ((indent (yapcad-dsl--calculate-indent)))
    (when indent
      (if (<= (current-column) (current-indentation))
          (indent-line-to indent)
        (save-excursion
          (indent-line-to indent))))))

(defun yapcad-dsl--calculate-indent ()
  "Calculate the indentation for the current line."
  (save-excursion
    (beginning-of-line)
    (if (bobp)
        0
      (let ((not-indented t)
            (cur-indent 0))
        ;; Check if we're at a dedent keyword
        (if (looking-at "^[ \t]*\\(else\\|elif\\)[ \t]*:")
            (progn
              (forward-line -1)
              (while (and (not (bobp))
                          (looking-at "^[ \t]*$"))
                (forward-line -1))
              (setq cur-indent (current-indentation)))
          ;; Normal indentation
          (forward-line -1)
          (while (and (not (bobp))
                      (looking-at "^[ \t]*$"))
            (forward-line -1))
          (setq cur-indent (current-indentation))
          ;; Check if previous line ends with colon
          (end-of-line)
          (if (looking-back ":[ \t]*\\(#.*\\)?$" (line-beginning-position))
              (setq cur-indent (+ cur-indent yapcad-dsl-indent-offset))))
        ;; Check if current line should dedent
        (save-excursion
          (forward-line 1)
          (beginning-of-line)
          (when (looking-at "^[ \t]*\\(return\\|emit\\|pass\\)\\b")
            nil))
        cur-indent))))

;;; Keymap

(defvar yapcad-dsl-mode-map
  (let ((map (make-sparse-keymap)))
    ;; Add any custom keybindings here
    map)
  "Keymap for `yapcad-dsl-mode'.")

;;; Mode Definition

;;;###autoload
(define-derived-mode yapcad-dsl-mode prog-mode "yapCAD-DSL"
  "Major mode for editing yapCAD DSL files.

\\{yapcad-dsl-mode-map}"
  :syntax-table yapcad-dsl-mode-syntax-table
  :group 'yapcad-dsl

  ;; Font lock
  (setq-local font-lock-defaults '(yapcad-dsl-font-lock-keywords))

  ;; Comments
  (setq-local comment-start "# ")
  (setq-local comment-end "")
  (setq-local comment-start-skip "#+ *")

  ;; Indentation
  (setq-local indent-line-function #'yapcad-dsl-indent-line)
  (setq-local indent-tabs-mode nil)
  (setq-local tab-width yapcad-dsl-indent-offset)

  ;; Electric indent
  (setq-local electric-indent-chars (cons ?: electric-indent-chars))

  ;; Paragraph handling
  (setq-local paragraph-start (concat "\\|[ \t]*$\\|" page-delimiter))
  (setq-local paragraph-separate "[ \t]*$\\|[ \t]*#"))

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.dsl\\'" . yapcad-dsl-mode))

(provide 'yapcad-dsl-mode)

;;; yapcad-dsl-mode.el ends here
