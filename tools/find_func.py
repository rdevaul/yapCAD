import ast
import inspect
import sys

def get_function_source(filename, function_name):
    """Prints the source code of a specified function."""
    try:
        with open(filename, "r") as source_file:
            source_code = source_file.read()
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return

    tree = ast.parse(source_code)
    found_func = False

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name == function_name:
                found_func = True

                # Extract the function body using line numbers from AST
                lines = source_code.splitlines()
                start_line = node.lineno - 1
                end_line = node.end_lineno

                # Print the extracted lines
                for line in lines[start_line:end_line]:
                    print(line)
                return

    if not found_func:
        print(f"Function '{function_name}' not found in {filename}.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python find_func.py <filename> <function_name>")
        sys.exit(1)

    filename = sys.argv[1]
    function_name = sys.argv[2]
    get_function_source(filename, function_name)
