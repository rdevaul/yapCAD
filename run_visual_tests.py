#!/usr/bin/env python3
"""
Runner for visual/interactive tests that require pyglet windows.
These tests are skipped during normal pytest runs.
Each test runs in a subprocess to prevent pyglet from terminating the runner.
"""

import sys
import os
import re
import glob
import subprocess
import ast

def find_visual_tests():
    """Find all tests marked with @pytest.mark.visual decorator"""
    visual_tests = []

    # Search for test files
    test_files = glob.glob("tests/test_*.py")

    for test_file in test_files:
        module_name = os.path.basename(test_file)[:-3]  # Remove .py

        with open(test_file, 'r') as f:
            content = f.read()

        try:
            tree = ast.parse(content)
        except SyntaxError:
            print(f"Warning: Could not parse {test_file}")
            continue

        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                class_name = node.name
                for item in node.body:
                    if _is_visual_test(item):
                        visual_tests.append((module_name, class_name, item.name))
            elif _is_visual_test(node):
                visual_tests.append((module_name, None, node.name))

    return visual_tests


def _is_visual_test(func_def):
    if not isinstance(func_def, ast.FunctionDef):
        return False
    if not func_def.name.startswith('test'):
        return False
    for decorator in func_def.decorator_list:
        if isinstance(decorator, ast.Attribute):
            if (hasattr(decorator.value, 'attr') and
                    decorator.value.attr == 'mark' and
                    decorator.attr == 'visual'):
                return True
    return False

def find_visual_tests_regex():
    """Alternative method using regex - fallback if AST parsing fails"""
    visual_tests = []
    test_files = glob.glob("tests/test_*.py")

    for test_file in test_files:
        module_name = os.path.basename(test_file)[:-3]

        with open(test_file, 'r') as f:
            content = f.read()

        # Find all occurrences of @pytest.mark.visual followed by def test_
        pattern = r'@pytest\.mark\.visual\s*\n\s*(?:class\s+(\w+)\s*:\s*)?def\s+(test\w+)'
        for class_name, method_name in re.findall(pattern, content):
            cls = class_name if class_name else None
            visual_tests.append((module_name, cls, method_name))

    return visual_tests

def run_test_subprocess(module_name, test_class, test_method):
    """Run a single test in a subprocess"""

    # Create a Python script to run in the subprocess
    if test_class:
        test_ref = f"{test_class}.{test_method}"
        resolve_code = f"""\ntest_cls = getattr(module, \"{test_class}\")\ntest_instance = test_cls()\ntest_func = getattr(test_instance, \"{test_method}\")\n"""
    else:
        test_ref = test_method
        resolve_code = f"""\ntest_func = getattr(module, \"{test_method}\")\n"""

    test_script = f"""
import sys
import os

# Add src to path
sys.path.insert(0, 'src')

# Set VISUALTEST environment variable
os.environ['VISUALTEST'] = 'true'

# Import the test module
import importlib.util
spec = importlib.util.spec_from_file_location("{module_name}", "tests/{module_name}.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

# Resolve test function
{resolve_code}

print("Running {module_name}.{test_ref}")
print("Close the window to continue...")

try:
    test_func()
    print("Test completed successfully")
except Exception as e:
    print(f"Test failed with error: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
"""

    # Run the test in a subprocess
    result = subprocess.run(
        [sys.executable, "-c", test_script],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__))
    )

    return result

def run_test_with_pytest(module_name, test_class, test_method):
    """Alternative: Run test using pytest in subprocess"""

    # Set up environment
    env = os.environ.copy()
    env['VISUALTEST'] = 'true'
    env['PYTHONPATH'] = './src'

    # Build pytest command - try to use activated venv if available
    if test_class:
        test_path = f"tests/{module_name}.py::{test_class}::{test_method}"
    else:
        test_path = f"tests/{module_name}.py::{test_method}"

    # Check if we're in a virtual environment and use it
    python_cmd = sys.executable
    if 'v_312' in python_cmd or 'VIRTUAL_ENV' in env:
        # Already using the right python
        cmd = [python_cmd, "-m", "pytest", test_path, "-xvs", "--tb=short"]
    else:
        # Try to activate the virtual environment
        venv_python = "./v_312/bin/python"
        if os.path.exists(venv_python):
            cmd = [venv_python, "-m", "pytest", test_path, "-xvs", "--tb=short"]
        else:
            cmd = [python_cmd, "-m", "pytest", test_path, "-xvs", "--tb=short"]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env,
        cwd=os.path.dirname(os.path.abspath(__file__))
    )

    return result

def main():
    """Main runner"""

    # Check for command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            print("Usage: python run_visual_tests.py [options] [test_pattern]")
            print("\nOptions:")
            print("  -h, --help     Show this help message")
            print("  --pytest       Use pytest to run tests (default: direct execution)")
            print("\nExamples:")
            print("  python run_visual_tests.py              # Run all visual tests")
            print("  python run_visual_tests.py test_geom    # Run only test_geom visual tests")
            print("  python run_visual_tests.py Face         # Run tests matching 'Face'")
            print("  python run_visual_tests.py --pytest     # Run all tests using pytest")
            sys.exit(0)

        use_pytest = False
        pattern = None

        for arg in sys.argv[1:]:
            if arg == '--pytest':
                use_pytest = True
            else:
                pattern = arg
    else:
        use_pytest = False
        pattern = None

    # Find all visual tests using AST parsing
    visual_tests = find_visual_tests()

    # If AST parsing didn't find any, try regex method
    if not visual_tests:
        print("No tests found via AST parsing, trying regex method...")
        visual_tests = find_visual_tests_regex()

    # Filter if pattern provided
    if pattern:
        filtered = []
        pattern_lower = pattern.lower()
        for module, cls, method in visual_tests:
            cls_name = cls or ''
            if (pattern_lower in module.lower() or
                    pattern_lower in cls_name.lower() or
                    pattern_lower in method.lower()):
                filtered.append((module, cls, method))
        visual_tests = filtered

    if not visual_tests:
        print("No visual tests found!")
        if pattern:
            print(f"(No tests matching pattern: {pattern})")
        else:
            print("Make sure tests are marked with @pytest.mark.visual")
        sys.exit(0)

    print("="*60)
    print("Visual Test Runner (Subprocess Isolated)")
    print("="*60)
    print(f"Found {len(visual_tests)} visual test(s) marked with @pytest.mark.visual")
    if pattern:
        print(f"Filtered by pattern: {pattern}")
    if use_pytest:
        print("Using pytest subprocess runner")
    else:
        print("Using direct execution runner")

    print("\nTests to run:")
    for i, (module, cls, method) in enumerate(visual_tests, 1):
        cls_display = cls if cls is not None else '<module>'
        print(f"  {i}. {module}::{cls_display}::{method}")

    print("\nEach test will open in a separate window.")
    print("Close each window to continue to the next test.")
    print("Press Ctrl+C to abort.\n")

    if sys.stdin.isatty():
        input("Press Enter to start...")
    else:
        print("Non-interactive environment detected; starting immediately.\n")

    passed = 0
    failed = 0
    skipped = 0

    for i, (module, cls, method) in enumerate(visual_tests, 1):
        print(f"\n{'='*60}")
        cls_display = cls if cls is not None else '<module>'
        print(f"Test {i}/{len(visual_tests)}: {module}::{cls_display}::{method}")
        print("="*60)

        try:
            # Run test in subprocess
            if use_pytest:
                result = run_test_with_pytest(module, cls, method)
            else:
                result = run_test_subprocess(module, cls, method)

            # Print output
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print("STDERR:", result.stderr)

            # Check result
            if result.returncode == 0:
                print(f"✓ Test passed")
                passed += 1
            else:
                print(f"✗ Test failed (exit code: {result.returncode})")
                failed += 1

                # Ask whether to continue
                if i < len(visual_tests):
                    response = input("Continue with next test? (y/n): ")
                    if response.lower() != 'y':
                        skipped = len(visual_tests) - i
                        break

        except KeyboardInterrupt:
            print("\n\nTests interrupted by user")
            skipped = len(visual_tests) - i
            break
        except Exception as e:
            print(f"Error running test: {e}")
            failed += 1
            if i < len(visual_tests):
                response = input("Continue with next test? (y/n): ")
                if response.lower() != 'y':
                    skipped = len(visual_tests) - i
                    break

    # Print summary
    print("\n" + "="*60)
    print("Test Run Summary")
    print("="*60)
    print(f"Passed:  {passed}")
    print(f"Failed:  {failed}")
    if skipped:
        print(f"Skipped: {skipped}")
    print(f"Total:   {len(visual_tests)}")
    print("="*60)

if __name__ == "__main__":
    main()
