#!/usr/bin/env python3
"""
Quick verification script to check that the environment is properly set up.
"""

import sys

def check_import(module_name, package_name=None):
    """Check if a module can be imported."""
    try:
        __import__(module_name)
        print(f"  ✓ {package_name or module_name}")
        return True
    except ImportError:
        print(f"  ✗ {package_name or module_name} - MISSING")
        return False

def main():
    print("=" * 70)
    print("HF TDOA Environment Verification")
    print("=" * 70)
    print()

    # Check Python version
    print("Python Version:")
    version = sys.version_info
    print(f"  Python {version.major}.{version.minor}.{version.micro}")
    if version.major == 3 and version.minor >= 12:
        print("  ✓ Version requirement met (>=3.12)")
    else:
        print("  ✗ Version requirement NOT met (need >=3.12)")
        return False
    print()

    # Check core dependencies
    print("Core Dependencies:")
    all_ok = True
    all_ok &= check_import("numpy")
    all_ok &= check_import("scipy")
    all_ok &= check_import("pandas")
    all_ok &= check_import("matplotlib")
    all_ok &= check_import("astropy")
    all_ok &= check_import("cartopy")
    all_ok &= check_import("geographiclib")
    all_ok &= check_import("pydub")
    all_ok &= check_import("tqdm")
    print()

    # Check Jupyter dependencies
    print("Jupyter Dependencies:")
    all_ok &= check_import("IPython", "ipykernel")
    all_ok &= check_import("jupyterlab")
    all_ok &= check_import("nbconvert")
    print()

    # Check testing dependencies
    print("Testing Dependencies:")
    all_ok &= check_import("PIL", "pillow")
    all_ok &= check_import("skimage", "scikit-image")
    print()

    # Check hf_tdoa package
    print("HF TDOA Package:")
    try:
        import hf_tdoa as tdoa
        print(f"  ✓ hf_tdoa")
        print(f"    Location: {tdoa.__file__}")

        # Check key components
        if hasattr(tdoa, 'PathInfo'):
            print("  ✓ PathInfo class available")
        else:
            print("  ✗ PathInfo class NOT available")
            all_ok = False

        if hasattr(tdoa, 'find_chirps'):
            print("  ✓ find_chirps function available")
        else:
            print("  ✗ find_chirps function NOT available")
            all_ok = False

    except ImportError as e:
        print(f"  ✗ hf_tdoa - NOT INSTALLED")
        print(f"    Error: {e}")
        print()
        print("To install, run:")
        print("  pip install -e .")
        all_ok = False

    print()
    print("=" * 70)
    if all_ok:
        print("✓ Environment is properly configured!")
        print()
        print("You can now:")
        print("  - Run notebooks: jupyter lab")
        print("  - Run tests: python test_figures.py")
        return True
    else:
        print("✗ Environment setup is incomplete")
        print()
        print("Please run:")
        print("  conda env create -f environment.yml")
        print("  conda activate hf-tdoa")
        print("  pip install -e .")
        return False
    print("=" * 70)

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
