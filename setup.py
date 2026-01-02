#!/usr/bin/env python3
"""
Setup script for installing hf_tdoa as a Python package.

This allows the hf_tdoa library to be installed in development mode,
making it accessible from anywhere on the system.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file for the long description
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "HF TDOA Analysis Library for Ionospheric Studies"

setup(
    name="hf_tdoa",
    version="1.0.1",
    description="HF TDOA Analysis Library for Ionospheric Studies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="HF TDOA Research Team",
    packages=find_packages(),
    python_requires=">=3.12",
    install_requires=[
        "numpy<2",
        "astropy>=6,<7",
        "matplotlib",
        "scipy",
        "pandas",
        "cartopy",
        "geographiclib",
        "pydub",
        "tqdm",
    ],
    extras_require={
        "test": [
            "pillow",
            "scikit-image",
            "nbconvert",
        ],
        "jupyter": [
            "ipykernel",
            "jupyterlab",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
    ],
    include_package_data=True,
    zip_safe=False,
)
