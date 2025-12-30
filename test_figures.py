#!/usr/bin/env python3
"""
Test script to verify figure reproducibility.

This script:
1. Executes each Jupyter notebook in the directory
2. Compares newly generated JPG figures to baseline figures in submitted_figures/
3. Reports differences using multiple metrics (pixel differences, structural similarity)
"""

import os
import sys
import subprocess
import tempfile
import shutil
import threading
import time
from pathlib import Path
from typing import Dict, List, Tuple
import json

try:
    import numpy as np
    from PIL import Image
    from skimage.metrics import structural_similarity as ssim
    from tqdm import tqdm
except ImportError as e:
    print(f"Error: Missing required package - {e}")
    print("\nPlease install required packages:")
    print("  pip install numpy pillow scikit-image tqdm")
    sys.exit(1)


class NotebookFigureTester:
    """Test runner for notebook figure generation and comparison."""

    def __init__(self, notebooks_dir: str = ".", submitted_dir: str = "submitted_figures"):
        """
        Initialize the tester.

        Args:
            notebooks_dir: Directory containing notebook subdirectories (default: current directory)
            submitted_dir: Directory containing baseline figures (default: submitted_figures)
        """
        self.notebooks_dir = Path(notebooks_dir).resolve()
        self.submitted_dir = Path(submitted_dir).resolve()

        # Mapping of notebook subdirectories to their notebooks and expected output figures
        # Format: 'subdir': {'notebook': 'name.ipynb', 'figures': ['fig1.jpg', ...]}
        self.notebook_outputs = {
            'fig_03': {'notebook': 'fig_03.ipynb', 'figures': ['fig_03.jpg']},
            'fig_05': {'notebook': 'fig_05.ipynb', 'figures': ['fig_05.jpg']},
            'fig_08': {'notebook': 'fig_08.ipynb', 'figures': ['fig_08.jpg']},
            'fig_09': {'notebook': 'fig_09.ipynb', 'figures': ['fig_09.jpg']},
            'fig_11': {'notebook': 'fig_11.ipynb', 'figures': ['fig_11.jpg']},
            'fig_12': {'notebook': 'fig_12.ipynb', 'figures': ['fig_12.jpg']},
            'fig_13': {'notebook': 'fig_13.ipynb', 'figures': ['fig_13.jpg']},
            'fig_15': {'notebook': 'fig_15.ipynb', 'figures': ['fig_15.jpg']}
        }

    def run_notebook(self, notebook_path: Path) -> Tuple[bool, str]:
        """
        Execute a Jupyter notebook using nbconvert.

        Args:
            notebook_path: Path to the notebook file

        Returns:
            Tuple of (success: bool, error_message: str)
        """
        try:
            print(f"  Executing notebook: {notebook_path.name}")

            # Start progress bar in a separate thread
            pbar = tqdm(total=100, desc="  Progress", bar_format='{desc}: {bar} {elapsed}', leave=False, dynamic_ncols=True)
            stop_event = threading.Event()

            def update_progress():
                """Update progress bar while notebook is running."""
                progress = 0
                while not stop_event.is_set() and progress < 95:
                    time.sleep(1)
                    # Increment by smaller amounts as we get further along
                    if progress < 30:
                        increment = 2
                    elif progress < 60:
                        increment = 1
                    else:
                        increment = 0.5
                    progress += increment
                    pbar.update(increment)

            progress_thread = threading.Thread(target=update_progress, daemon=True)
            progress_thread.start()

            # Execute notebook in place, allowing it to save outputs
            result = subprocess.run(
                [
                    'jupyter', 'nbconvert',
                    '--to', 'notebook',
                    '--execute',
                    '--inplace',
                    '--ExecutePreprocessor.timeout=600',  # 10 minute timeout
                    str(notebook_path)
                ],
                capture_output=True,
                text=True,
                cwd=self.notebooks_dir
            )

            # Stop progress bar
            stop_event.set()
            progress_thread.join(timeout=1)
            pbar.update(100 - pbar.n)  # Complete the bar
            pbar.close()

            if result.returncode != 0:
                error_msg = f"Notebook execution failed:\n{result.stderr}"
                return False, error_msg

            return True, ""

        except FileNotFoundError:
            return False, "jupyter command not found. Please ensure Jupyter is installed."
        except Exception as e:
            return False, f"Unexpected error: {str(e)}"

    def load_image_as_array(self, image_path: Path) -> np.ndarray:
        """
        Load an image file as a numpy array.

        Args:
            image_path: Path to the image file

        Returns:
            Numpy array of the image
        """
        img = Image.open(image_path)
        return np.array(img)

    def compare_images(self, img1_path: Path, img2_path: Path) -> Dict[str, float]:
        """
        Compare two images using multiple metrics.

        Args:
            img1_path: Path to first image (new)
            img2_path: Path to second image (baseline)

        Returns:
            Dictionary containing comparison metrics
        """
        # Load images
        img1 = self.load_image_as_array(img1_path)
        img2 = self.load_image_as_array(img2_path)

        # Check if dimensions match
        if img1.shape != img2.shape:
            return {
                'dimensions_match': False,
                'new_shape': img1.shape,
                'baseline_shape': img2.shape,
                'error': 'Image dimensions do not match'
            }

        # Calculate pixel-wise differences
        diff = np.abs(img1.astype(float) - img2.astype(float))

        # Mean absolute error (MAE)
        mae = np.mean(diff)

        # Root mean squared error (RMSE)
        rmse = np.sqrt(np.mean(diff ** 2))

        # Maximum absolute difference
        max_diff = np.max(diff)

        # Percentage of pixels that differ
        pixels_different = np.sum(diff > 0) / diff.size * 100

        # Structural Similarity Index (SSIM)
        # For color images, calculate SSIM for each channel and average
        if len(img1.shape) == 3:  # Color image
            ssim_value = np.mean([
                ssim(img1[:, :, i], img2[:, :, i], data_range=255)
                for i in range(img1.shape[2])
            ])
        else:  # Grayscale
            ssim_value = ssim(img1, img2, data_range=255)

        return {
            'dimensions_match': True,
            'shape': img1.shape,
            'mae': mae,
            'rmse': rmse,
            'max_diff': max_diff,
            'pixels_different_pct': pixels_different,
            'ssim': ssim_value,
            'identical': mae == 0
        }

    def format_comparison_result(self, figure_name: str, metrics: Dict) -> str:
        """
        Format comparison metrics as a readable string.

        Args:
            figure_name: Name of the figure being compared
            metrics: Dictionary of comparison metrics

        Returns:
            Formatted string report
        """
        lines = [f"\n  {figure_name}:"]

        if 'error' in metrics:
            lines.append(f"    ❌ ERROR: {metrics['error']}")
            if not metrics.get('dimensions_match', True):
                lines.append(f"       New shape: {metrics['new_shape']}")
                lines.append(f"       Baseline shape: {metrics['baseline_shape']}")
            return '\n'.join(lines)

        if metrics['identical']:
            lines.append("    ✓ IDENTICAL - Images match perfectly!")
        else:
            lines.append(f"    Shape: {metrics['shape']}")
            lines.append(f"    Mean Absolute Error: {metrics['mae']:.4f}")
            lines.append(f"    RMSE: {metrics['rmse']:.4f}")
            lines.append(f"    Max Pixel Difference: {metrics['max_diff']:.1f}")
            lines.append(f"    Pixels Different: {metrics['pixels_different_pct']:.2f}%")
            lines.append(f"    SSIM (Structural Similarity): {metrics['ssim']:.6f}")

            # Interpretation
            if metrics['ssim'] > 0.99:
                lines.append("    ✓ EXCELLENT - Negligible differences")
            elif metrics['ssim'] > 0.95:
                lines.append("    ⚠ GOOD - Minor differences detected")
            elif metrics['ssim'] > 0.90:
                lines.append("    ⚠ FAIR - Noticeable differences detected")
            else:
                lines.append("    ❌ POOR - Significant differences detected")

        return '\n'.join(lines)

    def run_tests(self) -> Dict:
        """
        Run all notebook tests and compare generated figures.

        Returns:
            Dictionary containing test results
        """
        results = {
            'notebooks_tested': 0,
            'notebooks_passed': 0,
            'notebooks_failed': 0,
            'figures_compared': 0,
            'figures_identical': 0,
            'details': {}
        }

        print("=" * 70)
        print("Figure Reproducibility Test")
        print("=" * 70)
        print(f"Notebooks directory: {self.notebooks_dir}")
        print(f"Baseline figures: {self.submitted_dir}")
        print()

        # Check that submitted_figures directory exists
        if not self.submitted_dir.exists():
            print(f"❌ ERROR: Baseline directory not found: {self.submitted_dir}")
            return results

        # Process each notebook subdirectory
        for subdir_name, config in self.notebook_outputs.items():
            notebook_name = config['notebook']
            expected_figures = config['figures']

            subdir_path = self.notebooks_dir / subdir_name
            notebook_path = subdir_path / notebook_name

            print(f"\n{'─' * 70}")
            print(f"Testing: {subdir_name}/{notebook_name}")
            print(f"{'─' * 70}")

            if not subdir_path.exists():
                print(f"  ⚠ SKIPPED - Subdirectory not found: {subdir_path}")
                results['details'][subdir_name] = {
                    'status': 'skipped',
                    'reason': 'subdirectory not found'
                }
                continue

            if not notebook_path.exists():
                print(f"  ⚠ SKIPPED - Notebook not found: {notebook_path}")
                results['details'][subdir_name] = {
                    'status': 'skipped',
                    'reason': 'notebook not found'
                }
                continue

            results['notebooks_tested'] += 1
            notebook_result = {
                'execution': {},
                'comparisons': {}
            }

            # Execute notebook
            success, error_msg = self.run_notebook(notebook_path)
            notebook_result['execution'] = {
                'success': success,
                'error': error_msg
            }

            if not success:
                print(f"  ❌ FAILED - {error_msg}")
                results['notebooks_failed'] += 1
                results['details'][subdir_name] = notebook_result
                continue

            print(f"  ✓ Notebook executed successfully")

            # Compare generated figures
            all_figures_good = True
            for figure_name in expected_figures:
                # Figures are generated in the subdirectory
                new_figure = subdir_path / figure_name
                baseline_figure = self.submitted_dir / figure_name

                results['figures_compared'] += 1

                if not new_figure.exists():
                    print(f"\n  {figure_name}:")
                    print(f"    ❌ ERROR: Generated figure not found: {new_figure}")
                    notebook_result['comparisons'][figure_name] = {
                        'error': 'generated figure not found'
                    }
                    all_figures_good = False
                    continue

                if not baseline_figure.exists():
                    print(f"\n  {figure_name}:")
                    print(f"    ❌ ERROR: Baseline figure not found: {baseline_figure}")
                    notebook_result['comparisons'][figure_name] = {
                        'error': 'baseline figure not found'
                    }
                    all_figures_good = False
                    continue

                # Compare images
                metrics = self.compare_images(new_figure, baseline_figure)
                notebook_result['comparisons'][figure_name] = metrics

                if metrics.get('identical', False):
                    results['figures_identical'] += 1

                # Print formatted results
                print(self.format_comparison_result(figure_name, metrics))

                # Check if comparison is acceptable (SSIM > 0.95 or identical)
                if not metrics.get('identical', False):
                    if 'ssim' not in metrics or metrics['ssim'] < 0.95:
                        all_figures_good = False

            # Update overall results
            if all_figures_good:
                results['notebooks_passed'] += 1
            else:
                results['notebooks_failed'] += 1

            results['details'][subdir_name] = notebook_result

        return results

    def print_summary(self, results: Dict):
        """
        Print a summary of test results.

        Args:
            results: Dictionary containing test results
        """
        print("\n" + "=" * 70)
        print("Test Summary")
        print("=" * 70)
        print(f"Notebooks tested: {results['notebooks_tested']}")
        print(f"  Passed: {results['notebooks_passed']}")
        print(f"  Failed: {results['notebooks_failed']}")
        print()
        print(f"Figures compared: {results['figures_compared']}")
        print(f"  Identical: {results['figures_identical']}")
        print(f"  Different: {results['figures_compared'] - results['figures_identical']}")
        print("=" * 70)

        # Exit code
        if results['notebooks_failed'] > 0:
            print("\n❌ TESTS FAILED")
            return 1
        elif results['notebooks_tested'] == 0:
            print("\n⚠ NO TESTS RUN")
            return 1
        else:
            print("\n✓ ALL TESTS PASSED")
            return 0


def main():
    """Main entry point for the test script."""
    # Parse command line arguments
    import argparse

    parser = argparse.ArgumentParser(
        description='Test figure reproducibility by running notebooks and comparing outputs'
    )
    parser.add_argument(
        '--notebooks-dir',
        default='.',
        help='Directory containing notebooks (default: current directory)'
    )
    parser.add_argument(
        '--submitted-dir',
        default='submitted_figures',
        help='Directory containing baseline figures (default: submitted_figures)'
    )
    parser.add_argument(
        '--json-output',
        help='Save detailed results to JSON file'
    )

    args = parser.parse_args()

    # Create tester and run tests
    tester = NotebookFigureTester(
        notebooks_dir=args.notebooks_dir,
        submitted_dir=args.submitted_dir
    )

    results = tester.run_tests()

    # Save JSON output if requested
    if args.json_output:
        with open(args.json_output, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nDetailed results saved to: {args.json_output}")

    # Print summary and exit
    exit_code = tester.print_summary(results)
    sys.exit(exit_code)


if __name__ == '__main__':
    main()
