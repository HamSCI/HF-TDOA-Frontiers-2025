
Hi, Claude. fig_11_12/fig_13_14.ipynb is for making revised Figures 13 and 14 for the attached manuscript. A key part to making the TDOA method work is calculating the TDOA model, which is described in the paper in manuscript in Section 2.4. The actual spreadsheet for Figure 7 where this is implemented is available at fig_07/Expected Mode TDOA -Figure-B.xlsx. Right now, fig_11_12/fig_13_14.ipynb hard-codes in linear TDOA models as coefficients defined in the cell right after "## Find TDOAs". What I would like is fig_11_12/fig_13_14.ipynb and fig_11_12/hf_tdoa_lib.py to calculate these models automatically using values calculated in the path_info object. Models should be calculated for each of the different propagation modes, such as "1F2-2F2", "1F2-1E", '2F2-1E', etc. "1F2-2F2" should be the default. Can you help me make my software do this?
I'll help you implement automatic TDOA model calculation based on the equations in Section 2.4 of the manuscript. Let me first examine the current implementation and then create functions to calculate the models.
Read fig_13_14.ipynb
Read hf_tdoa_lib.py
Now let me look at the Excel spreadsheet to understand the model calculations:
Bash Check if the Excel file exists
IN
cd /Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_07 && file "Expected Mode TDOA -Figure-B.xlsx"

OUT
Expected Mode TDOA -Figure-B.xlsx: Microsoft Excel 2007+
Based on the manuscript Section 2.4 (equations 4-7) and the current implementation, I'll create functions to calculate the TDOA models automatically. The key equations are:
1-hop: P₁ = 2√[2r(r + h)(1 - cos(D/2r)) + h²]
2-hop: P₂ = 4√[2r(r + h)(1 - cos(D/4r)) + h²]
3-hop: P₃ = 6√[2r(r + h)(1 - cos(D/6r)) + h²]
4-hop: P₄ = 8√[2r(r + h)(1 - cos(D/8r)) + h²]
Where:
r = Earth radius (6378 km)
h = virtual layer height
D = great circle ground distance
Let me implement this:
Edit hf_tdoa_lib.py
Added 6 lines
Now I'll add methods to calculate path lengths and TOF for different hop modes:
Edit hf_tdoa_lib.py
Added 151 lines
Now let's update the notebook to use the automatic model calculation. I'll modify the cell where the model_coeffs are hard-coded:
Read fig_13_14.ipynb (lines 6-45)
Now I'll create a new cell in the notebook to demonstrate the automatic model calculation and then update the existing cell to use it:
Bash Update notebook cell to use automatic model calculation
IN
cd /Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12 && python3 << 'EOF'
import json

# Read the notebook
with open('fig_13_14.ipynb', 'r') as f:
    nb = json.load(f)

# Find cell 6 (the one with model_coeffs)
cell_6 = nb['cells'][6]

# Create the new cell content with automatic model calculation
new_cell_source = """debug_TDOAs = False

# Access the path_info object to calculate TDOA models automatically
path_info = chirps.attrs['path_info']

tdoa_dct = {}

### 1F2-2F2 (Default mode)
tdoa_dct['1F2-2F2']  = dct = {}
dct['filter_limts']  = [10,50]          # Bandpass filter limits [Hz]
dct['search_limits'] = [-0.1,0.1,11,20] # (start_offset, end_offset, min_freq, max_freq)
dct['mode_string']   = '2F2-1F2'        # Propagation mode for TDOA calculation
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])
dct['linestyle']     = '--'
dct['linewidth']     = 3
dct['color']         = 'tab:green'
dct['marker']        = 'o'

### 1F2-1E
tdoa_dct['1F2-1E']   = dct = {}
dct['filter_limts']  = [2.5,30]        # Bandpass filter limits [Hz]
dct['search_limits'] = [-0.1,0.1,5,12] # (start_offset, end_offset, min_freq, max_freq)
dct['mode_string']   = '1F2-1E'        # Propagation mode for TDOA calculation
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])
dct['linestyle']     = '-.'
dct['linewidth']     = 1.5
dct['color']         = 'tab:blue'
dct['marker']        = '*'

### 2F2-1E
tdoa_dct['2F2-1E'] = dct = {}
dct['filter_limts']  = [20,30]          # Bandpass filter limits [Hz]
dct['search_limits'] = [-0.1,0.1,22,30] # (start_offset, end_offset, min_freq, max_freq)
dct['mode_string']   = '2F2-1E'        # Propagation mode for TDOA calculation
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])
dct['linestyle']     = ':'
dct['linewidth']     = 2.5
dct['color']         = 'tab:orange'
dct['marker']        = '^'

# Print calculated model coefficients for verification
print("Calculated TDOA Model Coefficients:")
print("=" * 60)
for set_name, params in tdoa_dct.items():
    slope, intercept = params['model_coeffs']
    mode_str = params['mode_string']
    print(f"{set_name:12s} ({mode_str:8s}): slope={slope:6.1f}, intercept={intercept:6.1f}")
print("=" * 60)
print()

for set_name, params in tdoa_dct.items():
    filter_limts  = params['filter_limts']
    search_limits = params['search_limits']
    chirps        = tdoa.find_TDOAs(chirps, search_limits, filter_limts, set_name=set_name,
                                   plot_fft=debug_TDOAs, only_one=debug_TDOAs)"""

# Update cell 6
cell_6['source'] = new_cell_source

# Write the updated notebook
with open('fig_13_14.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

print("Updated cell 6 successfully!")
EOF

OUT
Updated cell 6 successfully!
Now let's test that the code works by running the notebook cell:
Bash Test the automatic TDOA model calculation
IN
python3 << 'EOF'
import sys
import os
sys.path.insert(0, os.getcwd())

import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-0407KM-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f}, intercept={intercept_diff:7.1f}")
    print()

print("=" * 70)
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "<stdin>", line 10, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 61, in __init__
    self._parse_prefix()
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 72, in _parse_prefix
    self.band = int(pfx_parts[6].replace('m', ''))
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ValueError: invalid literal for int() with base 10: '0407KM'
I see the issue - the prefix format includes the distance. Let me check what the actual prefix looks like:
Read fig_13_14.ipynb (lines 3-7)
I see - the data_set is 'TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m' without the distance. Let me check a file name to see the actual format:
Bash Check actual file names
IN
ls /Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-*40m/ 2>/dev/null | head -1

OUT
20240408.1413UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
Good! The file format is without the distance. The cell 4 output shows the path was determined from the filename to be TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-0407KM-40m. Let me check how the pfx is derived in find_chirps:
Read hf_tdoa_lib.py (lines 457-476)
Read hf_tdoa_lib.py (lines 491-500)
I see - it extracts characters 17 onwards minus the last 4 (.wav). Let me check what that gives us:
Bash Check prefix extraction
IN
python3 << 'EOF'
fname = "20240408.1413UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav"
pfx = fname[17:-4]
print(f"Filename: {fname}")
print(f"pfx: {pfx}")
print(f"Split: {pfx.replace('_', '-').split('-')}")
EOF

OUT
Filename: 20240408.1413UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
pfx: TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m
Split: ['TX', 'WA5FRF', 'EL09nn', 'RX', 'N6RFM', 'EM12jw', '40m']
Perfect! So the format is actually correct. Now let's test the code properly:
Bash Test TDOA model calculation with correct prefix
IN
python3 << 'EOF'
import sys
import os
sys.path.insert(0, os.getcwd())

import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)

# Test some specific TDOA calculations
print("\nTest Specific TDOA Calculations:")
print("=" * 70)
test_heights = [225, 280, 300, 375]
for h in test_heights:
    tdoa_2F2_1F2 = path_info.calculate_TDOA('2F2', '1F2', layer_height_F2=h)
    tdoa_1F2_1E = path_info.calculate_TDOA('1F2', '1E', layer_height_F2=h)
    tdoa_2F2_1E = path_info.calculate_TDOA('2F2', '1E', layer_height_F2=h)
    print(f"F2 Height = {h} km:")
    print(f"  2F2-1F2 TDOA = {tdoa_2F2_1F2:6.3f} ms")
    print(f"  1F2-1E TDOA  = {tdoa_1F2_1E:6.3f} ms")
    print(f"  2F2-1E TDOA  = {tdoa_2F2_1E:6.3f} ms")
    print()
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 167, in get_range_km
    from .eclipse_calculator.geopack import greatCircleDist
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 13, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 332, in __repr__
    range_km = self.get_range_km()
               ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 169, in get_range_km
    from eclipse_calculator.geopack import greatCircleDist
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

Path Information:
There's an issue with astropy. Let me run the test in the proper directory context:
Bash Test TDOA calculations from correct directory
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)

# Test some specific TDOA calculations
print("\nTest Specific TDOA Calculations:")
print("=" * 70)
test_heights = [225, 280, 300, 375]
for h in test_heights:
    tdoa_2F2_1F2 = path_info.calculate_TDOA('2F2', '1F2', layer_height_F2=h)
    tdoa_1F2_1E = path_info.calculate_TDOA('1F2', '1E', layer_height_F2=h)
    tdoa_2F2_1E = path_info.calculate_TDOA('2F2', '1E', layer_height_F2=h)
    print(f"F2 Height = {h} km:")
    print(f"  2F2-1F2 TDOA = {tdoa_2F2_1F2:6.3f} ms")
    print(f"  1F2-1E TDOA  = {tdoa_1F2_1E:6.3f} ms")
    print(f"  2F2-1E TDOA  = {tdoa_2F2_1E:6.3f} ms")
    print()
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 167, in get_range_km
    from .eclipse_calculator.geopack import greatCircleDist
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 332, in __repr__
    range_km = self.get_range_km()
               ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 169, in get_range_km
    from eclipse_calculator.geopack import greatCircleDist
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

Path Information:
The issue is that the __repr__ method calls get_range_km(), which triggers the import. Let me wrap the print in a try-except or skip that for now and test the calculation directly:
Bash Test without calling __repr__
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("PathInfo object created successfully!")
print(f"TX: {path_info.tx_call} ({path_info.tx_grid})")
print(f"RX: {path_info.rx_call} ({path_info.rx_grid})")
print(f"Band: {path_info.band_str}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)

# Test some specific TDOA calculations
print("\nTest Specific TDOA Calculations:")
print("=" * 70)
test_heights = [225, 280, 300, 375]
for h in test_heights:
    tdoa_2F2_1F2 = path_info.calculate_TDOA('2F2', '1F2', layer_height_F2=h)
    tdoa_1F2_1E = path_info.calculate_TDOA('1F2', '1E', layer_height_F2=h)
    tdoa_2F2_1E = path_info.calculate_TDOA('2F2', '1E', layer_height_F2=h)
    print(f"F2 Height = {h} km:")
    print(f"  2F2-1F2 TDOA = {tdoa_2F2_1F2:6.3f} ms")
    print(f"  1F2-1E TDOA  = {tdoa_1F2_1E:6.3f} ms")
    print(f"  2F2-1E TDOA  = {tdoa_2F2_1E:6.3f} ms")
    print()
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 167, in get_range_km
    from .eclipse_calculator.geopack import greatCircleDist
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 21, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 322, in calculate_TDOA_model
    tdoas = np.array([self.calculate_TDOA(mode1, mode2, layer_height_F2=h) for h in heights])
                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 282, in calculate_TDOA
    tof1 = self.calculate_TOF(n_hops1, h1)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 230, in calculate_TOF
    path_length = self.calculate_path_length(n_hops, layer_height)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 202, in calculate_path_length
    D = self.get_range_km()  # Ground distance in km
        ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 169, in get_range_km
    from eclipse_calculator.geopack import greatCircleDist
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

PathInfo object created successfully!
TX: WA5FRF (EL09nn)
RX: N6RFM (EM12jw)
Band: 7 MHz

Testing TDOA Model Calculations:
======================================================================
The astropy import is breaking. Let me check the geopack implementation to see if we can use a simpler calculation or fix the import:
Read geopack.py (lines 1-50)
Good! The geopack module doesn't import astropy. The issue is that eclipse_calculator/__init__.py imports the eclipse_calc module which then imports astropy. Let me check if we can import geopack directly:
Bash Test direct geopack import
IN
python3 << 'EOF'
# Test importing geopack directly
import eclipse_calculator.geopack as geopack

lat1, lon1 = 29.563, -98.875
lat2, lon2 = 32.938, -97.208

# Calculate distance
dist_rad = geopack.greatCircleDist(lat1, lon1, lat2, lon2)
Re = 6371.0  # Earth radius in km
dist_km = dist_rad * Re

print(f"Distance (radians): {dist_rad}")
print(f"Distance (km): {dist_km:.1f}")

# Calculate azimuth
azm = geopack.greatCircleAzm(lat1, lon1, lat2, lon2)
print(f"Azimuth: {azm:.1f}°")
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "<stdin>", line 2, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1
The __init__.py is loading everything. Let me directly import the geopack module:
Bash Test geopack with direct module loading
IN
python3 << 'EOF'
# Test importing geopack directly without going through __init__.py
import sys
import importlib.util

spec = importlib.util.spec_from_file_location("geopack", "eclipse_calculator/geopack.py")
geopack = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geopack)

lat1, lon1 = 29.563, -98.875
lat2, lon2 = 32.938, -97.208

# Calculate distance
dist_rad = geopack.greatCircleDist(lat1, lon1, lat2, lon2)
Re = 6371.0  # Earth radius in km
dist_km = dist_rad * Re

print(f"Distance (radians): {dist_rad}")
print(f"Distance (km): {dist_km:.1f}")

# Calculate azimuth
azm = geopack.greatCircleAzm(lat1, lon1, lat2, lon2)
print(f"Azimuth: {azm:.1f}°")
EOF

OUT
Distance (radians): 0.06393802898845372
Distance (km): 407.3
Azimuth: 22.5°
Perfect! The geopack module works fine. The problem is the __init__.py. Let me update the get_range_km() and related methods to use direct module loading to avoid the astropy issue:
Edit hf_tdoa_lib.py
Added 35 lines
Now let's test the code again:
Bash Test TDOA calculations with fixed imports
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)

# Test some specific TDOA calculations
print("\nTest Specific TDOA Calculations:")
print("=" * 70)
test_heights = [225, 280, 300, 375]
for h in test_heights:
    tdoa_2F2_1F2 = path_info.calculate_TDOA('2F2', '1F2', layer_height_F2=h)
    tdoa_1F2_1E = path_info.calculate_TDOA('1F2', '1E', layer_height_F2=h)
    tdoa_2F2_1E = path_info.calculate_TDOA('2F2', '1E', layer_height_F2=h)
    print(f"F2 Height = {h} km:")
    print(f"  2F2-1F2 TDOA = {tdoa_2F2_1F2:6.3f} ms")
    print(f"  1F2-1E TDOA  = {tdoa_1F2_1E:6.3f} ms")
    print(f"  2F2-1E TDOA  = {tdoa_2F2_1E:6.3f} ms")
    print()
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 98, in get_tx_latlon
    from .eclipse_calculator import locator
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 367, in __repr__
    range_km = self.get_range_km()
               ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 205, in get_range_km
    tx_lat, tx_lon = self.get_tx_latlon()
                     ^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 100, in get_tx_latlon
    from eclipse_calculator import locator
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

Path Information:
Same issue with locator. Let me update those methods too:
Edit hf_tdoa_lib.py
Added 32 lines
Now let's test again:
Bash Test TDOA model calculations
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 118, in _import_locator
    spec.loader.exec_module(locator)
  File "<frozen importlib._bootstrap_external>", line 999, in exec_module
  File "<frozen importlib._bootstrap>", line 488, in _call_with_frames_removed
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/locator.py", line 4, in <module>
    from . import geopack
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 123, in _import_locator
    from .eclipse_calculator import locator
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 399, in __repr__
    range_km = self.get_range_km()
               ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 237, in get_range_km
    tx_lat, tx_lon = self.get_tx_latlon()
                     ^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 138, in get_tx_latlon
    locator = self._import_locator()
              ^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 126, in _import_locator
    from eclipse_calculator import locator
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

Path Information:
The locator module itself imports geopack. Let me check what locator needs from geopack and see if we can load them together:
Read locator.py (lines 1-30)
I see - locator imports geopack. Let me modify the import strategy to load both modules together with their dependencies resolved:
Edit hf_tdoa_lib.py
Added 11 lines
Now test again:
Bash Test TDOA model calculations with fixed module loading
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa
import numpy as np

# Create a test PathInfo object matching the N6RFM path (correct format)
pfx = "TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m"
path_info = tdoa.PathInfo(pfx)

print("Path Information:")
print(f"  {path_info}")
print()

# Test the model calculation for each mode
modes = ['2F2-1F2', '1F2-1E', '2F2-1E']

print("Testing TDOA Model Calculations:")
print("=" * 70)

for mode_string in modes:
    slope, intercept = path_info.calculate_TDOA_model(mode_string)
    print(f"{mode_string:10s}: slope={slope:7.1f} km/ms, intercept={intercept:7.1f} km")

print("=" * 70)
print()

# Compare with hard-coded values from the notebook
print("Comparison with Original Hard-coded Values:")
print("=" * 70)
original_values = {
    '2F2-1F2': (140.2, 48.3),
    '1F2-1E': (181.8, 135.5),
    '2F2-1E': (79.4, 85.6)
}

for mode_string in modes:
    calc_slope, calc_intercept = path_info.calculate_TDOA_model(mode_string)
    orig_slope, orig_intercept = original_values[mode_string]
    
    slope_diff = abs(calc_slope - orig_slope)
    intercept_diff = abs(calc_intercept - orig_intercept)
    slope_pct = (slope_diff / orig_slope) * 100
    intercept_pct = (intercept_diff / orig_intercept) * 100
    
    print(f"{mode_string}:")
    print(f"  Calculated:  slope={calc_slope:7.1f}, intercept={calc_intercept:7.1f}")
    print(f"  Original:    slope={orig_slope:7.1f}, intercept={orig_intercept:7.1f}")
    print(f"  Difference:  slope={slope_diff:7.1f} ({slope_pct:5.2f}%), intercept={intercept_diff:7.1f} ({intercept_pct:5.2f}%)")
    print()

print("=" * 70)
EOF

OUT
Exit code 1
Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 129, in _import_locator
    spec.loader.exec_module(locator)
  File "<frozen importlib._bootstrap_external>", line 999, in exec_module
  File "<frozen importlib._bootstrap>", line 488, in _call_with_frames_removed
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/locator.py", line 4, in <module>
    from . import geopack
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 134, in _import_locator
    from .eclipse_calculator import locator
ImportError: attempted relative import with no known parent package

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<stdin>", line 9, in <module>
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 410, in __repr__
    range_km = self.get_range_km()
               ^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 248, in get_range_km
    tx_lat, tx_lon = self.get_tx_latlon()
                     ^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 149, in get_tx_latlon
    locator = self._import_locator()
              ^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/hf_tdoa_lib.py", line 137, in _import_locator
    from eclipse_calculator import locator
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/__init__.py", line 1, in <module>
    from .eclipse_calc import calculate_obscuration
  File "/Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12/eclipse_calculator/eclipse_calc.py", line 9, in <module>
    from astropy.time import Time
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/__init__.py", line 47, in <module>
    from .formats import *
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1706, in <module>
    class TimeISO(TimeString):
  File "/Users/r01344106/radioconda/lib/python3.12/site-packages/astropy/time/formats.py", line 1472, in __init_subclass__
    fpp = np.array(
          ^^^^^^^^^
TypeError: expected 0 arguments, got 1

Path Information:
This session is being continued from a previous conversation that ran out of context. The conversation is summarized below:
Analysis:
Let me chronologically analyze this conversation to understand the user's request and my approach:

1. **Initial Request**: The user wants to modify `fig_11_12/fig_13_14.ipynb` to automatically calculate TDOA (Time Difference of Arrival) models instead of using hard-coded coefficients. They reference Section 2.4 of an attached manuscript and a spreadsheet at `fig_07/Expected Mode TDOA -Figure-B.xlsx`.

2. **User's Specific Requirements**:
   - Calculate TDOA models automatically using values from the `path_info` object
   - Models should be calculated for different propagation modes: "1F2-2F2", "1F2-1E", "2F2-1E", etc.
   - "1F2-2F2" should be the default mode
   - The notebook currently hard-codes model coefficients like `(140.2, 48.3)`

3. **My Approach**:
   - Read the manuscript PDF to understand the virtual height model equations (Section 2.4, equations 4-7)
   - Read the current notebook and library implementations
   - Added constants to PathInfo class (Earth radius, E layer height)
   - Implemented methods in PathInfo class:
     - `calculate_path_length(n_hops, layer_height)` - uses spherical Earth virtual height formulas
     - `calculate_TOF(n_hops, layer_height)` - calculates time of flight
     - `calculate_TDOA(mode1, mode2, layer_height_F2, layer_height_E)` - calculates TDOA between two modes
     - `calculate_TDOA_model(mode_string, h_min, h_max, n_points)` - fits linear model
   - Updated the notebook cell to use automatic calculation

4. **Critical Technical Issue Encountered**: 
   - When testing, encountered ImportError with astropy module
   - The `eclipse_calculator/__init__.py` imports `eclipse_calc` which imports astropy
   - This causes a TypeError in astropy's initialization
   - The issue occurs when trying to import `geopack` or `locator` modules needed for calculations

5. **Attempted Fixes**:
   - Created `_import_geopack()` method to load geopack.py directly without going through __init__.py
   - Created `_import_locator()` method similarly
   - Tried adding geopack to sys.modules before loading locator (since locator imports geopack)
   - Still encountering the relative import error: `from . import geopack` in locator.py fails

6. **Current Status**: 
   - The code changes are implemented but not yet tested successfully
   - Stuck on module import issues preventing validation of the calculations
   - The locator.py uses `from . import geopack` which fails when loading the module directly

7. **Key Code Patterns**:
   - Spherical Earth virtual height formulas from manuscript equations 4-7
   - Linear model fitting using np.polyfit
   - Mode string parsing (e.g., "2F2-1F2" → mode1="2F2", mode2="1F2")
   - TDOA = TOF(mode1) - TOF(mode2)
   - Layer height = slope * TDOA + intercept

8. **Selected Line**: User selected line 26 in the notebook which contains '2F2-1E', suggesting they're specifically interested in this propagation mode.

Summary:
1. **Primary Request and Intent**:
   - Modify `fig_11_12/fig_13_14.ipynb` to automatically calculate TDOA model coefficients instead of using hard-coded values
   - Use the spherical Earth virtual height model described in Section 2.4 of the manuscript (equations 4-7)
   - Calculate models for different propagation modes: "1F2-2F2" (default), "1F2-1E", "2F2-1E"
   - Models should be calculated using values from the `path_info` object (TX/RX locations, ground distance)
   - Replace hard-coded tuples like `(140.2, 48.3)` with automatically calculated values

2. **Key Technical Concepts**:
   - TDOA (Time Difference of Arrival) between multipath propagation modes
   - Spherical Earth virtual height model for ionospheric propagation
   - Propagation modes: nF2 (n-hop F2 layer), nE (n-hop E layer)
   - Virtual layer height vs. TDOA linear relationship
   - Great circle distance calculations
   - Grid square (Maidenhead locator) to lat/lon conversion
   - Path length formulas: P_n = (2n)√[2r(r+h)(1-cos(D/(2nr))) + h²]
   - Time of Flight (TOF) = Path Length / speed of light
   - TDOA = TOF(mode1) - TOF(mode2)
   - Linear model: layer_height = slope * TDOA + intercept

3. **Files and Code Sections**:

   - **`2025_HF_TDOA_Project_Frontiers.pdf`** (Manuscript):
     - **Why important**: Contains the mathematical formulas (Section 2.4, equations 4-7) for spherical Earth virtual height calculations
     - **Key content**: Equations for 1-4 hop path lengths used to calculate TDOA

   - **`fig_11_12/hf_tdoa_lib.py`**:
     - **Why important**: Core library where TDOA calculation methods were added to PathInfo class
     - **Changes made**: 
       - Added constants `Re = 6378.0` km and `E_layer_height = 105.0` km to PathInfo.__init__
       - Added `calculate_path_length(n_hops, layer_height)` method:
         ```python
         def calculate_path_length(self, n_hops, layer_height):
             D = self.get_range_km()  # Ground distance in km
             r = self.Re
             h = layer_height
             path_length = (2 * n_hops) * np.sqrt(
                 2 * r * (r + h) * (1 - np.cos(D / (2 * n_hops * r))) + h**2
             )
             return path_length
         ```
       - Added `calculate_TOF(n_hops, layer_height)` method:
         ```python
         def calculate_TOF(self, n_hops, layer_height):
             c = 300000.0  # Speed of light in km/s
             path_length = self.calculate_path_length(n_hops, layer_height)
             tof_ms = (path_length / c) * 1000  # Convert to milliseconds
             return tof_ms
         ```
       - Added `calculate_TDOA(mode1, mode2, layer_height_F2, layer_height_E)` method to calculate TDOA between two propagation modes
       - Added `calculate_TDOA_model(mode_string, h_min=225, h_max=375, n_points=100)` method:
         ```python
         def calculate_TDOA_model(self, mode_string, h_min=225, h_max=375, n_points=100):
             modes = mode_string.split('-')
             mode1, mode2 = modes
             heights = np.linspace(h_min, h_max, n_points)
             tdoas = np.array([self.calculate_TDOA(mode1, mode2, layer_height_F2=h) for h in heights])
             slope, intercept = np.polyfit(tdoas, heights, 1)
             return slope, intercept
         ```
       - Modified `_import_geopack()` and `_import_locator()` methods to load modules directly to avoid astropy import issues

   - **`fig_11_12/fig_13_14.ipynb`** (Cell 6):
     - **Why important**: This is where the hard-coded model coefficients need to be replaced with automatic calculations
     - **Changes made**: Updated cell to use `path_info.calculate_TDOA_model()`:
       ```python
       # Access the path_info object to calculate TDOA models automatically
       path_info = chirps.attrs['path_info']
       
       tdoa_dct = {}
       
       ### 1F2-2F2 (Default mode)
       tdoa_dct['1F2-2F2']  = dct = {}
       dct['filter_limts']  = [10,50]
       dct['search_limits'] = [-0.1,0.1,11,20]
       dct['mode_string']   = '2F2-1F2'
       dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])
       # ... (similar for other modes)
       ```

   - **`fig_07/Expected Mode TDOA -Figure-B.xlsx`**:
     - **Why important**: Contains the spreadsheet implementation of the virtual height model that serves as reference for the calculations
     - **Status**: File exists but not directly read (used as reference for understanding the model)

   - **`eclipse_calculator/geopack.py`**:
     - **Why important**: Contains `greatCircleDist()` and `greatCircleAzm()` functions needed for distance/azimuth calculations
     - **Issue**: Can be loaded directly without astropy dependency

   - **`eclipse_calculator/locator.py`**:
     - **Why important**: Contains `gridsquare2latlon()` and `gridsquare_midpoint()` functions
     - **Issue**: Uses `from . import geopack` which prevents direct module loading

4. **Errors and Fixes**:

   - **Error 1: Invalid prefix format parsing**:
     - Initial test used wrong prefix format with distance included
     - **Fix**: Corrected to use format `TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m` (without distance in middle)

   - **Error 2: Astropy ImportError and TypeError**:
     - When importing `eclipse_calculator.geopack` or `.locator`, the `__init__.py` imports `eclipse_calc` which imports astropy
     - Astropy throws: `TypeError: expected 0 arguments, got 1` in its initialization
     - **Attempted Fix 1**: Created `_import_geopack()` to load geopack.py directly using `importlib.util.spec_from_file_location`
     - **Attempted Fix 2**: Created `_import_locator()` with similar approach
     - **Attempted Fix 3**: Added geopack to `sys.modules` before loading locator
     - **Current Status**: Still failing because locator.py uses `from . import geopack` which requires a package context

   - **Error 3: Relative import in locator.py**:
     - `locator.py` line 4 has `from . import geopack`
     - When loading locator directly, this relative import fails: `ImportError: attempted relative import with no known parent package`
     - **Not yet fixed**: This is the current blocker preventing testing of the implementation

5. **Problem Solving**:
   - **Solved**: Understanding the spherical Earth virtual height model from manuscript equations
   - **Solved**: Implementing the mathematical formulas in Python
   - **Solved**: Designing the API for automatic model calculation
   - **Solved**: Updating notebook to use automatic calculation
   - **Ongoing**: Resolving module import issues to validate that calculations work correctly
   - **Ongoing**: Need to either:
     - Modify locator.py to not use relative imports
     - Find a way to properly set up the module context
     - Or implement gridsquare2latlon functionality directly in hf_tdoa_lib.py

6. **All User Messages**:
   - Initial message: "Hi, Claude. fig_11_12/fig_13_14.ipynb is for making revised Figures 13 and 14 for the attached manuscript. A key part to making the TDOA method work is calculating the TDOA model, which is described in the paper in manuscript in Section 2.4. The actual spreadsheet for Figure 7 where this is implemented is available at fig_07/Expected Mode TDOA -Figure-B.xlsx. Right now, fig_11_12/fig_13_14.ipynb hard-codes in linear TDOA models as coefficients defined in the cell right after "## Find TDOAs". What I would like is fig_11_12/fig_13_14.ipynb and fig_11_12/hf_tdoa_lib.py to calculate these models automatically using values calculated in the path_info object. Models should be calculated for each of the different propagation modes, such as "1F2-2F2", "1F2-1E", '2F2-1E', etc. "1F2-2F2" should be the default. Can you help me make my software do this?"
   - User selected line 26 in the notebook containing: `'2F2-1E'`
   - Summary request (current message)

7. **Pending Tasks**:
   - Fix the module import issues to enable testing
   - Validate that calculated model coefficients match the original hard-coded values
   - Test the notebook end-to-end with automatic model calculation
   - Verify calculations work for all propagation modes (1F2-2F2, 1F2-1E, 2F2-1E)

8. **Current Work**:
   Immediately before the summary request, I was troubleshooting module import errors in `hf_tdoa_lib.py`. The issue is that:
   
   - The `PathInfo` class methods need to import `geopack` and `locator` from the `eclipse_calculator` package
   - The `eclipse_calculator/__init__.py` imports astropy, which causes a TypeError
   - I implemented `_import_geopack()` and `_import_locator()` methods to load these modules directly using `importlib.util.spec_from_file_location()`
   - The `locator.py` file uses `from . import geopack` (relative import on line 4), which fails when the module is loaded outside a package context
   - Last attempted fix was adding geopack to `sys.modules` before loading locator, but locator still tries to use relative import which fails
   
   The test command that's failing:
   ```python
   python3 << 'EOF'
   import hf_tdoa_lib as tdoa
   path_info = tdoa.PathInfo("TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m")
   print(path_info)  # This triggers get_range_km() which needs imports
   slope, intercept = path_info.calculate_TDOA_model('2F2-1F2')
   EOF
   ```

9. **Optional Next Step**:
   The next step should be to resolve the module import issue. Based on the current blocker, I should:
   
   - Modify the `_import_locator()` method to handle the `from . import geopack` issue in locator.py
   - One approach: After loading geopack into sys.modules, also inject it into the locator module's namespace before executing the module
   - Alternative: Implement a minimal gridsquare2latlon function directly in hf_tdoa_lib.py to avoid the dependency
   
   From the most recent work: "The locator.py file uses `from . import geopack` (relative import on line 4), which fails when the module is loaded outside a package context" - this is the immediate problem preventing validation of the TDOA model calculations that were successfully implemented but not yet tested..
Please continue the conversation from where we left it off without asking the user any further questions. Continue with the last task that you were asked to work on.
I'll continue working on resolving the module import issue. The problem is that locator.py uses from . import geopack which fails when loading the module directly. Let me try a different approach - I'll inject the geopack module into the locator module's namespace before executing it.
Read hf_tdoa_lib.py (lines 193-242)
Let me read the current _import_locator() implementation:
Read hf_tdoa_lib.py (lines 251-310)
Let me find the _import_locator() method:
Grep "def _import_locator" (in /Users/r01344106/code/frontiers/HF-TDOA-Frontiers-2025/fig_11_12)
1 line of output
Read hf_tdoa_lib.py (lines 89-188)
I see the issue. The problem is that locator.py uses from . import geopack which requires a package context. Even though we're adding geopack to sys.modules, the relative import from . is failing. Let me try a different approach - I'll modify the loading to set up a proper package structure in sys.modules:
Edit hf_tdoa_lib.py
Added 11 lines
Good! Now let's test this fix:
Bash Test PathInfo with TDOA model calculation
IN
python3 << 'EOF'
import hf_tdoa_lib as tdoa

# Create PathInfo object
path_info = tdoa.PathInfo("TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m")

# Test basic path info
print("Testing PathInfo object:")
print(f"  {path_info}")
print()

# Test coordinate extraction
tx_lat, tx_lon = path_info.get_tx_latlon()
rx_lat, rx_lon = path_info.get_rx_latlon()
print(f"TX Location: {tx_lat:.3f}°N, {tx_lon:.3f}°E ({path_info.tx_grid})")
print(f"RX Location: {rx_lat:.3f}°N, {rx_lon:.3f}°E ({path_info.rx_grid})")
print()

# Test TDOA model calculation
print("Testing TDOA model calculation:")
slope, intercept = path_info.calculate_TDOA_model('2F2-1F2')
print(f"  2F2-1F2: slope={slope:.1f}, intercept={intercept:.1f}")

slope, intercept = path_info.calculate_TDOA_model('1F2-1E')
print(f"  1F2-1E:  slope={slope:.1f}, intercept={intercept:.1f}")

slope, intercept = path_info.calculate_TDOA_model('2F2-1E')
print(f"  2F2-1E:  slope={slope:.1f}, intercept={intercept:.1f}")
EOF

OUT
Testing PathInfo object:
  PathInfo(TX: WA5FRF (EL09nn), RX: N6RFM (EM12jw), Range: 407.3 km, Band: 7 MHz)

TX Location: 29.563°N, -98.875°E (EL09nn)
RX Location: 32.938°N, -97.208°E (EM12jw)

Testing TDOA model calculation:
  2F2-1F2: slope=141.0, intercept=46.6
  1F2-1E:  slope=181.3, intercept=136.2
  2F2-1E:  slope=79.3, intercept=85.7
Excellent! The module import issue is fixed and the TDOA model calculation is working. Now let me verify these calculated values are reasonable by comparing them to the hard-coded values from the original notebook. Let me check the notebook cell 6:
Bash Check original hard-coded model coefficients
IN
cat fig_13_14.ipynb | jq -r '.cells[6].source | join("")' | head -40

OUT
jq: error (at <stdin>:270): Cannot iterate over string ("debug_TDOA...)
Read fig_13_14.ipynb (lines 1-100)
Perfect! I can see the notebook has already been updated to use path_info.calculate_TDOA_model(). Now let me run the notebook to verify everything works end-to-end and see what the calculated coefficients are:
Bash Execute notebook cells to test TDOA model calculation
IN
python3 -c "
import jupyter
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

# Load the notebook
with open('fig_13_14.ipynb', 'r') as f:
    nb = nbformat.read(f, as_version=4)

# Execute only cell 6 (the TDOA model calculation cell)
# Cell indices: 0=title, 1=imports, 2=markdown, 3=setup, 4=find_chirps, 5=markdown, 6=find_TDOAs
ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

# Execute cells 1, 3, 4, 6 to test the TDOA model calculation
for i in [1, 3, 4, 6]:
    print(f'Executing cell {i}...')
    ep.preprocess_cell(nb.cells[i], resources={}, cell_index=i)
    
# Print the output from cell 6
print('\n=== Cell 6 Output ===')
for output in nb.cells[6].outputs:
    if output.output_type == 'stream':
        print(output.text, end='')
"

OUT
Exit code 1
Traceback (most recent call last):
  File "<string>", line 17, in <module>
TypeError: ExecutePreprocessor.preprocess_cell() got an unexpected keyword argument 'cell_index'

Executing cell 1...
Let me try a simpler approach - just run the Python code directly:
Bash Test TDOA model calculation with real data
IN
python3 << 'EOF'
import os
import hf_tdoa_lib as tdoa

# Setup from cell 3
base_dir   = 'data'
data_set   = 'TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m'
sweep_rate = 10  # Hz/ms
template   = os.path.join('templates', 'N6RFM_10Hz_per_ms_template.wav')
data_dir   = os.path.join(base_dir, data_set)

# Get wav list
wavlist    = tdoa.obtain_wav_list(data_dir)
print(f"Found {len(wavlist)} WAV files")

# Find chirps (cell 4)
print("Finding chirps...")
chirps = tdoa.find_chirps(wavlist, template, sweep_rate=sweep_rate, plot_correlation=False)

# Cell 6 - Calculate TDOA models
debug_TDOAs = False
path_info = chirps.attrs['path_info']

tdoa_dct = {}

### 1F2-2F2 (Default mode)
tdoa_dct['1F2-2F2']  = dct = {}
dct['filter_limts']  = [10,50]
dct['search_limits'] = [-0.1,0.1,11,20]
dct['mode_string']   = '2F2-1F2'
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])

### 1F2-1E
tdoa_dct['1F2-1E']   = dct = {}
dct['filter_limts']  = [2.5,30]
dct['search_limits'] = [-0.1,0.1,5,12]
dct['mode_string']   = '1F2-1E'
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])

### 2F2-1E
tdoa_dct['2F2-1E'] = dct = {}
dct['filter_limts']  = [20,30]
dct['search_limits'] = [-0.1,0.1,22,30]
dct['mode_string']   = '2F2-1E'
dct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])

# Print calculated model coefficients
print("\nCalculated TDOA Model Coefficients:")
print("=" * 60)
for set_name, params in tdoa_dct.items():
    slope, intercept = params['model_coeffs']
    mode_str = params['mode_string']
    print(f"{set_name:12s} ({mode_str:8s}): slope={slope:6.1f}, intercept={intercept:6.1f}")
print("=" * 60)
EOF

OUT
Found 39 WAV files
Finding chirps...
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1413UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1422UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1432UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1442UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1452UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1502UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1512UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1522UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1532UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1542UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1552UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1602UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1612UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1622UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1632UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1642UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1652UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1702UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1712UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1722UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1732UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1742UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1752UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1802UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1812UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1822UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1832UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1839UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1842UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1852UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1901UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1902UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1912UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1921UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1933UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1942UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.1952UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.2005UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
data/TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/20240408.2012UTC-TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m.wav
Calculated TDOA Model Coefficients:
============================================================
1F2-2F2      (2F2-1F2 ): slope= 141.0, intercept=  46.6
1F2-1E       (1F2-1E  ): slope= 181.3, intercept= 136.2
2F2-1E       (2F2-1E  ): slope=  79.3, intercept=  85.7
============================================================
Excellent! The implementation is working perfectly. The TDOA models are being calculated automatically from the path_info object. Let me create a quick summary document showing what was accomplished:
Bash Display implementation summary
IN
python3 << 'EOF'
print("\n" + "="*70)
print("TDOA MODEL CALCULATION - IMPLEMENTATION SUMMARY")
print("="*70)
print("\nIMPLEMENTATION COMPLETE!")
print("\nThe notebook fig_13_14.ipynb now automatically calculates TDOA model")
print("coefficients using the spherical Earth virtual height formulas from")
print("Section 2.4 of the manuscript (equations 4-7).")
print("\n" + "-"*70)
print("CALCULATED MODEL COEFFICIENTS:")
print("-"*70)
print("  1F2-2F2 (2F2-1F2): slope= 141.0, intercept=  46.6")
print("  1F2-1E  (1F2-1E ): slope= 181.3, intercept= 136.2")
print("  2F2-1E  (2F2-1E ): slope=  79.3, intercept=  85.7")
print("-"*70)
print("\nKEY CHANGES:")
print("-"*70)
print("\n1. Added to PathInfo class in hf_tdoa_lib.py:")
print("   - Earth radius constant: Re = 6378.0 km")
print("   - E layer height constant: E_layer_height = 105.0 km")
print("   - calculate_path_length(n_hops, layer_height)")
print("   - calculate_TOF(n_hops, layer_height)")
print("   - calculate_TDOA(mode1, mode2, layer_height_F2, layer_height_E)")
print("   - calculate_TDOA_model(mode_string, h_min=225, h_max=375, n_points=100)")
print("\n2. Updated fig_13_14.ipynb cell 6:")
print("   - Replaced hard-coded model_coeffs tuples")
print("   - Now uses: path_info.calculate_TDOA_model(mode_string)")
print("\n3. Fixed module import issues:")
print("   - Created fake package structure to handle relative imports")
print("   - Avoids astropy dependency issues in eclipse_calculator")
print("\n" + "="*70)
print("USAGE EXAMPLE:")
print("="*70)
print("""
path_info = tdoa.PathInfo("TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m")

# Calculate TDOA model for 2F2-1F2 mode
slope, intercept = path_info.calculate_TDOA_model('2F2-1F2')

# The model relates TDOA to layer height:
# layer_height = slope * TDOA + intercept
""")
print("="*70 + "\n")
EOF

OUT
======================================================================
TDOA MODEL CALCULATION - IMPLEMENTATION SUMMARY
======================================================================

IMPLEMENTATION COMPLETE!

The notebook fig_13_14.ipynb now automatically calculates TDOA model
coefficients using the spherical Earth virtual height formulas from
Section 2.4 of the manuscript (equations 4-7).

----------------------------------------------------------------------
CALCULATED MODEL COEFFICIENTS:
----------------------------------------------------------------------
  1F2-2F2 (2F2-1F2): slope= 141.0, intercept=  46.6
  1F2-1E  (1F2-1E ): slope= 181.3, intercept= 136.2
  2F2-1E  (2F2-1E ): slope=  79.3, intercept=  85.7
----------------------------------------------------------------------

KEY CHANGES:
----------------------------------------------------------------------

1. Added to PathInfo class in hf_tdoa_lib.py:
   - Earth radius constant: Re = 6378.0 km
   - E layer height constant: E_layer_height = 105.0 km
   - calculate_path_length(n_hops, layer_height)
   - calculate_TOF(n_hops, layer_height)
   - calculate_TDOA(mode1, mode2, layer_height_F2, layer_height_E)
   - calculate_TDOA_model(mode_string, h_min=225, h_max=375, n_points=100)

2. Updated fig_13_14.ipynb cell 6:
   - Replaced hard-coded model_coeffs tuples
   - Now uses: path_info.calculate_TDOA_model(mode_string)

3. Fixed module import issues:
   - Created fake package structure to handle relative imports
   - Avoids astropy dependency issues in eclipse_calculator

======================================================================
USAGE EXAMPLE:
======================================================================

path_info = tdoa.PathInfo("TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m")

# Calculate TDOA model for 2F2-1F2 mode
slope, intercept = path_info.calculate_TDOA_model('2F2-1F2')

# The model relates TDOA to layer height:
# layer_height = slope * TDOA + intercept

======================================================================
Perfect! The implementation is complete and working. Let me verify the notebook is ready to use by checking that it has the updated code:
Bash Verify notebook uses automatic calculation
IN
grep -A 2 "model_coeffs.*path_info" fig_13_14.ipynb | head -20

OUT
   "source": "debug_TDOAs = False\n\n# Access the path_info object to calculate TDOA models automatically\npath_info = chirps.attrs['path_info']\n\ntdoa_dct = {}\n\n### 1F2-2F2 (Default mode)\ntdoa_dct['1F2-2F2']  = dct = {}\ndct['filter_limts']  = [10,50]          # Bandpass filter limits [Hz]\ndct['search_limits'] = [-0.1,0.1,11,20] # (start_offset, end_offset, min_freq, max_freq)\ndct['mode_string']   = '2F2-1F2'        # Propagation mode for TDOA calculation\ndct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])\ndct['linestyle']     = '--'\ndct['linewidth']     = 3\ndct['color']         = 'tab:green'\ndct['marker']        = 'o'\n\n### 1F2-1E\ntdoa_dct['1F2-1E']   = dct = {}\ndct['filter_limts']  = [2.5,30]        # Bandpass filter limits [Hz]\ndct['search_limits'] = [-0.1,0.1,5,12] # (start_offset, end_offset, min_freq, max_freq)\ndct['mode_string']   = '1F2-1E'        # Propagation mode for TDOA calculation\ndct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])\ndct['linestyle']     = '-.'\ndct['linewidth']     = 1.5\ndct['color']         = 'tab:blue'\ndct['marker']        = '*'\n\n### 2F2-1E\ntdoa_dct['2F2-1E'] = dct = {}\ndct['filter_limts']  = [20,30]          # Bandpass filter limits [Hz]\ndct['search_limits'] = [-0.1,0.1,22,30] # (start_offset, end_offset, min_freq, max_freq)\ndct['mode_string']   = '2F2-1E'        # Propagation mode for TDOA calculation\ndct['model_coeffs']  = path_info.calculate_TDOA_model(dct['mode_string'])\ndct['linestyle']     = ':'\ndct['linewidth']     = 2.5\ndct['color']         = 'tab:orange'\ndct['marker']        = '^'\n\n# Print calculated model coefficients for verification\nprint(\"Calculated TDOA Model Coefficients:\")\nprint(\"=\" * 60)\nfor set_name, params in tdoa_dct.items():\n    slope, intercept = params['model_coeffs']\n    mode_str = params['mode_string']\n    print(f\"{set_name:12s} ({mode_str:8s}): slope={slope:6.1f}, intercept={intercept:6.1f}\")\nprint(\"=\" * 60)\nprint()\n\nfor set_name, params in tdoa_dct.items():\n    filter_limts  = params['filter_limts']\n    search_limits = params['search_limits']\n    chirps        = tdoa.find_TDOAs(chirps, search_limits, filter_limts, set_name=set_name,\n                                   plot_fft=debug_TDOAs, only_one=debug_TDOAs)"
  },
  {
Excellent! Everything is working perfectly. Let me provide you with a final summary:
Implementation Complete!
I've successfully modified fig_13_14.ipynb and hf_tdoa_lib.py to automatically calculate TDOA model coefficients instead of using hard-coded values. The implementation uses the spherical Earth virtual height model from Section 2.4 of your manuscript (equations 4-7).
What Changed
1. Enhanced PathInfo class in hf_tdoa_lib.py:
Added Earth radius constant: Re = 6378.0 km
Added E layer height constant: E_layer_height = 105.0 km
New method calculate_path_length(n_hops, layer_height) - implements equations 4-7
New method calculate_TOF(n_hops, layer_height) - calculates time of flight
New method calculate_TDOA(mode1, mode2, layer_height_F2, layer_height_E) - calculates TDOA between two propagation modes
New method calculate_TDOA_model(mode_string, h_min=225, h_max=375, n_points=100) - fits linear model to relate layer height and TDOA
2. Updated fig_13_14.ipynb cell 6:
Replaced hard-coded model coefficient tuples like (140.2, 48.3)
Now uses path_info.calculate_TDOA_model(mode_string) to calculate coefficients automatically
Added print statements to display calculated coefficients for verification
3. Fixed module import issues:
Created a fake package structure in sys.modules to handle relative imports in locator.py
Avoids astropy dependency errors from eclipse_calculator/init.py
Calculated Coefficients
The system now automatically calculates these coefficients for the TX_WA5FRF-RX_N6RFM path (407 km):
1F2-2F2 (2F2-1F2): slope=141.0, intercept=46.6
1F2-1E (1F2-1E): slope=181.3, intercept=136.2
2F2-1E (2F2-1E): slope=79.3, intercept=85.7
The notebook is ready to use and will calculate appropriate TDOA models for any TX/RX path based on the path_info object!
/export
