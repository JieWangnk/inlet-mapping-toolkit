# Inlet Mapping Toolkit

**From 1D Flow Rate to 3D Velocity Profiles for OpenFOAM**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-12-blue.svg)](https://openfoam.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-green.svg)](https://python.org/)

Transform your 1D inlet flow rate data (CSV) into realistic 3D velocity profiles for cardiovascular CFD simulations.

## Quick Start

```bash
# Generate wall-distance profile from flow rate CSV
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --output boundaryData/

# Constant flow with parabolic profile
python map_inlet.py inlet.stl --flowrate 5.0 --profile parabolic --output boundaryData/

# Wall-distance with custom exponent
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --exponent 1.5 --output boundaryData/
```

## Features

- **Multiple Profile Types**:
  - `plug` - Uniform velocity (simplest)
  - `parabolic` - Classic Poiseuille profile
  - `wall_distance` - Distance-to-wall based (irregular geometries)

- **Flexible Input**:
  - Time-varying flow rate from CSV
  - Constant flow rate
  - Velocity data (auto-converted)

- **Auto-Detection**:
  - Inlet boundary extraction
  - Normal direction (inward/outward)
  - Area calculation

- **OpenFOAM Ready**:
  - Generates `boundaryData/<patch>/` structure
  - Compatible with `timeVaryingMappedFixedValue` BC

## Installation

```bash
git clone https://github.com/JieWangnk/inlet-mapping-toolkit.git
cd inlet-mapping-toolkit

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install dependencies
python3 -m pip install -r requirements.txt
```

### Dependencies
- Python 3.8+
- NumPy
- SciPy (for KDTree optimization)
- trimesh (for STL processing)

## Usage

### Basic Usage

```bash
# Wall-distance profile from CSV (recommended for irregular inlets)
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance

# Parabolic profile with constant flow
python map_inlet.py inlet.stl --flowrate 5.0 --profile parabolic

# Custom exponent (n=2 is parabolic-like, n=1/7 is turbulent)
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --exponent 2.0
```

### Advanced Options

```bash
# Specify output directory
python map_inlet.py inlet.stl flowrate.csv --output /path/to/case/constant/boundaryData/inlet

# Generate points file for existing mesh
python map_inlet.py inlet.stl --generate-points-only --output boundaryData/

# Flip normal direction
python map_inlet.py inlet.stl flowrate.csv --flip-normal

# Verbose output
python map_inlet.py inlet.stl flowrate.csv --verbose
```

### CSV Format

```csv
time,flowrate
0.0,0.0
0.1,5.2
0.2,4.8
...
```

Supported columns:
- `time` (required): Time in seconds
- `flowrate`: Flow rate in L/min (auto-detected) or m³/s
- `velocity`: Mean velocity in m/s

## Profile Types

### Plug Profile
Uniform velocity across the inlet. Simple but unrealistic.
```
U(r) = U_mean
```

### Parabolic Profile (Poiseuille)
Classic fully-developed pipe flow profile.
```
U(r) = 2 * U_mean * (1 - (r/R)²)
```

### Wall-Distance Profile
Based on distance from each point to the inlet boundary. Best for irregular geometries.
```
U(d) = U_0 * (d/d_max)^n
```
- `n=1`: Linear
- `n=2`: Quadratic (parabolic-like)
- `n=1/7`: Turbulent power law

## Output Structure

```
boundaryData/
└── inlet/
    ├── points          # Face centre coordinates
    ├── 0.000000/
    │   └── U           # Velocity vectors at t=0
    ├── 0.001000/
    │   └── U           # Velocity vectors at t=0.001
    └── ...
```

## Integration with OpenFOAM

### Boundary Condition (0/U)
```cpp
inlet
{
    type            timeVaryingMappedFixedValue;
    offset          (0 0 0);
    setAverage      false;
}
```

### For Constant Profiles
Even constant profiles use `timeVaryingMappedFixedValue` with a single timestep at t=0.

## Examples

### Example 1: Patient-Specific Irregular Inlet
```bash
python map_inlet.py valve_inlet.stl flowrate.csv \
    --profile wall_distance \
    --exponent 2.0 \
    --output case/constant/boundaryData/inlet
```

### Example 2: Steady-State Validation
```bash
python map_inlet.py inlet.stl \
    --flowrate 6.0 \
    --profile parabolic \
    --output case/constant/boundaryData/inlet
```

### Example 3: Time-Varying Wall-Distance
```bash
python map_inlet.py aorta_inlet.stl cardiac_cycle.csv \
    --profile wall_distance \
    --exponent 2.0 \
    --output case/constant/boundaryData/inlet
```

## API Usage

```python
from inlet_mapper import InletMapper

# Create mapper
mapper = InletMapper(
    stl_path="inlet.stl",
    csv_path="flowrate.csv",
    profile="wall_distance",
    exponent=2.0
)

# Generate boundaryData
mapper.generate(output_dir="boundaryData/inlet")

# Get profile statistics
stats = mapper.get_statistics()
print(f"Max velocity: {stats['max_velocity']:.2f} m/s")
print(f"Mean velocity: {stats['mean_velocity']:.2f} m/s")
```

## Validation

The toolkit has been validated against:
- Analytical parabolic solutions
- Published cardiovascular CFD benchmarks
- OpenFOAM built-in profiles

## Citation

```bibtex
@software{inlet_mapping_toolkit_2025,
  title={Inlet Mapping Toolkit: 1D to 3D Velocity Profiles for OpenFOAM},
  author={Jie Wang},
  year={2025},
  url={https://github.com/JieWangnk/inlet-mapping-toolkit}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

---

**Transform your inlet boundary conditions from simple to realistic!**
