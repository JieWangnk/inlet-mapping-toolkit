# Inlet Mapping Toolkit

**From 1D Flow Rate to 3D Velocity Profiles for OpenFOAM**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-12-blue.svg)](https://openfoam.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-green.svg)](https://python.org/)

Transform 1D flow rate data into realistic 3D velocity profiles for cardiovascular CFD simulations.

## Installation

```bash
git clone https://github.com/JieWangnk/inlet-mapping-toolkit.git
cd inlet-mapping-toolkit
pip install -r requirements.txt
```

**Dependencies:** Python 3.8+, NumPy, SciPy, trimesh

## Quick Start

```bash
# Time-varying flow from CSV
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance

# Constant flow rate (5 L/min)
python map_inlet.py inlet.stl --flowrate 5.0 --profile parabolic

# Using existing mesh points
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance \
    --points-file constant/boundaryData/inlet/points \
    --output constant/boundaryData/inlet
```

## Profile Types

| Profile | Description | Best For |
|---------|-------------|----------|
| `plug` | Uniform velocity | Testing |
| `parabolic` | Classic Poiseuille (1 - r²/R²) | Circular inlets |
| `wall_distance` | Based on distance to wall | Irregular geometries |

## Examples

The `examples/` folder contains sample files. Example outputs are included:

```bash
# Constant flow - using STL face centres (for testing/visualization)
python map_inlet.py ./examples/inlet.stl --flowrate 4.7 --profile wall_distance --output boundaryData_stl/

# Constant flow - using mesh face centres (recommended for production)
python map_inlet.py ./examples/inlet.stl --flowrate 4.7 --profile wall_distance --points-file ./examples/points --output boundaryData_mesh/

# Time-varying flow from CSV
python map_inlet.py ./examples/inlet.stl ./examples/flowrate.csv --profile wall_distance --points-file ./examples/points --output boundaryData_timevarying_mesh/
```

- `boundaryData_stl/` - constant flow, mapped to STL triangle centres
- `boundaryData_mesh/` - constant flow, mapped to mesh face centres
- `boundaryData_timevarying_mesh/` - time-varying flow from CSV, mapped to mesh face centres

```bash
# Wall-distance with custom exponent (1/7 power law for turbulent)
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --exponent 0.143

# Preview statistics without generating files
python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --preview
```

## Input Formats

### STL File
- Unit: **meters** (ensure your STL is in meters, not mm)
- Use `--scale` option if conversion needed (e.g., `--scale 0.001` for mm to m)

### CSV File

```csv
time,flowrate
0.0,0.0
0.1,5.2
0.2,4.8
```

- `time`: seconds
- `flowrate`: L/min (auto-detected) or m³/s

## Output

Generates OpenFOAM `boundaryData` format:

```
boundaryData/inlet/
├── points
├── 0.000000/U
├── 0.001000/U
└── ...
```

Place in `constant/boundaryData/inlet/` for your OpenFOAM case.

## OpenFOAM Setup

Add to `0/U`:

```cpp
inlet
{
    type            timeVaryingMappedFixedValue;
    offset          (0 0 0);
    setAverage      false;
}
```

## Limitations

### Centerline Gradient

The wall-distance power law profile U(d) = U₀ × (d/d_max)^n has a known limitation: it does not produce zero velocity gradient at the centerline for circular inlets. This is a recognised deficiency of power law formulations ([Kahine et al., 2021](https://www.mdpi.com/2311-5521/6/10/369)).

For this toolkit's intended use (patient-specific cardiovascular geometries), the priority is correct near-wall behaviour for wall shear stress calculations. For ideal circular pipes, classical radial formulations may be more appropriate.

### Turbulent Inlet Conditions

This toolkit generates **mean velocity profiles only**. For proper turbulent simulations, you also need turbulent fluctuations. Consider:

- **Synthetic turbulence methods** - `turbulentDFSEMInlet` or `turbulentDigitalFilterInlet` in OpenFOAM
- **Turbulence quantities** - estimate k and ε from turbulence intensity:
  ```cpp
  inlet
  {
      type            turbulentIntensityKineticEnergyInlet;
      intensity       0.05;  // 5% turbulence intensity
      value           uniform 0.1;
  }
  ```

For cardiovascular flows, turbulence intensity of 1-5% is typical.

## Options

| Option | Description |
|--------|-------------|
| `--profile` | `plug`, `parabolic`, `wall_distance` |
| `--flowrate` | Constant flow rate in L/min |
| `--output` | Output directory (default: `boundaryData/inlet`) |
| `--exponent` | Wall-distance exponent (default: 2.0) |
| `--scale` | STL scale factor to meters (e.g., 0.001 for mm) |
| `--points-file` | Use existing mesh face centres |
| `--flip-normal` | Flip inlet normal direction |
| `--preview` | Show statistics without generating files |
| `--verbose` | Show detailed output |

## Python API

```python
from inlet_mapper import InletMapper

mapper = InletMapper("inlet.stl", "flowrate.csv", profile="wall_distance")
mapper.generate("boundaryData/inlet")
```

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
