# Examples

## Example Files

- `example_flowrate.csv` - Sample cardiac flow waveform (1 second cycle)

## Usage Examples

### 1. Basic Wall-Distance Profile

```bash
# Using the example flow waveform
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile wall_distance \
    --output boundaryData/inlet
```

### 2. Constant Flow with Parabolic Profile

```bash
# 5 L/min constant flow
python map_inlet.py your_inlet.stl \
    --flowrate 5.0 \
    --profile parabolic \
    --output boundaryData/inlet
```

### 3. Womersley Profile for Pulsatile Flow

```bash
# Full Womersley solution with 10 harmonics
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile womersley \
    --viscosity 3.5e-6 \
    --n-harmonics 10 \
    --output boundaryData/inlet
```

### 4. Custom Wall-Distance Exponent

```bash
# Linear (n=1) for turbulent-like profile
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile wall_distance \
    --exponent 1.0 \
    --output boundaryData/inlet

# 1/7 power law for turbulent pipe flow
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile wall_distance \
    --exponent 0.143 \
    --output boundaryData/inlet
```

### 5. Blunted Profile

```bash
# 30% flat core
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile blunted \
    --core-fraction 0.3 \
    --output boundaryData/inlet
```

### 6. Preview Without Generating Files

```bash
# Show statistics only
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile wall_distance \
    --preview
```

### 7. Using Existing Mesh Points

```bash
# Map to existing OpenFOAM mesh face centres
python map_inlet.py your_inlet.stl examples/example_flowrate.csv \
    --profile wall_distance \
    --points-file case/constant/boundaryData/inlet/points \
    --output case/constant/boundaryData/inlet
```

## Output Structure

After running, you'll have:

```
boundaryData/inlet/
├── points              # Face centre coordinates
├── 0.000000/
│   └── U               # Velocity at t=0.0s
├── 0.050000/
│   └── U               # Velocity at t=0.05s
├── 0.100000/
│   └── U               # Velocity at t=0.1s
└── ...
```

## OpenFOAM Boundary Condition

In your `0/U` file:

```cpp
inlet
{
    type            timeVaryingMappedFixedValue;
    offset          (0 0 0);
    setAverage      false;
}
```

The `boundaryData/inlet/` directory should be placed in `constant/boundaryData/inlet/`.
