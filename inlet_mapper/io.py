"""
Input/Output Utilities
======================

CSV reading and OpenFOAM boundaryData writing.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List
import os


class CSVReader:
    """
    Read flow rate data from CSV files.

    Supports various CSV formats:
    - time,flowrate (L/min or m³/s)
    - time,velocity (m/s)

    Auto-detects units based on magnitude.

    Parameters
    ----------
    csv_path : str or Path
        Path to CSV file.
    """

    def __init__(self, csv_path: str):
        self.csv_path = Path(csv_path)

        if not self.csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")

        self._load_data()

    def _load_data(self) -> None:
        """Load and parse CSV data."""
        # Check for header
        with open(self.csv_path, 'r') as f:
            first_line = f.readline()
            has_header = any(c.isalpha() for c in first_line)

        skiprows = 1 if has_header else 0

        # Load data
        data = np.genfromtxt(
            self.csv_path,
            delimiter=',',
            skip_header=skiprows
        )

        if data.ndim < 2 or data.shape[1] < 2:
            raise ValueError("CSV must have at least 2 columns (time, value)")

        self.times = data[:, 0]
        self.values = data[:, 1]

        # Determine data type and units
        self._detect_units(has_header, first_line if has_header else None)

    def _detect_units(self, has_header: bool, header_line: Optional[str]) -> None:
        """Auto-detect data type and units."""
        # Default assumptions
        self.data_type = 'flowrate'
        self.units = 'L/min'

        # Check header for clues
        if header_line:
            header_lower = header_line.lower()
            if 'velocity' in header_lower or 'm/s' in header_lower:
                self.data_type = 'velocity'
                self.units = 'm/s'
            elif 'm3/s' in header_lower or 'm³/s' in header_lower:
                self.units = 'm3/s'

        # Auto-detect based on magnitude
        max_val = np.max(np.abs(self.values))
        if self.data_type == 'flowrate':
            if max_val < 0.01:
                # Likely already in m³/s
                self.units = 'm3/s'
            else:
                # Likely in L/min
                self.units = 'L/min'

    def get_flow_rates_m3s(self, inlet_area: Optional[float] = None) -> np.ndarray:
        """
        Get flow rates in m³/s.

        Parameters
        ----------
        inlet_area : float, optional
            Inlet area (m²), required if data_type is 'velocity'.

        Returns
        -------
        np.ndarray
            Flow rates in m³/s.
        """
        if self.data_type == 'velocity':
            if inlet_area is None:
                raise ValueError("inlet_area required for velocity data")
            return self.values * inlet_area

        # Flow rate data
        if self.units == 'L/min':
            return self.values * 1e-3 / 60.0  # L/min -> m³/s
        else:
            return self.values  # Already m³/s

    def get_cardiac_period(self) -> float:
        """
        Get cardiac period (cycle duration).

        Returns
        -------
        float
            Period in seconds.
        """
        return self.times[-1] - self.times[0]

    def summary(self) -> dict:
        """Get data summary."""
        flow_rates = self.get_flow_rates_m3s() if self.data_type == 'flowrate' else None

        return {
            "csv_path": str(self.csv_path),
            "n_timesteps": len(self.times),
            "time_range": [self.times[0], self.times[-1]],
            "data_type": self.data_type,
            "units": self.units,
            "value_range": [np.min(self.values), np.max(self.values)],
            "cardiac_period": self.get_cardiac_period(),
        }


class BoundaryDataWriter:
    """
    Write OpenFOAM boundaryData format.

    Creates directory structure:
    ```
    output_dir/
    ├── points          # Face centre coordinates
    ├── 0.000000/
    │   └── U           # Velocity at t=0
    ├── 0.001000/
    │   └── U           # Velocity at t=0.001
    └── ...
    ```

    Parameters
    ----------
    output_dir : str or Path
        Output directory path.
    """

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def write_points(self, points: np.ndarray) -> None:
        """
        Write points file.

        Parameters
        ----------
        points : np.ndarray
            Point coordinates (N, 3).
        """
        filepath = self.output_dir / "points"
        n = len(points)

        with open(filepath, 'w') as f:
            f.write(f"{n}\n(\n")
            for pt in points:
                f.write(f"({pt[0]:.10e} {pt[1]:.10e} {pt[2]:.10e})\n")
            f.write(")\n")

    def write_velocity(
        self,
        time: float,
        velocities: np.ndarray
    ) -> None:
        """
        Write velocity file for a single timestep.

        Parameters
        ----------
        time : float
            Time value.
        velocities : np.ndarray
            Velocity vectors (N, 3).
        """
        time_dir = self.output_dir / f"{time:.6f}"
        time_dir.mkdir(exist_ok=True)

        filepath = time_dir / "U"
        n = len(velocities)

        # Fast writing using join
        lines = [f"{n}", "("]
        lines.extend(
            f"({v[0]:.6e} {v[1]:.6e} {v[2]:.6e})"
            for v in velocities
        )
        lines.append(")")

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines) + '\n')

    def write_all_velocities(
        self,
        times: np.ndarray,
        velocities_list: List[np.ndarray],
        show_progress: bool = True
    ) -> None:
        """
        Write velocity files for all timesteps.

        Parameters
        ----------
        times : np.ndarray
            Time values.
        velocities_list : list of np.ndarray
            Velocity vectors for each timestep.
        show_progress : bool
            Print progress updates.
        """
        n_times = len(times)

        for i, (t, vel) in enumerate(zip(times, velocities_list)):
            self.write_velocity(t, vel)

            if show_progress and i % max(1, n_times // 10) == 0:
                print(f"  Writing timestep {i+1}/{n_times} (t={t:.4f}s)")

    def clean_old_timesteps(self) -> None:
        """Remove old timestep directories."""
        for item in self.output_dir.iterdir():
            if item.is_dir():
                # Check if directory name is numeric (timestep)
                try:
                    float(item.name)
                    import shutil
                    shutil.rmtree(item)
                except ValueError:
                    pass


class PointsFileReader:
    """
    Read existing OpenFOAM points file.

    Use this when mapping to an existing mesh's face centres
    rather than STL face centres.
    """

    def __init__(self, points_path: str):
        self.points_path = Path(points_path)

        if not self.points_path.exists():
            raise FileNotFoundError(f"Points file not found: {points_path}")

        self._load_points()

    def _load_points(self) -> None:
        """Parse OpenFOAM points file format."""
        import re

        with open(self.points_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        # First line is count
        self.n_points = int(lines[0])

        # Find start of data
        start_idx = 2 if lines[1] == '(' else 1

        # Parse points
        points = []
        for i in range(self.n_points):
            line = re.sub(r'[()]', '', lines[i + start_idx])
            xyz = [float(p) for p in line.split()]
            points.append(xyz)

        self.points = np.array(points, dtype=float)

    def get_points(self) -> np.ndarray:
        """Get point coordinates."""
        return self.points
