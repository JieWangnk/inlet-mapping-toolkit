"""
Core Inlet Mapper
=================

Main class that orchestrates the inlet mapping process.
"""

import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Any

from .geometry import InletGeometry
from .profiles import get_profile, BaseProfile, WomersleyProfile
from .io import CSVReader, BoundaryDataWriter, PointsFileReader


class InletMapper:
    """
    Map 1D flow rate data to 3D velocity profiles for OpenFOAM.

    This is the main interface for the inlet mapping toolkit.
    It coordinates geometry loading, profile calculation, and output writing.

    Parameters
    ----------
    stl_path : str
        Path to inlet STL file.
    csv_path : str, optional
        Path to flow rate CSV file. Required for time-varying profiles.
    profile : str
        Profile type: 'plug', 'parabolic', 'wall_distance', 'womersley', 'blunted'
    flowrate : float, optional
        Constant flow rate (L/min) if csv_path not provided.
    scale_factor : float
        STL scale factor to convert to meters. Default 1.0.
    flip_normal : bool
        Flip inlet normal direction.
    **kwargs : dict
        Profile-specific parameters (exponent, viscosity, etc.)

    Example
    -------
    >>> mapper = InletMapper(
    ...     stl_path="inlet.stl",
    ...     csv_path="flowrate.csv",
    ...     profile="wall_distance",
    ...     exponent=2.0
    ... )
    >>> mapper.generate("boundaryData/inlet")
    """

    def __init__(
        self,
        stl_path: str,
        csv_path: Optional[str] = None,
        profile: str = "parabolic",
        flowrate: Optional[float] = None,
        scale_factor: float = 1.0,
        flip_normal: bool = False,
        points_file: Optional[str] = None,
        **kwargs
    ):
        self.stl_path = stl_path
        self.csv_path = csv_path
        self.profile_type = profile
        self.constant_flowrate = flowrate
        self.scale_factor = scale_factor
        self.flip_normal = flip_normal
        self.points_file = points_file
        self.profile_kwargs = kwargs

        # Validate inputs
        if csv_path is None and flowrate is None:
            raise ValueError("Either csv_path or flowrate must be provided")

        # Load geometry
        print(f"Loading inlet geometry from: {stl_path}")
        self.geometry = InletGeometry(stl_path, scale_factor)

        if flip_normal:
            self.geometry.flip_normal()
            print("  Normal direction flipped")

        print(self.geometry)

        # Load flow data if provided
        self.csv_data = None
        if csv_path:
            print(f"\nLoading flow data from: {csv_path}")
            self.csv_data = CSVReader(csv_path)
            print(f"  {self.csv_data.summary()}")

        # Load existing points file if provided
        self.mesh_points = None
        if points_file:
            print(f"\nLoading mesh points from: {points_file}")
            reader = PointsFileReader(points_file)
            self.mesh_points = reader.get_points()
            print(f"  Loaded {len(self.mesh_points)} face centres")

        # Create profile generator
        self.profile = get_profile(
            profile,
            self.geometry,
            **kwargs
        )
        print(f"\nProfile type: {profile}")

    def generate(
        self,
        output_dir: str,
        clean_existing: bool = True,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Generate boundaryData output.

        Parameters
        ----------
        output_dir : str
            Output directory path.
        clean_existing : bool
            Remove existing timestep directories.
        verbose : bool
            Print progress messages.

        Returns
        -------
        dict
            Generation statistics.
        """
        print(f"\n{'='*60}")
        print("Generating Inlet Boundary Data")
        print(f"{'='*60}")

        # Setup writer
        writer = BoundaryDataWriter(output_dir)

        if clean_existing:
            writer.clean_old_timesteps()
            if verbose:
                print("Cleaned existing timestep directories")

        # Get sampling points (mesh points or STL face centres)
        if self.mesh_points is not None:
            points = self.mesh_points
            # Estimate face areas (uniform)
            face_areas = np.full(len(points), self.geometry.area / len(points))
            if verbose:
                print(f"Using {len(points)} mesh face centres")
        else:
            points = self.geometry.get_face_centres()
            face_areas = self.geometry.get_face_areas()
            if verbose:
                print(f"Using {len(points)} STL face centres")

        # Write points file
        writer.write_points(points)
        if verbose:
            print(f"Wrote points file ({len(points)} points)")

        # Generate velocity data
        if self.csv_data:
            stats = self._generate_time_varying(
                writer, points, face_areas, verbose
            )
        else:
            stats = self._generate_constant(
                writer, points, face_areas, verbose
            )

        print(f"\n{'='*60}")
        print("Generation Complete")
        print(f"{'='*60}")
        print(f"Output directory: {output_dir}")

        return stats

    def _generate_constant(
        self,
        writer: BoundaryDataWriter,
        points: np.ndarray,
        face_areas: np.ndarray,
        verbose: bool
    ) -> Dict[str, Any]:
        """Generate constant flow profile (single timestep at t=0)."""
        # Convert flowrate from L/min to m³/s
        flow_rate_m3s = self.constant_flowrate * 1e-3 / 60.0

        if verbose:
            print(f"\nConstant flow rate: {self.constant_flowrate:.2f} L/min "
                  f"({flow_rate_m3s:.6e} m³/s)")

        # Compute velocities
        velocities = self.profile.compute_velocities(
            points,
            flow_rate_m3s,
            face_areas,
            self.geometry.normal,
            **self.profile_kwargs
        )

        # Write single timestep
        writer.write_velocity(0.0, velocities)

        # Statistics
        magnitudes = np.linalg.norm(velocities, axis=1)
        stats = {
            "n_timesteps": 1,
            "n_points": len(points),
            "flow_rate_m3s": flow_rate_m3s,
            "min_velocity": float(np.min(magnitudes)),
            "max_velocity": float(np.max(magnitudes)),
            "mean_velocity": float(np.mean(magnitudes)),
        }

        if verbose:
            print(f"\nVelocity statistics:")
            print(f"  Min:  {stats['min_velocity']:.4f} m/s")
            print(f"  Max:  {stats['max_velocity']:.4f} m/s")
            print(f"  Mean: {stats['mean_velocity']:.4f} m/s")

        return stats

    def _generate_time_varying(
        self,
        writer: BoundaryDataWriter,
        points: np.ndarray,
        face_areas: np.ndarray,
        verbose: bool
    ) -> Dict[str, Any]:
        """Generate time-varying velocity profiles."""
        times = self.csv_data.times
        flow_rates = self.csv_data.get_flow_rates_m3s(self.geometry.area)
        period = self.csv_data.get_cardiac_period()

        if verbose:
            print(f"\nTime-varying data:")
            print(f"  Timesteps: {len(times)}")
            print(f"  Period: {period:.4f} s")
            print(f"  Flow rate range: [{np.min(flow_rates):.6e}, "
                  f"{np.max(flow_rates):.6e}] m³/s")

        # For Womersley profile, compute Fourier coefficients
        fourier_coeffs = None
        if isinstance(self.profile, WomersleyProfile):
            fourier_coeffs = self.profile.compute_fourier_coefficients(
                times, flow_rates
            )
            if verbose:
                print(f"  Computed {len(fourier_coeffs)} Fourier harmonics")

        # Pre-compute shape factors for non-Womersley profiles
        # (they don't depend on time)
        if not isinstance(self.profile, WomersleyProfile):
            shape_factors = self.profile.compute_shape_factors(
                points, **self.profile_kwargs
            )
            integrated = np.sum(shape_factors * face_areas)
            velocity_directions = np.outer(shape_factors, self.geometry.normal)

        # Generate velocity for each timestep
        all_magnitudes = []

        if verbose:
            print(f"\nGenerating {len(times)} timesteps...")

        for i, t in enumerate(times):
            Q = flow_rates[i]

            if isinstance(self.profile, WomersleyProfile):
                # Womersley needs time-dependent shape factors
                velocities = self.profile.compute_velocities(
                    points, Q, face_areas, self.geometry.normal,
                    time=t, period=period, fourier_coeffs=fourier_coeffs,
                    **self.profile_kwargs
                )
            else:
                # Other profiles: just scale pre-computed shape factors
                U_0 = abs(Q) / integrated if integrated > 0 else abs(Q) / self.geometry.area
                velocities = U_0 * velocity_directions

            writer.write_velocity(t, velocities)

            magnitudes = np.linalg.norm(velocities, axis=1)
            all_magnitudes.append(magnitudes)

            # Progress update
            if verbose and i % max(1, len(times) // 10) == 0:
                print(f"  t={t:.4f}s: Q={Q:.4e} m³/s, "
                      f"U_max={np.max(magnitudes):.4f} m/s")

        # Statistics
        all_magnitudes = np.array(all_magnitudes)
        stats = {
            "n_timesteps": len(times),
            "n_points": len(points),
            "time_range": [float(times[0]), float(times[-1])],
            "period": period,
            "flow_rate_range": [float(np.min(flow_rates)), float(np.max(flow_rates))],
            "min_velocity": float(np.min(all_magnitudes)),
            "max_velocity": float(np.max(all_magnitudes)),
            "mean_velocity": float(np.mean(all_magnitudes)),
        }

        if verbose:
            print(f"\nOverall velocity statistics:")
            print(f"  Min:  {stats['min_velocity']:.4f} m/s")
            print(f"  Max:  {stats['max_velocity']:.4f} m/s")
            print(f"  Mean: {stats['mean_velocity']:.4f} m/s")

        return stats

    def get_statistics(self) -> Dict[str, Any]:
        """
        Get geometry and data statistics without generating output.

        Returns
        -------
        dict
            Statistics dictionary.
        """
        stats = {
            "geometry": self.geometry.summary(),
            "profile": self.profile_type,
        }

        if self.csv_data:
            stats["csv_data"] = self.csv_data.summary()

        if self.constant_flowrate:
            stats["constant_flowrate_Lmin"] = self.constant_flowrate
            stats["constant_flowrate_m3s"] = self.constant_flowrate * 1e-3 / 60.0

        return stats

    def preview(
        self,
        n_points: int = 100,
        time: float = 0.0
    ) -> Dict[str, np.ndarray]:
        """
        Preview velocity profile without writing files.

        Parameters
        ----------
        n_points : int
            Number of preview points.
        time : float
            Time for preview (for time-varying profiles).

        Returns
        -------
        dict
            Preview data with 'points', 'velocities', 'magnitudes'.
        """
        # Sample points
        if self.mesh_points is not None:
            indices = np.linspace(0, len(self.mesh_points)-1, n_points, dtype=int)
            points = self.mesh_points[indices]
            face_areas = np.full(n_points, self.geometry.area / len(self.mesh_points))
        else:
            all_points = self.geometry.get_face_centres()
            all_areas = self.geometry.get_face_areas()
            indices = np.linspace(0, len(all_points)-1, n_points, dtype=int)
            points = all_points[indices]
            face_areas = all_areas[indices]

        # Get flow rate
        if self.csv_data:
            # Interpolate flow rate at specified time
            flow_rate = np.interp(
                time,
                self.csv_data.times,
                self.csv_data.get_flow_rates_m3s(self.geometry.area)
            )
        else:
            flow_rate = self.constant_flowrate * 1e-3 / 60.0

        # Compute velocities
        velocities = self.profile.compute_velocities(
            points, flow_rate, face_areas, self.geometry.normal,
            **self.profile_kwargs
        )

        return {
            "points": points,
            "velocities": velocities,
            "magnitudes": np.linalg.norm(velocities, axis=1),
            "flow_rate_m3s": flow_rate,
            "time": time,
        }
