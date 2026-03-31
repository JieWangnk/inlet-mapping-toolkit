"""
Velocity Profile Generators
===========================

Different velocity profile shapes for inlet boundary conditions.

Available profiles:
- PlugProfile: Uniform velocity
- ParabolicProfile: Classic Poiseuille
- WallDistanceProfile: Distance-to-wall based
- WomersleyProfile: Pulsatile with frequency effects
- BluntedProfile: Flat core + parabolic wall
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Optional, Tuple


class BaseProfile(ABC):
    """Abstract base class for velocity profiles."""

    @abstractmethod
    def compute_shape_factors(
        self,
        points: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """
        Compute shape factors f(x) for each point.

        Shape factors satisfy:
        - f(wall) = 0 (no-slip)
        - f(center) = 1 (maximum)

        Parameters
        ----------
        points : np.ndarray
            Query points (N, 3).
        **kwargs : dict
            Profile-specific parameters.

        Returns
        -------
        np.ndarray
            Shape factors in [0, 1] for each point (N,).
        """
        pass

    def compute_velocities(
        self,
        points: np.ndarray,
        flow_rate: float,
        face_areas: np.ndarray,
        normal: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """
        Compute velocity vectors for given flow rate.

        Parameters
        ----------
        points : np.ndarray
            Query points (N, 3).
        flow_rate : float
            Target volumetric flow rate (m³/s).
        face_areas : np.ndarray
            Face areas for each point (N,).
        normal : np.ndarray
            Inlet normal vector (3,).
        **kwargs : dict
            Profile-specific parameters.

        Returns
        -------
        np.ndarray
            Velocity vectors (N, 3).
        """
        shape_factors = self.compute_shape_factors(points, **kwargs)

        # Compute scaling factor to match flow rate
        # Q = sum(f_i * U_0 * dA_i) => U_0 = Q / sum(f_i * dA_i)
        integrated = np.sum(shape_factors * face_areas)

        if integrated <= 0 or np.isnan(integrated):
            # Fallback to uniform if shape factors are zero or NaN
            U_0 = abs(flow_rate) / np.sum(face_areas)
            shape_factors = np.ones(len(points))
        else:
            U_0 = abs(flow_rate) / integrated

        # Velocity = U_0 * shape_factor * normal
        velocities = U_0 * np.outer(shape_factors, normal)

        return velocities


class PlugProfile(BaseProfile):
    """
    Uniform (plug) velocity profile.

    U(x) = U_mean for all points.

    Simple but unrealistic - use for testing only.
    """

    def compute_shape_factors(
        self,
        points: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """Return uniform shape factors of 1.0."""
        return np.ones(len(points))


class ParabolicProfile(BaseProfile):
    """
    Parabolic (Poiseuille) velocity profile.

    U(r) = U_max * (1 - (r/R)²)

    Classic fully-developed laminar pipe flow profile.
    Best for circular inlets with fully-developed flow.

    Parameters
    ----------
    centroid : np.ndarray
        Inlet centroid for radial distance calculation.
    radius : float
        Inlet radius (or equivalent radius for non-circular).
    """

    def __init__(self, centroid: np.ndarray, radius: float):
        self.centroid = centroid
        self.radius = radius

    def compute_shape_factors(
        self,
        points: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """
        Compute parabolic shape factors.

        f(r) = 1 - (r/R)²
        """
        radial_distances = np.linalg.norm(points - self.centroid, axis=1)
        r_normalized = np.minimum(radial_distances / self.radius, 1.0)
        return np.maximum(0, 1.0 - r_normalized ** 2)


class WallDistanceProfile(BaseProfile):
    """
    Distance-to-wall velocity profile.

    U(d) = U_0 * (d/d_max)^n

    Best for irregular inlet geometries where circular/elliptical
    assumptions don't apply.

    Parameters
    ----------
    boundary_vertices : np.ndarray
        Boundary edge vertices (M, 3).
    exponent : float
        Shape exponent (default 2.0 for quadratic).
        - n=1.0: Linear
        - n=2.0: Quadratic (parabolic-like)
        - n=1/7: Turbulent power law
    """

    def __init__(
        self,
        boundary_vertices: np.ndarray,
        exponent: float = 2.0
    ):
        self.boundary_vertices = boundary_vertices
        self.exponent = exponent
        self._tree = None

    def _get_tree(self):
        """Lazy initialization of KDTree."""
        if self._tree is None:
            try:
                from scipy.spatial import cKDTree
                self._tree = cKDTree(self.boundary_vertices)
            except ImportError:
                pass
        return self._tree

    def compute_shape_factors(
        self,
        points: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """
        Compute wall-distance shape factors.

        f(d) = (d/d_max)^n
        """
        # Calculate distances to boundary
        tree = self._get_tree()
        if tree is not None:
            distances, _ = tree.query(points, k=1)
        else:
            # Fallback
            distances = np.array([
                np.min(np.linalg.norm(self.boundary_vertices - pt, axis=1))
                for pt in points
            ])

        # Normalize by maximum distance
        d_max = np.max(distances)
        if d_max <= 0:
            return np.ones(len(points))

        d_normalized = np.minimum(distances / d_max, 1.0)

        # Apply power law
        return d_normalized ** self.exponent


class WomersleyProfile(BaseProfile):
    """
    Womersley velocity profile for pulsatile flow.

    Analytical solution for oscillatory flow in a circular tube.
    Captures frequency-dependent effects where higher frequencies
    produce flatter profiles.

    Parameters
    ----------
    centroid : np.ndarray
        Inlet centroid.
    radius : float
        Inlet radius.
    kinematic_viscosity : float
        Kinematic viscosity (m²/s). Default 3.5e-6 for blood.
    n_harmonics : int
        Number of Fourier harmonics to include.
    """

    def __init__(
        self,
        centroid: np.ndarray,
        radius: float,
        kinematic_viscosity: float = 3.5e-6,
        n_harmonics: int = 10
    ):
        self.centroid = centroid
        self.radius = radius
        self.nu = kinematic_viscosity
        self.n_harmonics = n_harmonics

    def compute_shape_factors(
        self,
        points: np.ndarray,
        time: float = 0.0,
        period: float = 1.0,
        fourier_coeffs: Optional[np.ndarray] = None,
        **kwargs
    ) -> np.ndarray:
        """
        Compute Womersley shape factors.

        For steady component, returns parabolic profile.
        For oscillating components, applies Womersley solution.

        Parameters
        ----------
        time : float
            Current time (s).
        period : float
            Cardiac period (s).
        fourier_coeffs : np.ndarray, optional
            Complex Fourier coefficients of flow waveform.
            If None, uses parabolic profile.
        """
        radial_distances = np.linalg.norm(points - self.centroid, axis=1)
        r_normalized = np.minimum(radial_distances / self.radius, 1.0)

        if fourier_coeffs is None:
            # Fallback to parabolic
            return np.maximum(0, 1.0 - r_normalized ** 2)

        # Angular frequency
        omega = 2 * np.pi / period

        # Initialize with steady component (parabolic)
        shape = 2.0 * (1.0 - r_normalized ** 2)

        # Add oscillating components
        try:
            from scipy.special import jv  # Bessel function

            for k in range(1, min(len(fourier_coeffs), self.n_harmonics + 1)):
                if abs(fourier_coeffs[k]) < 1e-10:
                    continue

                # Womersley number for this harmonic
                alpha_k = self.radius * np.sqrt(k * omega / self.nu)

                if alpha_k < 0.1:
                    # Low Womersley: quasi-steady parabolic
                    wom_factor = 1.0 - r_normalized ** 2
                else:
                    # Full Womersley solution
                    i32 = complex(1, 1) / np.sqrt(2) * alpha_k

                    # J0(i^(3/2) * alpha * r/R) / J0(i^(3/2) * alpha)
                    bessel_ratio = (
                        jv(0, i32 * r_normalized) /
                        jv(0, i32)
                    )
                    wom_factor = np.real(1.0 - bessel_ratio)

                # Add harmonic contribution
                phase = k * omega * time
                amplitude = abs(fourier_coeffs[k])
                shape += amplitude * wom_factor * np.cos(phase + np.angle(fourier_coeffs[k]))

        except ImportError:
            print("Warning: scipy not available for Womersley profile, using parabolic")

        # Normalize to [0, 1]
        shape_max = np.max(shape)
        if shape_max > 0:
            shape = shape / shape_max

        return np.maximum(0, shape)

    def compute_fourier_coefficients(
        self,
        times: np.ndarray,
        flow_rates: np.ndarray
    ) -> np.ndarray:
        """
        Compute Fourier coefficients from flow waveform.

        Parameters
        ----------
        times : np.ndarray
            Time points.
        flow_rates : np.ndarray
            Flow rate values.

        Returns
        -------
        np.ndarray
            Complex Fourier coefficients.
        """
        n = len(times)
        coeffs = np.fft.fft(flow_rates) / n
        return coeffs[:self.n_harmonics + 1]


class BluntedProfile(BaseProfile):
    """
    Blunted velocity profile with flat core.

    Flat (plug) velocity in the core region, parabolic decay near walls.
    Models transitional or mildly turbulent inlet conditions.

    Parameters
    ----------
    boundary_vertices : np.ndarray
        Boundary edge vertices.
    core_fraction : float
        Fraction of d_max where core starts (default 0.3).
        Points with d > core_fraction * d_max have uniform velocity.
    """

    def __init__(
        self,
        boundary_vertices: np.ndarray,
        core_fraction: float = 0.3
    ):
        self.boundary_vertices = boundary_vertices
        self.core_fraction = core_fraction
        self._tree = None

    def _get_tree(self):
        """Lazy initialization of KDTree."""
        if self._tree is None:
            try:
                from scipy.spatial import cKDTree
                self._tree = cKDTree(self.boundary_vertices)
            except ImportError:
                pass
        return self._tree

    def compute_shape_factors(
        self,
        points: np.ndarray,
        **kwargs
    ) -> np.ndarray:
        """
        Compute blunted shape factors.

        f(d) = 1.0                      for d >= core_dist
        f(d) = 1 - (1 - d/core_dist)²   for d < core_dist
        """
        # Calculate distances to boundary
        tree = self._get_tree()
        if tree is not None:
            distances, _ = tree.query(points, k=1)
        else:
            distances = np.array([
                np.min(np.linalg.norm(self.boundary_vertices - pt, axis=1))
                for pt in points
            ])

        d_max = np.max(distances)
        if d_max <= 0:
            return np.ones(len(points))

        core_dist = self.core_fraction * d_max

        # Initialize all to 1.0 (core value)
        shape = np.ones(len(points))

        # Apply parabolic decay near walls
        near_wall = distances < core_dist
        if np.any(near_wall):
            local_r = 1.0 - (distances[near_wall] / core_dist)
            shape[near_wall] = np.maximum(0, 1.0 - local_r ** 2)

        return shape


def get_profile(
    profile_type: str,
    geometry,
    **kwargs
) -> BaseProfile:
    """
    Factory function to create profile instance.

    Parameters
    ----------
    profile_type : str
        Profile type: 'plug', 'parabolic', 'wall_distance', 'womersley', 'blunted'
    geometry : InletGeometry
        Inlet geometry object.
    **kwargs : dict
        Profile-specific parameters.

    Returns
    -------
    BaseProfile
        Profile instance.
    """
    profile_type = profile_type.lower().strip()

    if profile_type == 'plug':
        return PlugProfile()

    elif profile_type in ('parabolic', 'poiseuille'):
        return ParabolicProfile(
            centroid=geometry.centroid,
            radius=geometry.equivalent_radius
        )

    elif profile_type in ('wall_distance', 'distance_wall', 'wall-distance'):
        return WallDistanceProfile(
            boundary_vertices=geometry.boundary_vertices,
            exponent=kwargs.get('exponent', 2.0)
        )

    elif profile_type == 'womersley':
        return WomersleyProfile(
            centroid=geometry.centroid,
            radius=geometry.equivalent_radius,
            kinematic_viscosity=kwargs.get('viscosity', 3.5e-6),
            n_harmonics=kwargs.get('n_harmonics', 10)
        )

    elif profile_type in ('blunted', 'blunted_parabolic'):
        return BluntedProfile(
            boundary_vertices=geometry.boundary_vertices,
            core_fraction=kwargs.get('core_fraction', 0.3)
        )

    else:
        raise ValueError(
            f"Unknown profile type: {profile_type}. "
            f"Available: plug, parabolic, wall_distance, womersley, blunted"
        )
