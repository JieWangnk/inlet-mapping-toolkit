"""
Inlet Geometry Processing
=========================

Handles STL loading, boundary extraction, and geometric calculations.
"""

import numpy as np
from typing import Tuple, Optional, List
from collections import Counter
from pathlib import Path


class InletGeometry:
    """
    Process inlet STL geometry for velocity profile generation.

    Extracts:
    - Centroid and normal vector
    - Surface area
    - Boundary edge vertices (for wall-distance profiles)
    - Face centres (sampling points)

    Parameters
    ----------
    stl_path : str or Path
        Path to inlet STL file.
    scale_factor : float, optional
        Scale factor to convert STL units to meters. Default 1.0 (already in meters).

    Example
    -------
    >>> geom = InletGeometry("inlet.stl")
    >>> print(f"Area: {geom.area:.6e} m²")
    >>> print(f"Normal: {geom.normal}")
    """

    def __init__(self, stl_path: str, scale_factor: float = 1.0):
        self.stl_path = Path(stl_path)
        self.scale_factor = scale_factor

        if not self.stl_path.exists():
            raise FileNotFoundError(f"STL file not found: {stl_path}")

        # Load mesh
        self._load_mesh()

        # Calculate properties
        self._calculate_properties()

    def _load_mesh(self) -> None:
        """Load STL mesh using trimesh."""
        try:
            import trimesh
        except ImportError:
            raise ImportError(
                "trimesh is required for STL processing. "
                "Install with: pip install trimesh"
            )

        self.mesh = trimesh.load(self.stl_path)

        # Apply scale factor
        if self.scale_factor != 1.0:
            self.mesh.apply_scale(self.scale_factor)

    def _calculate_properties(self) -> None:
        """Calculate geometric properties."""
        # Centroid
        self.centroid = np.mean(self.mesh.vertices, axis=0)

        # Area (sum of triangle areas)
        self.area = self.mesh.area

        # Normal (area-weighted average of face normals)
        face_areas = self.mesh.area_faces
        face_normals = self.mesh.face_normals
        weighted_normal = np.sum(face_normals * face_areas[:, np.newaxis], axis=0)
        self.normal = weighted_normal / np.linalg.norm(weighted_normal)

        # Equivalent radius (circular approximation)
        self.equivalent_radius = np.sqrt(self.area / np.pi)

        # Extract boundary vertices
        self.boundary_vertices = self._extract_boundary()

    def _extract_boundary(self) -> np.ndarray:
        """
        Extract boundary edge vertices from the mesh.

        Boundary edges are edges that belong to only one triangle face
        (i.e., they form the open edge of the inlet surface).

        Returns
        -------
        np.ndarray
            Coordinates of boundary vertices (N, 3).
        """
        # Find all edges
        edges = []
        for face in self.mesh.faces:
            edges.append(tuple(sorted([face[0], face[1]])))
            edges.append(tuple(sorted([face[1], face[2]])))
            edges.append(tuple(sorted([face[2], face[0]])))

        # Count edge occurrences
        edge_counts = Counter(edges)

        # Boundary edges appear only once
        boundary_edges = [e for e, count in edge_counts.items() if count == 1]

        if not boundary_edges:
            # Fallback: use peripheral vertices (mesh may be closed)
            print("Warning: No boundary edges found, using peripheral vertices")
            distances = np.linalg.norm(self.mesh.vertices - self.centroid, axis=1)
            threshold = 0.85 * np.max(distances)
            return self.mesh.vertices[distances > threshold]

        # Get unique boundary vertex indices
        boundary_indices = set()
        for e in boundary_edges:
            boundary_indices.add(e[0])
            boundary_indices.add(e[1])

        return self.mesh.vertices[list(boundary_indices)]

    def get_face_centres(self) -> np.ndarray:
        """
        Get face centre coordinates.

        Returns
        -------
        np.ndarray
            Face centre coordinates (N, 3).
        """
        # Each face centre is the average of its 3 vertices
        vertices = self.mesh.vertices
        faces = self.mesh.faces

        centres = (
            vertices[faces[:, 0]] +
            vertices[faces[:, 1]] +
            vertices[faces[:, 2]]
        ) / 3.0

        return centres

    def resample_surface(self, n_points: int = 200) -> np.ndarray:
        """
        Generate well-distributed sample points across the inlet surface.

        Uses trimesh to sample points uniformly over the surface area,
        ensuring good radial distribution even for coarse STL meshes.

        Parameters
        ----------
        n_points : int
            Number of sample points to generate.

        Returns
        -------
        np.ndarray
            Resampled point coordinates (n_points, 3).
        """
        points, _ = self.mesh.sample(n_points, return_index=True)
        return points

    def get_face_areas(self) -> np.ndarray:
        """
        Get individual face areas.

        Returns
        -------
        np.ndarray
            Face areas (N,).
        """
        return self.mesh.area_faces

    def compute_wall_distances(self, points: np.ndarray) -> np.ndarray:
        """
        Compute distance from each point to the nearest boundary vertex.

        Uses KDTree for O(N log M) performance.

        Parameters
        ----------
        points : np.ndarray
            Query points (N, 3).

        Returns
        -------
        np.ndarray
            Distance to nearest boundary for each point (N,).
        """
        try:
            from scipy.spatial import cKDTree
            tree = cKDTree(self.boundary_vertices)
            distances, _ = tree.query(points, k=1)
        except ImportError:
            # Fallback: vectorized numpy (slower)
            print("Warning: scipy not available, using slower distance calculation")
            distances = np.array([
                np.min(np.linalg.norm(self.boundary_vertices - pt, axis=1))
                for pt in points
            ])

        return distances

    def compute_radial_distances(self, points: np.ndarray) -> np.ndarray:
        """
        Compute radial distance from centroid for each point.

        Parameters
        ----------
        points : np.ndarray
            Query points (N, 3).

        Returns
        -------
        np.ndarray
            Radial distances (N,).
        """
        return np.linalg.norm(points - self.centroid, axis=1)

    def flip_normal(self) -> None:
        """Flip the normal direction."""
        self.normal = -self.normal

    def should_flip_normal(self, reference_point: Optional[np.ndarray] = None) -> bool:
        """
        Determine if normal should be flipped based on reference point.

        If reference_point is downstream of the inlet (flow direction),
        the normal should point toward it.

        Parameters
        ----------
        reference_point : np.ndarray, optional
            Reference point (e.g., outlet centroid).

        Returns
        -------
        bool
            True if normal should be flipped.
        """
        if reference_point is None:
            return False

        flow_direction = reference_point - self.centroid
        flow_direction = flow_direction / np.linalg.norm(flow_direction)

        # If normal points opposite to flow direction, flip it
        return np.dot(self.normal, flow_direction) < 0

    def summary(self) -> dict:
        """
        Get geometry summary.

        Returns
        -------
        dict
            Geometry properties.
        """
        return {
            "stl_path": str(self.stl_path),
            "centroid": self.centroid.tolist(),
            "normal": self.normal.tolist(),
            "area_m2": self.area,
            "equivalent_radius_m": self.equivalent_radius,
            "n_faces": len(self.mesh.faces),
            "n_vertices": len(self.mesh.vertices),
            "n_boundary_vertices": len(self.boundary_vertices),
        }

    def __repr__(self) -> str:
        return (
            f"InletGeometry(\n"
            f"  path={self.stl_path.name},\n"
            f"  area={self.area:.6e} m²,\n"
            f"  radius={self.equivalent_radius*1000:.2f} mm,\n"
            f"  normal={self.normal}\n"
            f")"
        )
