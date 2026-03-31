"""
Tests for inlet mapping profiles and coordinate mismatch handling.

Covers the bug where parabolic (and other) profiles produced all-zero
velocities when --points-file coordinates didn't match the STL geometry.
"""

import numpy as np
import pytest
import tempfile
import shutil
import os

from inlet_mapper.profiles import (
    BaseProfile,
    PlugProfile,
    ParabolicProfile,
    WallDistanceProfile,
    BluntedProfile,
)


# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

class ZeroProfile(BaseProfile):
    """Profile that always returns zero shape factors."""
    def compute_shape_factors(self, points, **kw):
        return np.zeros(len(points))


class NaNProfile(BaseProfile):
    """Profile that always returns NaN shape factors."""
    def compute_shape_factors(self, points, **kw):
        return np.full(len(points), np.nan)


SAMPLE_POINTS = np.array([[0, 0, 0], [0.003, 0, 0], [0, 0.003, 0]], dtype=float)
SAMPLE_AREAS = np.array([1e-5, 1e-5, 1e-5])
SAMPLE_NORMAL = np.array([0.0, 0.0, -1.0])


# ----------------------------------------------------------------
# Unit tests: individual profile shape factors
# ----------------------------------------------------------------

class TestPlugProfile:
    def test_always_ones(self):
        sf = PlugProfile().compute_shape_factors(SAMPLE_POINTS)
        np.testing.assert_array_equal(sf, np.ones(3))


class TestParabolicProfile:
    def test_matching_coordinates(self):
        """Points within radius should produce non-zero shape factors."""
        profile = ParabolicProfile(centroid=np.array([0, 0, 0]), radius=0.01)
        sf = profile.compute_shape_factors(SAMPLE_POINTS)
        assert np.all(sf > 0), f"Expected all > 0, got {sf}"

    def test_centre_is_maximum(self):
        """Point at centroid should have shape factor = 1."""
        profile = ParabolicProfile(centroid=np.array([0, 0, 0]), radius=0.01)
        sf = profile.compute_shape_factors(np.array([[0, 0, 0]]))
        assert sf[0] == pytest.approx(1.0)

    def test_boundary_is_zero(self):
        """Point at r = R should have shape factor = 0."""
        profile = ParabolicProfile(centroid=np.array([0, 0, 0]), radius=0.01)
        sf = profile.compute_shape_factors(np.array([[0.01, 0, 0]]))
        assert sf[0] == pytest.approx(0.0)

    def test_coordinate_mismatch_produces_zeros(self):
        """Offset centroid causes all shape factors to be zero."""
        profile = ParabolicProfile(centroid=np.array([100, 200, 300]), radius=0.01)
        sf = profile.compute_shape_factors(SAMPLE_POINTS)
        assert np.all(sf == 0)


class TestWallDistanceProfile:
    def test_shape_variation(self):
        """Wall-distance profile should produce varying shape factors."""
        boundary = np.array([
            [0.005, 0, 0], [-0.005, 0, 0],
            [0, 0.005, 0], [0, -0.005, 0],
        ])
        profile = WallDistanceProfile(boundary_vertices=boundary, exponent=2.0)
        points = np.array([[0, 0, 0], [0.004, 0, 0]])
        sf = profile.compute_shape_factors(points)
        # Centre should have higher shape factor than near-wall point
        assert sf[0] > sf[1]

    def test_graceful_degradation_with_mismatch(self):
        """Boundary far from points still produces non-zero output."""
        boundary = np.array([[100, 200, 300], [100.01, 200, 300],
                             [100, 200.01, 300], [100, 200, 300.01]])
        profile = WallDistanceProfile(boundary_vertices=boundary, exponent=2.0)
        sf = profile.compute_shape_factors(SAMPLE_POINTS)
        assert np.all(sf > 0), f"Expected all > 0, got {sf}"


class TestBluntedProfile:
    def test_core_region_is_one(self):
        """Points far from wall should have shape factor = 1."""
        boundary = np.array([
            [0.01, 0, 0], [-0.01, 0, 0],
            [0, 0.01, 0], [0, -0.01, 0],
        ])
        profile = BluntedProfile(boundary_vertices=boundary, core_fraction=0.3)
        sf = profile.compute_shape_factors(np.array([[0, 0, 0]]))
        assert sf[0] == pytest.approx(1.0)


# ----------------------------------------------------------------
# Unit tests: compute_velocities fallback
# ----------------------------------------------------------------

class TestComputeVelocitiesFallback:
    def test_zero_shape_factors_fallback(self):
        """When all shape factors are zero, fallback to uniform non-zero velocities."""
        profile = ZeroProfile()
        vels = profile.compute_velocities(SAMPLE_POINTS, 1e-4, SAMPLE_AREAS, SAMPLE_NORMAL)
        mags = np.linalg.norm(vels, axis=1)
        assert np.all(mags > 0), f"Zero fallback failed: {mags}"

    def test_nan_shape_factors_fallback(self):
        """When shape factors are NaN, fallback to uniform non-zero velocities."""
        profile = NaNProfile()
        vels = profile.compute_velocities(SAMPLE_POINTS, 1e-4, SAMPLE_AREAS, SAMPLE_NORMAL)
        mags = np.linalg.norm(vels, axis=1)
        assert np.all(mags > 0), f"NaN fallback failed: {mags}"
        assert not np.any(np.isnan(vels)), "Velocities should not contain NaN"

    def test_normal_case_unchanged(self):
        """Normal (non-zero) shape factors should not trigger fallback."""
        profile = PlugProfile()
        vels = profile.compute_velocities(SAMPLE_POINTS, 1e-4, SAMPLE_AREAS, SAMPLE_NORMAL)
        mags = np.linalg.norm(vels, axis=1)
        assert np.all(mags > 0)
        # All magnitudes should be equal for plug
        np.testing.assert_allclose(mags, mags[0])


# ----------------------------------------------------------------
# Integration tests: InletMapper with --points-file mismatch
# ----------------------------------------------------------------

@pytest.fixture
def example_stl():
    """Path to the example inlet STL."""
    path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'inlet.stl')
    if not os.path.exists(path):
        pytest.skip("Example inlet.stl not found")
    return path


@pytest.fixture
def example_csv():
    """Path to the example flowrate CSV."""
    path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'flowrate.csv')
    if not os.path.exists(path):
        pytest.skip("Example flowrate.csv not found")
    return path


@pytest.fixture
def example_points():
    """Path to the example points file."""
    path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'points')
    if not os.path.exists(path):
        pytest.skip("Example points file not found")
    return path


@pytest.fixture
def offset_points(example_points):
    """Create a points file with coordinates offset by 500m (simulating mm/m mismatch)."""
    from inlet_mapper.io import PointsFileReader

    reader = PointsFileReader(example_points)
    pts = reader.get_points() + np.array([500, 300, 100])

    tmpdir = tempfile.mkdtemp()
    path = os.path.join(tmpdir, 'points_offset')
    with open(path, 'w') as f:
        f.write(f"{len(pts)}\n(\n")
        for p in pts:
            f.write(f"({p[0]:.10e} {p[1]:.10e} {p[2]:.10e})\n")
        f.write(")\n")

    yield path
    shutil.rmtree(tmpdir)


@pytest.fixture
def output_dir():
    """Temporary output directory."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)


class TestCoordinateMismatchIntegration:
    """Test that all profiles produce non-zero output despite coordinate mismatch."""

    @pytest.mark.parametrize("profile", ["plug", "parabolic", "wall_distance", "womersley", "blunted"])
    def test_mismatch_produces_nonzero(self, profile, example_stl, example_csv, offset_points, output_dir):
        from inlet_mapper import InletMapper

        out = os.path.join(output_dir, f"out_{profile}", "inlet")
        mapper = InletMapper(
            stl_path=example_stl,
            csv_path=example_csv,
            profile=profile,
            points_file=offset_points,
        )
        stats = mapper.generate(output_dir=out, verbose=False)
        assert stats["max_velocity"] > 0, (
            f"{profile} produced zero velocity with coordinate mismatch"
        )

    @pytest.mark.parametrize("profile", ["plug", "parabolic", "wall_distance", "womersley", "blunted"])
    def test_matching_coords_produce_nonzero(self, profile, example_stl, example_csv, example_points, output_dir):
        from inlet_mapper import InletMapper

        out = os.path.join(output_dir, f"out_{profile}", "inlet")
        mapper = InletMapper(
            stl_path=example_stl,
            csv_path=example_csv,
            profile=profile,
            points_file=example_points,
        )
        stats = mapper.generate(output_dir=out, verbose=False)
        assert stats["max_velocity"] > 0, (
            f"{profile} produced zero velocity with matching coordinates"
        )

    @pytest.mark.parametrize("profile", ["parabolic", "womersley"])
    def test_mismatch_close_to_baseline(self, profile, example_stl, example_csv, example_points, offset_points, output_dir):
        """Centroid-based profiles should produce similar results after correction."""
        from inlet_mapper import InletMapper

        out_base = os.path.join(output_dir, f"base_{profile}", "inlet")
        mapper_base = InletMapper(
            stl_path=example_stl, csv_path=example_csv,
            profile=profile, points_file=example_points,
        )
        stats_base = mapper_base.generate(output_dir=out_base, verbose=False)

        out_mis = os.path.join(output_dir, f"mis_{profile}", "inlet")
        mapper_mis = InletMapper(
            stl_path=example_stl, csv_path=example_csv,
            profile=profile, points_file=offset_points,
        )
        stats_mis = mapper_mis.generate(output_dir=out_mis, verbose=False)

        ratio = stats_mis["max_velocity"] / stats_base["max_velocity"]
        assert 0.9 < ratio < 1.1, (
            f"{profile} mismatch ratio {ratio:.3f} too far from 1.0 "
            f"(baseline={stats_base['max_velocity']:.4f}, "
            f"mismatch={stats_mis['max_velocity']:.4f})"
        )
