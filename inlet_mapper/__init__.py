"""
Inlet Mapping Toolkit
=====================

Transform 1D flow rate data into 3D velocity profiles for OpenFOAM.

Supported profile types:
- plug: Uniform velocity
- parabolic: Classic Poiseuille profile
- wall_distance: Distance-to-wall based (for irregular geometries)
- womersley: Pulsatile flow with frequency-dependent effects
- blunted: Flat core with parabolic near-wall

Example:
    >>> from inlet_mapper import InletMapper
    >>> mapper = InletMapper("inlet.stl", "flowrate.csv", profile="wall_distance")
    >>> mapper.generate("boundaryData/inlet")
"""

from .core import InletMapper
from .profiles import (
    PlugProfile,
    ParabolicProfile,
    WallDistanceProfile,
    WomersleyProfile,
    BluntedProfile,
)
from .geometry import InletGeometry
from .io import CSVReader, BoundaryDataWriter

__version__ = "1.0.0"
__author__ = "Jie Wang"

__all__ = [
    "InletMapper",
    "PlugProfile",
    "ParabolicProfile",
    "WallDistanceProfile",
    "WomersleyProfile",
    "BluntedProfile",
    "InletGeometry",
    "CSVReader",
    "BoundaryDataWriter",
]
