#!/usr/bin/env python3
"""
Inlet Mapping Toolkit - Command Line Interface
==============================================

Transform 1D flow rate data into 3D velocity profiles for OpenFOAM.

Usage:
    # Wall-distance profile from CSV
    python map_inlet.py inlet.stl flowrate.csv --profile wall_distance

    # Constant flow with parabolic profile
    python map_inlet.py inlet.stl --flowrate 5.0 --profile parabolic

    # Womersley profile for pulsatile flow
    python map_inlet.py inlet.stl flowrate.csv --profile womersley --viscosity 3.5e-6

Author: Jie Wang
Date: December 2025
"""

import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Map 1D flow rate to 3D velocity profiles for OpenFOAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Wall-distance profile from CSV (best for irregular inlets)
  python map_inlet.py inlet.stl flowrate.csv --profile wall_distance

  # Parabolic profile with constant flow
  python map_inlet.py inlet.stl --flowrate 5.0 --profile parabolic

  # Womersley profile for pulsatile flow
  python map_inlet.py inlet.stl flowrate.csv --profile womersley

  # Custom wall-distance exponent
  python map_inlet.py inlet.stl flowrate.csv --profile wall_distance --exponent 1.5

  # Use existing mesh points file
  python map_inlet.py inlet.stl flowrate.csv --points-file boundaryData/inlet/points

Profile Types:
  plug          - Uniform velocity (simplest, unrealistic)
  parabolic     - Classic Poiseuille profile (circular inlets)
  wall_distance - Distance-to-wall based (irregular geometries)
  womersley     - Pulsatile with frequency effects (advanced)
  blunted       - Flat core + parabolic wall (transitional flow)
        """
    )

    # Required arguments
    parser.add_argument(
        "stl",
        help="Path to inlet STL file"
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default=None,
        help="Path to flow rate CSV file (optional if --flowrate provided)"
    )

    # Profile selection
    parser.add_argument(
        "-p", "--profile",
        choices=["plug", "parabolic", "wall_distance", "womersley", "blunted"],
        default="parabolic",
        help="Velocity profile type (default: parabolic)"
    )

    # Flow specification
    parser.add_argument(
        "-f", "--flowrate",
        type=float,
        default=None,
        help="Constant flow rate in L/min (alternative to CSV)"
    )

    # Output
    parser.add_argument(
        "-o", "--output",
        default="boundaryData/inlet",
        help="Output directory (default: boundaryData/inlet)"
    )

    # Profile-specific options
    parser.add_argument(
        "--exponent",
        type=float,
        default=2.0,
        help="Exponent for wall_distance profile (default: 2.0)"
    )
    parser.add_argument(
        "--viscosity",
        type=float,
        default=3.5e-6,
        help="Kinematic viscosity for womersley profile (default: 3.5e-6 m²/s)"
    )
    parser.add_argument(
        "--core-fraction",
        type=float,
        default=0.3,
        help="Core fraction for blunted profile (default: 0.3)"
    )
    parser.add_argument(
        "--n-harmonics",
        type=int,
        default=10,
        help="Number of Fourier harmonics for womersley (default: 10)"
    )

    # Geometry options
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="STL scale factor to convert to meters (default: 1.0)"
    )
    parser.add_argument(
        "--flip-normal",
        action="store_true",
        help="Flip inlet normal direction"
    )
    parser.add_argument(
        "--points-file",
        default=None,
        help="Use existing mesh points file instead of STL face centres"
    )

    # Output options
    parser.add_argument(
        "--generate-points-only",
        action="store_true",
        help="Only generate points file (no velocity data)"
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Don't remove existing timestep directories"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Minimal output"
    )

    # Info commands
    parser.add_argument(
        "--info",
        action="store_true",
        help="Show geometry and data info without generating output"
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Preview velocity statistics without generating full output"
    )

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.stl).exists():
        print(f"Error: STL file not found: {args.stl}")
        sys.exit(1)

    if args.csv is None and args.flowrate is None:
        print("Error: Either CSV file or --flowrate must be provided")
        sys.exit(1)

    if args.csv and not Path(args.csv).exists():
        print(f"Error: CSV file not found: {args.csv}")
        sys.exit(1)

    # Import here to allow --help without dependencies
    try:
        from inlet_mapper import InletMapper
    except ImportError as e:
        print(f"Error importing inlet_mapper: {e}")
        print("\nMake sure you're running from the toolkit directory or have it installed.")
        sys.exit(1)

    # Build profile kwargs
    profile_kwargs = {
        "exponent": args.exponent,
        "viscosity": args.viscosity,
        "core_fraction": args.core_fraction,
        "n_harmonics": args.n_harmonics,
    }

    # Create mapper
    try:
        mapper = InletMapper(
            stl_path=args.stl,
            csv_path=args.csv,
            profile=args.profile,
            flowrate=args.flowrate,
            scale_factor=args.scale,
            flip_normal=args.flip_normal,
            points_file=args.points_file,
            **profile_kwargs
        )
    except Exception as e:
        print(f"Error creating mapper: {e}")
        sys.exit(1)

    # Handle info/preview modes
    if args.info:
        import json
        stats = mapper.get_statistics()
        print("\n" + "="*60)
        print("INLET MAPPING INFO")
        print("="*60)
        print(json.dumps(stats, indent=2, default=str))
        return

    if args.preview:
        preview = mapper.preview()
        print("\n" + "="*60)
        print("VELOCITY PREVIEW")
        print("="*60)
        print(f"Flow rate: {preview['flow_rate_m3s']:.6e} m³/s")
        print(f"Min velocity:  {preview['magnitudes'].min():.4f} m/s")
        print(f"Max velocity:  {preview['magnitudes'].max():.4f} m/s")
        print(f"Mean velocity: {preview['magnitudes'].mean():.4f} m/s")
        return

    if args.generate_points_only:
        from inlet_mapper.io import BoundaryDataWriter
        writer = BoundaryDataWriter(args.output)

        if args.points_file:
            from inlet_mapper.io import PointsFileReader
            reader = PointsFileReader(args.points_file)
            points = reader.get_points()
        else:
            points = mapper.geometry.get_face_centres()

        writer.write_points(points)
        print(f"Wrote points file to {args.output}/points ({len(points)} points)")
        return

    # Generate full output
    try:
        stats = mapper.generate(
            output_dir=args.output,
            clean_existing=not args.no_clean,
            verbose=not args.quiet
        )

        if not args.quiet:
            print(f"\n✅ Successfully generated inlet boundary data")
            print(f"   Output: {args.output}")
            print(f"   Points: {stats['n_points']}")
            print(f"   Timesteps: {stats['n_timesteps']}")

    except Exception as e:
        print(f"Error generating output: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
