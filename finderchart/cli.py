"""Command-line interface:  python -m finderchart "M 57" --fov 15"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

from . import core
from .core import ChartOptions, FinderChartError


def parse_size(text: str) -> tuple[float, float]:
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*[x×*]\s*(\d+(?:\.\d+)?)\s*", text or "")
    if not m:
        raise argparse.ArgumentTypeError("use WIDTHxHEIGHT in arcminutes, for example 21.5x14.3")
    return float(m.group(1)), float(m.group(2))


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9+-]+", "_", text).strip("_") or "target"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="finderchart",
        description="Make a finder chart for a sky target from a survey image (NASA SkyView) or your own FITS file.",
        epilog='Examples:\n  finderchart "NGC 869" --fov 45\n'
               '  finderchart "02h20m30s +57d07m00s" --fov 45 --survey dss2b\n'
               '  finderchart "SS Cyg" --fov 20 --instrument 21.5x14.3 --gaia 13 --labels 12\n'
               '  finderchart --batch targets.csv --fov 30 -o charts/',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("target", nargs="?", help="object name (resolved by SIMBAD) or coordinates, e.g. '02h20m30s +57d07m00s' or '35.125 57.117'")
    p.add_argument("--fov", type=float, default=30.0, help="field of view in arcminutes (default 30)")
    p.add_argument("--survey", default="dss2r", help=f"SkyView survey: {', '.join(core.SURVEYS)} or any SkyView name (default dss2r)")
    p.add_argument("--fits", help="use this FITS image instead of downloading one (needs a WCS in its header)")
    p.add_argument("--pixels", type=int, default=800, help="image size to request from SkyView (default 800)")
    p.add_argument("-o", "--output", help="output file (.png, .pdf or .svg) or, with --batch, a directory")
    p.add_argument("--title", help="chart title (default: the target)")
    p.add_argument("--stretch", default="asinh", choices=["asinh", "sqrt", "log", "linear"], help="contrast stretch (default asinh)")
    p.add_argument("--no-invert", action="store_true", help="white stars on black instead of black on white")
    p.add_argument("--no-grid", action="store_true", help="hide the RA/Dec grid")
    p.add_argument("--no-reticle", action="store_true", help="hide the target marker")
    p.add_argument("--mirror", action="store_true", help="flip east–west, as seen through a star diagonal")
    p.add_argument("--rotate180", action="store_true", help="rotate 180°, as seen through a Newtonian or inverting eyepiece")
    p.add_argument("--instrument", type=parse_size, metavar="WxH", help="draw your camera's field, in arcmin, e.g. 21.5x14.3")
    p.add_argument("--pa", type=float, default=0.0, help="position angle of the camera field, degrees east of north (default 0)")
    p.add_argument("--circle", type=float, action="append", default=[], metavar="ARCMIN", help="draw a circle of this radius (repeatable)")
    p.add_argument("--gaia", type=float, metavar="GMAG", help="mark Gaia DR3 stars brighter than this G magnitude as comparison stars")
    p.add_argument("--labels", type=int, default=15, help="number of Gaia stars to label (default 15)")
    p.add_argument("--epoch", type=float, help="move Gaia stars to this epoch with their proper motions (default: today)")
    p.add_argument("--stars-csv", help="also write the labelled Gaia stars to this CSV file")
    p.add_argument("--batch", help="CSV with a column 'target' (optional: 'fov', 'title') to make many charts")
    p.add_argument("--dpi", type=int, default=150, help="output resolution (default 150)")
    return p


def chart_for(target: str, args, output: Path, fov: float | None = None, title: str | None = None) -> Path:
    fov = fov or args.fov
    center = core.resolve_target(target)
    if args.fits:
        data, wcs, _ = core.load_fits(args.fits)
        data, wcs = core.cutout(data, wcs, center, fov)
        survey_label = f"image: {Path(args.fits).name}"
    else:
        hdu = core.fetch_survey_image(center, fov, args.survey, args.pixels)[0]
        data, wcs = np.asarray(hdu.data, dtype=float), WCS(hdu.header).celestial
        survey_label = core.SURVEYS.get(args.survey.lower(), args.survey)
    stars = []
    if args.gaia is not None:
        epoch = args.epoch if args.epoch is not None else _decimal_year()
        stars = core.query_gaia(center, fov / np.sqrt(2), args.gaia, epoch=epoch)
        stars = core.label_stars(stars, center, args.labels)
        if args.stars_csv:
            core.write_star_table(stars, args.stars_csv)
    opts = ChartOptions(
        title=title or args.title or target, invert=not args.no_invert, stretch=args.stretch, grid=not args.no_grid,
        reticle=not args.no_reticle, mirror=args.mirror, rotate180=args.rotate180, instrument=args.instrument,
        instrument_pa=args.pa, circles_arcmin=args.circle, stars=stars, survey_label=survey_label, dpi=args.dpi)
    fig = core.make_chart(data, wcs, center, opts)
    return core.save_chart(fig, output, args.dpi)


def _decimal_year() -> float:
    now = datetime.now(timezone.utc)
    start = datetime(now.year, 1, 1, tzinfo=timezone.utc)
    return now.year + (now - start).total_seconds() / (365.25 * 86400)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.batch:
            outdir = Path(args.output or "charts")
            with open(args.batch, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            if not rows or "target" not in rows[0]:
                raise FinderChartError("The batch CSV needs a header row with a 'target' column.")
            failures = 0
            for row in rows:
                name = row["target"].strip()
                if not name:
                    continue
                try:
                    fov = float(row["fov"]) if row.get("fov") else None
                    out = chart_for(name, args, outdir / f"{slug(name)}.png", fov, row.get("title") or None)
                    print(f"saved {out}")
                except FinderChartError as exc:
                    failures += 1
                    print(f"skipped {name}: {exc}", file=sys.stderr)
            return 1 if failures else 0
        if not args.target:
            build_parser().print_usage(sys.stderr)
            print("finderchart: give a target (or --batch FILE)", file=sys.stderr)
            return 2
        out = Path(args.output or f"finder_{slug(args.target)}.png")
        print(f"saved {chart_for(args.target, args, out)}")
        return 0
    except FinderChartError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
