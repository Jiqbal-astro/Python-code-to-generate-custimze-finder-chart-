"""Finder charts from sky surveys or local FITS images.

All sky positions, labels, compass arrows and scale bars come from the image's
World Coordinate System (WCS), so they stay correct for any field, declination,
rotation or mirroring.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.patches import Circle, Polygon

SURVEYS = {
    "dss2r": "DSS2 Red",
    "dss2b": "DSS2 Blue",
    "dss2ir": "DSS2 IR",
    "dss": "DSS",
    "2massj": "2MASS-J",
    "2massk": "2MASS-K",
    "sdssr": "SDSSr",
}
NICE_SCALES_ARCSEC = [5, 10, 15, 30, 60, 120, 300, 600, 900, 1200, 1800, 3600, 7200]


class FinderChartError(RuntimeError):
    """An error with a message meant for the user."""


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
def resolve_target(target: str) -> SkyCoord:
    """Parse coordinates ("02h20m30s +57d07m00s", "02:20:30 +57:07:00",
    "35.125 57.1167") or resolve an object name through Sesame (SIMBAD/NED/VizieR)."""
    text = target.strip()
    nums = re.fullmatch(r"([+-]?\d+(?:\.\d+)?)[\s,]+([+-]?\d+(?:\.\d+)?)", text)
    try:
        if nums:
            return SkyCoord(float(nums.group(1)) * u.deg, float(nums.group(2)) * u.deg, frame="icrs")
        if re.search(r"\d+[h:\s]\d+.*[+-]?\d+[d:°\s]\d+", text) and not re.search(r"[a-gi-z]{3,}", text, re.I):
            return SkyCoord(text, unit=(u.hourangle, u.deg), frame="icrs")
    except ValueError as exc:
        raise FinderChartError(f"Could not read the coordinates '{target}': {exc}") from exc
    try:
        return SkyCoord.from_name(text)
    except Exception as exc:
        raise FinderChartError(
            f"Could not resolve '{target}'. Check the name, or give coordinates such as "
            f"'02h20m30s +57d07m00s' or '35.125 57.117'.") from exc


def format_coord(c: SkyCoord) -> str:
    return (f"{c.ra.to_string(unit=u.hour, sep=':', precision=1, pad=True)} "
            f"{c.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)}")


# --------------------------------------------------------------------------
# Images
# --------------------------------------------------------------------------
def fetch_survey_image(center: SkyCoord, fov_arcmin: float, survey: str = "dss2r", pixels: int = 800) -> fits.HDUList:
    """Download an image from NASA SkyView."""
    name = SURVEYS.get(survey.lower(), survey)
    try:
        from astroquery.skyview import SkyView
        images = SkyView.get_images(position=center, survey=[name], width=fov_arcmin * u.arcmin,
                                    height=fov_arcmin * u.arcmin, pixels=pixels)
    except Exception as exc:
        raise FinderChartError(f"Could not download a {name} image from SkyView: {exc}") from exc
    if not images:
        raise FinderChartError(f"SkyView returned no {name} image for this position.")
    return images[0]


def load_fits(path: str | Path) -> tuple[np.ndarray, WCS, fits.Header]:
    """Read the first 2-D image with a celestial WCS from a FITS file."""
    if not Path(path).is_file():
        raise FinderChartError(f"FITS file not found: {path}")
    try:
        hdul = fits.open(path)
    except OSError as exc:
        raise FinderChartError(f"Could not read {Path(path).name} as a FITS file: {exc}") from exc
    with hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:
                    data = data[0]
                wcs = WCS(hdu.header).celestial
                if not wcs.has_celestial:
                    raise FinderChartError(f"{Path(path).name} has no sky coordinates (WCS) in its header.")
                return data, wcs, hdu.header.copy()
    raise FinderChartError(f"{Path(path).name} contains no 2-D image.")


def cutout(data: np.ndarray, wcs: WCS, center: SkyCoord, fov_arcmin: float):
    """Cut a square field around `center` from a larger image."""
    try:
        c = Cutout2D(data, center, (fov_arcmin * u.arcmin, fov_arcmin * u.arcmin), wcs=wcs, mode="partial", fill_value=np.nan)
    except Exception as exc:
        raise FinderChartError(f"The target is not inside the image: {exc}") from exc
    return c.data, c.wcs


def stretch(data: np.ndarray, kind: str = "asinh", low: float = 0.5, high: float = 99.7) -> np.ndarray:
    """Scale to 0–1 between two percentiles, then apply a contrast stretch."""
    finite = np.isfinite(data)
    if not finite.any():
        raise FinderChartError("The image contains no valid pixels.")
    lo, hi = np.percentile(data[finite], [low, high])
    if hi <= lo:
        hi = lo + 1
    x = np.clip((np.nan_to_num(data, nan=lo) - lo) / (hi - lo), 0, 1)
    if kind == "asinh":
        return np.arcsinh(10 * x) / np.arcsinh(10)
    if kind == "sqrt":
        return np.sqrt(x)
    if kind == "log":
        return np.log10(1 + 999 * x) / 3
    if kind == "linear":
        return x
    raise FinderChartError(f"Unknown stretch '{kind}'. Use asinh, sqrt, log or linear.")


# --------------------------------------------------------------------------
# Geometry from the WCS
# --------------------------------------------------------------------------
def pixel_scale_arcsec(wcs: WCS) -> float:
    return float(np.mean(proj_plane_pixel_scales(wcs)) * 3600)


def direction_vectors(wcs: WCS, center: SkyCoord) -> dict:
    """Unit vectors (in pixel coordinates) pointing north and east from `center`."""
    x0, y0 = wcs.world_to_pixel(center)
    step = 5 * pixel_scale_arcsec(wcs) * u.arcsec
    out = {}
    for name, pa in (("N", 0), ("E", 90)):
        x1, y1 = wcs.world_to_pixel(center.directional_offset_by(pa * u.deg, step))
        v = np.array([x1 - x0, y1 - y0], dtype=float)
        out[name] = v / np.hypot(*v)
    return out


def nice_scale_length(fov_arcsec: float) -> float:
    """A round length (arcsec) about a fifth of the field."""
    target = fov_arcsec / 5
    return min(NICE_SCALES_ARCSEC, key=lambda s: abs(np.log(s / target)))


def scale_label(arcsec: float) -> str:
    if arcsec >= 3600:
        return f"{round(arcsec / 3600, 2):g}°"
    if arcsec >= 60:
        return f"{round(arcsec / 60, 1):g}′"
    return f"{round(arcsec):g}″"


def rectangle_corners(center: SkyCoord, width_arcmin: float, height_arcmin: float, pa_deg: float = 0.0) -> SkyCoord:
    """Sky corners of a detector field of `width` x `height` rotated by position angle
    `pa_deg` (east of north) around `center`."""
    w, h = width_arcmin / 2, height_arcmin / 2
    pa = np.radians(pa_deg)
    corners = []
    for dx, dy in ((-w, -h), (w, -h), (w, h), (-w, h)):   # dx towards east, dy towards north
        east = dx * np.cos(pa) + dy * np.sin(pa)
        north = -dx * np.sin(pa) + dy * np.cos(pa)
        sep = np.hypot(east, north) * u.arcmin
        ang = np.degrees(np.arctan2(east, north)) * u.deg
        corners.append(center.directional_offset_by(ang, sep))
    return SkyCoord(corners)


# --------------------------------------------------------------------------
# Comparison stars (Gaia DR3)
# --------------------------------------------------------------------------
@dataclass
class Star:
    source_id: str
    ra: float
    dec: float
    g: float
    bp_rp: Optional[float] = None
    label: str = ""
    sep_arcmin: float = 0.0
    pa_deg: float = 0.0


def query_gaia(center: SkyCoord, radius_arcmin: float, g_limit: float, max_stars: int = 200,
               epoch: Optional[float] = None) -> list[Star]:
    """Gaia DR3 stars brighter than `g_limit` within `radius_arcmin`, optionally moved
    to `epoch` (decimal year) with their proper motions."""
    try:
        from astroquery.gaia import Gaia
        q = (f"SELECT TOP {int(max_stars)} source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp "
             f"FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
             f"CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {radius_arcmin / 60})) "
             f"AND phot_g_mean_mag < {g_limit} ORDER BY phot_g_mean_mag")
        table = Gaia.launch_job(q).get_results()
    except Exception as exc:
        raise FinderChartError(f"Could not query Gaia DR3: {exc}") from exc
    stars = []
    for row in table:
        ra, dec = float(row["ra"]), float(row["dec"])
        if epoch is not None and np.isfinite(row["pmra"]) and np.isfinite(row["pmdec"]):
            dt = epoch - 2016.0                                  # Gaia DR3 reference epoch
            ra += float(row["pmra"]) * dt / 3.6e6 / np.cos(np.radians(dec))
            dec += float(row["pmdec"]) * dt / 3.6e6
        bp = row["bp_rp"]
        stars.append(Star(str(row["source_id"]), ra, dec, float(row["phot_g_mean_mag"]),
                          None if np.ma.is_masked(bp) or not np.isfinite(bp) else float(bp)))
    return stars


def label_stars(stars: Sequence[Star], center: SkyCoord, n_labels: int, exclude_arcsec: float = 3.0) -> list[Star]:
    """Number the brightest stars (excluding the target itself) and record their
    separation and position angle from the target. AAVSO-style labels use the
    magnitude without its decimal point (G = 11.43 -> "114")."""
    out = []
    for s in stars:
        c = SkyCoord(s.ra * u.deg, s.dec * u.deg)
        s.sep_arcmin = float(center.separation(c).arcmin)
        s.pa_deg = float(center.position_angle(c).deg)
        if s.sep_arcmin * 60 > exclude_arcsec:
            out.append(s)
    out.sort(key=lambda s: s.g)
    for s in out[:n_labels]:
        s.label = f"{round(s.g * 10):d}"
    return out


def write_star_table(stars: Sequence[Star], path: str | Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["label", "gaia_dr3_source_id", "ra_deg", "dec_deg", "G_mag", "BP_RP", "sep_arcmin", "pa_deg"])
        for s in stars:
            if s.label:
                w.writerow([s.label, s.source_id, f"{s.ra:.6f}", f"{s.dec:.6f}", f"{s.g:.3f}",
                            "" if s.bp_rp is None else f"{s.bp_rp:.3f}", f"{s.sep_arcmin:.3f}", f"{s.pa_deg:.1f}"])


# --------------------------------------------------------------------------
# The chart
# --------------------------------------------------------------------------
@dataclass
class ChartOptions:
    title: str = ""
    invert: bool = True                 # dark stars on a white background (prints well)
    stretch: str = "asinh"
    grid: bool = True
    reticle: bool = True
    mirror: bool = False                # east–west flip, as seen through a star diagonal
    rotate180: bool = False             # as seen through a Newtonian or inverting eyepiece
    instrument: Optional[tuple[float, float]] = None   # detector field, arcmin (width, height)
    instrument_pa: float = 0.0
    circles_arcmin: Sequence[float] = field(default_factory=list)
    stars: Sequence[Star] = field(default_factory=list)
    survey_label: str = ""
    size_inches: float = 8.0
    dpi: int = 150


def make_chart(data: np.ndarray, wcs: WCS, center: SkyCoord, opts: ChartOptions):
    """Draw the finder chart. Returns the matplotlib figure."""
    img = stretch(data, opts.stretch)
    if opts.invert:
        img = 1 - img
    img = np.ma.masked_where(~np.isfinite(data), img)          # outside the survey image
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad("#d9dde3" if opts.invert else "#2a2f38")
    ink = "black" if opts.invert else "white"
    accent = "#c81e1e" if opts.invert else "#ff6b6b"
    fig = plt.figure(figsize=(opts.size_inches, opts.size_inches * 1.08))
    ax = fig.add_axes([0.13, 0.11, 0.82, 0.78], projection=wcs)
    ax.imshow(img, origin="lower", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    ny, nx = data.shape
    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(-0.5, ny - 0.5)
    # mirror = flip east–west; rotate180 = flip both axes; both together = flip north–south only
    if opts.mirror != opts.rotate180:
        ax.invert_xaxis()
    if opts.rotate180:
        ax.invert_yaxis()

    ra, dec = ax.coords[0], ax.coords[1]
    ra.set_axislabel("Right ascension (ICRS)")
    dec.set_axislabel("Declination (ICRS)")
    ra.set_major_formatter("hh:mm:ss")
    dec.set_major_formatter("dd:mm")
    ra.set_ticks(number=5)
    dec.set_ticks(number=5)
    if opts.grid:
        ax.coords.grid(color=ink, alpha=0.18, linestyle=":", linewidth=0.8)

    scale = pixel_scale_arcsec(wcs)
    fov_x, fov_y = nx * scale, ny * scale
    x0, y0 = wcs.world_to_pixel(center)

    # target reticle: four ticks with a gap, so the target itself stays visible
    if opts.reticle:
        gap, length = 0.025 * nx, 0.06 * nx
        for vx, vy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ax.plot([x0 + vx * gap, x0 + vx * (gap + length)], [y0 + vy * gap, y0 + vy * (gap + length)],
                    color=accent, lw=2, solid_capstyle="butt")

    # instrument field and circles
    if opts.instrument:
        corners = rectangle_corners(center, opts.instrument[0], opts.instrument[1], opts.instrument_pa)
        px, py = wcs.world_to_pixel(corners)
        ax.add_patch(Polygon(np.column_stack([px, py]), closed=True, fill=False, ec=accent, lw=1.4, ls="--"))
        ax.text(px[3], py[3], f" {opts.instrument[0]:g}′ × {opts.instrument[1]:g}′", color=accent, fontsize=9,
                va="bottom", ha="left")
    for r in opts.circles_arcmin:
        ax.add_patch(Circle((x0, y0), r * 60 / scale, fill=False, ec=accent, lw=1, ls=":"))

    # comparison stars
    for s in opts.stars:
        if not s.label:
            continue
        sx, sy = wcs.world_to_pixel(SkyCoord(s.ra * u.deg, s.dec * u.deg))
        if not (0 <= sx < nx and 0 <= sy < ny):
            continue
        ax.add_patch(Circle((sx, sy), 0.014 * nx, fill=False, ec="#1f5fbf", lw=1.1))
        ax.text(sx + 0.018 * nx, sy + 0.012 * ny, s.label, color="#1f5fbf", fontsize=8.5, fontweight="bold")

    # compass from the WCS, drawn in the top-left corner of the image (in screen space)
    vec = direction_vectors(wcs, center)
    flip_x = -1 if ax.xaxis_inverted() else 1
    flip_y = -1 if ax.yaxis_inverted() else 1
    ox, oy, L = 0.12, 0.84, 0.09
    for name, v in vec.items():
        dx, dy = v[0] * flip_x * L, v[1] * flip_y * L
        ax.annotate("", xy=(ox + dx, oy + dy), xytext=(ox, oy), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=ink, lw=1.6))
        ax.text(ox + dx * 1.3, oy + dy * 1.3, name, transform=ax.transAxes, color=ink,
                ha="center", va="center", fontsize=12, fontweight="bold")

    # scale bar, bottom right
    length = nice_scale_length(fov_x)
    frac = (length / scale) / nx
    bx1, by = 0.95, 0.06
    ax.plot([bx1 - frac, bx1], [by, by], transform=ax.transAxes, color=ink, lw=2.5, solid_capstyle="butt")
    ax.text(bx1 - frac / 2, by + 0.02, scale_label(length), transform=ax.transAxes, color=ink, ha="center", fontsize=10)

    # titles
    title = opts.title or format_coord(center)
    fig.text(0.13, 0.955, title, fontsize=15, fontweight="bold", ha="left", va="center")
    fov_text = f"{scale_label(round(fov_x))} × {scale_label(round(fov_y))}"
    orient = {(False, False): "standard view", (True, False): "mirrored east–west (star diagonal)",
              (False, True): "rotated 180° (Newtonian / inverting eyepiece)", (True, True): "mirrored north–south"}[
        (bool(opts.mirror), bool(opts.rotate180))]
    sub = f"{format_coord(center)} (ICRS)   |   field {fov_text}   |   {orient}"
    if opts.survey_label:
        sub += f"   |   {opts.survey_label}"
    fig.text(0.13, 0.922, sub, fontsize=9, color="#444", ha="left", va="center")
    missing = float(np.mean(~np.isfinite(data)))
    if missing > 0.01:
        fig.text(0.13, 0.055, f"Grey area ({missing:.0%} of the field) lies outside the image.", fontsize=8.5, color="#a33")
    if opts.stars and any(s.label for s in opts.stars):
        fig.text(0.13, 0.03, "Blue labels: Gaia DR3 G magnitude × 10 (114 = 11.4 mag)", fontsize=8.5, color="#1f5fbf")
    fig.text(0.95, 0.03, f"Made {datetime.now(timezone.utc):%Y-%m-%d} with finderchart", fontsize=7.5, color="#888", ha="right")
    return fig


def save_chart(fig, path: str | Path, dpi: int = 150) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return path
