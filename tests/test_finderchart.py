"""Offline tests. Network services (SkyView, Gaia, Sesame) are replaced by fakes;
the image tests use a real Digitized Sky Survey cutout with its WCS."""
import csv
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from finderchart import cli, core

DSS = Path(__file__).with_name("horsehead_dss.fits")
CENTER = SkyCoord("05h40m59s -02d27m30s")


def synthetic_wcs(rotation_deg=0.0, dec=57.1, scale_arcsec=2.0, n=300):
    """A standard (north up, east left) TAN WCS, optionally rotated."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [35.125, dec]
    w.wcs.crpix = [n / 2, n / 2]
    s, t = scale_arcsec / 3600, np.radians(rotation_deg)
    w.wcs.cd = [[-s * np.cos(t), s * np.sin(t)], [s * np.sin(t), s * np.cos(t)]]
    return w


# ------------------------------------------------------------------ targets
@pytest.mark.parametrize("text", ["02h20m30s +57d07m00s", "02:20:30 +57:07:00", "35.125 57.11667", "35.125, 57.11667"])
def test_resolve_coordinates(text):
    c = core.resolve_target(text)
    assert c.ra.deg == pytest.approx(35.125, abs=1e-4)
    assert c.dec.deg == pytest.approx(57.11667, abs=1e-4)


def test_resolve_name_uses_sesame(monkeypatch):
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(lambda cls, name: SkyCoord(283.396 * u.deg, 33.029 * u.deg)))
    assert core.resolve_target("M 57").ra.deg == pytest.approx(283.396)


def test_unknown_name_gives_helpful_error(monkeypatch):
    def fail(cls, name):
        raise ValueError("not found")
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(fail))
    with pytest.raises(core.FinderChartError, match="Could not resolve"):
        core.resolve_target("Not A Real Object")


# --------------------------------------------------------------- geometry
def test_real_dss_image_is_north_up_east_left():
    data, wcs, _ = core.load_fits(DSS)
    v = core.direction_vectors(wcs, CENTER)
    assert v["N"][1] > 0.99            # north points up
    assert v["E"][0] < -0.99           # east points LEFT (the original script drew it to the right)


def test_ra_span_includes_cos_dec_factor():
    # The original script labelled RA as a linear offset in degrees with no cos(Dec):
    # at Dec +57 that is wrong by a factor of 1.84. The WCS gets it right.
    w = synthetic_wcs(dec=57.1167, scale_arcsec=3.375, n=800)     # 45 arcmin field
    left = w.pixel_to_world(0, 400)
    right = w.pixel_to_world(799, 400)
    dra = ((left.ra - right.ra).wrap_at(180 * u.deg)).deg
    assert dra > 0                                              # RA increases to the left
    assert dra == pytest.approx(0.75 / np.cos(np.radians(57.1167)), rel=0.01)


@pytest.mark.parametrize("rot", [0, 30, 135])
def test_compass_follows_rotated_wcs(rot):
    w = synthetic_wcs(rotation_deg=rot)
    c = w.pixel_to_world(150, 150)
    v = core.direction_vectors(w, c)
    assert np.degrees(np.arccos(np.clip(np.dot(v["N"], v["E"]), -1, 1))) == pytest.approx(90, abs=0.5)
    # east is 90 degrees counter-clockwise from north on a standard (non-mirrored) image
    assert v["N"][0] * v["E"][1] - v["N"][1] * v["E"][0] > 0


def test_scale_bar_and_labels():
    assert core.nice_scale_length(45 * 60) == 600            # 10' for a 45' field
    assert core.nice_scale_length(12 * 60) == 120
    assert core.scale_label(600) == "10′" and core.scale_label(30) == "30″" and core.scale_label(3600) == "1°"


@pytest.mark.parametrize("pa", [0, 30, 90])
def test_instrument_rectangle(pa):
    corners = core.rectangle_corners(CENTER, 20, 10, pa)
    seps = CENTER.separation(corners).arcmin
    assert np.allclose(seps, np.hypot(10, 5), rtol=1e-3)
    # corners 1 and 2 bound the east side of the field: its midpoint lies at PA + 90
    mid = SkyCoord(ra=np.mean(corners[1:3].ra.deg) * u.deg, dec=np.mean(corners[1:3].dec.deg) * u.deg)
    assert ((CENTER.position_angle(mid).deg - (pa + 90)) + 180) % 360 - 180 == pytest.approx(0, abs=0.5)
    assert CENTER.separation(mid).arcmin == pytest.approx(10, rel=1e-3)


def test_stretch_range_and_errors():
    data = np.random.default_rng(1).normal(100, 10, (50, 50))
    for kind in ("asinh", "sqrt", "log", "linear"):
        s = core.stretch(data, kind)
        assert s.min() >= 0 and s.max() <= 1
    with pytest.raises(core.FinderChartError):
        core.stretch(data, "bogus")
    with pytest.raises(core.FinderChartError):
        core.stretch(np.full((5, 5), np.nan))


def test_cutout_outside_image_is_reported():
    data, wcs, _ = core.load_fits(DSS)
    with pytest.raises(core.FinderChartError, match="not inside"):
        core.cutout(data, wcs, SkyCoord("12h00m00s +30d00m00s"), 5)


def test_load_fits_without_wcs(tmp_path):
    p = tmp_path / "plain.fits"
    fits.PrimaryHDU(np.zeros((10, 10))).writeto(p)
    with pytest.raises(core.FinderChartError, match="no sky coordinates"):
        core.load_fits(p)


# ------------------------------------------------------------------- Gaia
def fake_gaia_table():
    return Table(rows=[
        (1, CENTER.ra.deg, CENTER.dec.deg, 0.0, 0.0, 9.0, 0.8),                     # the target itself
        (2, CENTER.ra.deg + 0.01, CENTER.dec.deg, 100.0, -50.0, 11.43, 0.6),
        (3, CENTER.ra.deg, CENTER.dec.deg + 0.02, 0.0, 0.0, 12.07, np.nan),
        (4, CENTER.ra.deg - 0.02, CENTER.dec.deg - 0.01, 0.0, 0.0, 13.5, 1.1),
    ], names=("source_id", "ra", "dec", "pmra", "pmdec", "phot_g_mean_mag", "bp_rp"))


class FakeJob:
    def get_results(self):
        return fake_gaia_table()


def test_gaia_labels_proper_motion_and_csv(monkeypatch, tmp_path):
    from astroquery.gaia import Gaia
    monkeypatch.setattr(Gaia, "launch_job", lambda q: FakeJob())
    stars = core.query_gaia(CENTER, 10, 14, epoch=2026.0)
    moved = stars[1]
    assert moved.dec == pytest.approx(CENTER.dec.deg - 50 * 10 / 3.6e6)          # 10 yr of pmdec
    labelled = core.label_stars(stars, CENTER, 2)
    assert [s.label for s in labelled[:2]] == ["114", "121"]                      # target excluded
    assert labelled[2].label == ""
    out = tmp_path / "stars.csv"
    core.write_star_table(labelled, out)
    rows = list(csv.DictReader(open(out)))
    assert [r["label"] for r in rows] == ["114", "121"] and rows[1]["BP_RP"] == ""
    assert float(rows[1]["pa_deg"]) == pytest.approx(0, abs=1)                     # star 3 is due north


# --------------------------------------------------------------- the chart
def test_chart_from_real_dss_image(tmp_path):
    data, wcs, _ = core.load_fits(DSS)
    data, wcs = core.cutout(data, wcs, CENTER, 7)
    opts = core.ChartOptions(title="Horsehead", instrument=(5, 3), instrument_pa=20, circles_arcmin=[1, 2])
    path = core.save_chart(core.make_chart(data, wcs, CENTER, opts), tmp_path / "c.png")
    assert path.stat().st_size > 50_000


@pytest.mark.parametrize("mirror,rot", [(True, False), (False, True), (True, True)])
def test_mirrored_and_rotated_compass(tmp_path, mirror, rot):
    data, wcs, _ = core.load_fits(DSS)
    fig = core.make_chart(data, wcs, CENTER, core.ChartOptions(mirror=mirror, rotate180=rot))
    ax = fig.axes[0]
    # find the E label: on a mirrored chart east must appear on the RIGHT of the compass origin
    e = [t for t in ax.texts if t.get_text() == "E"][0]
    x_e = e.get_position()[0]
    assert (x_e > 0.12) == (mirror != rot)
    core.save_chart(fig, tmp_path / "m.png")


# -------------------------------------------------------------------- CLI
def test_cli_with_local_fits(tmp_path):
    out = tmp_path / "hh.pdf"
    rc = cli.main(["05h40m59s -02d27m30s", "--fits", str(DSS), "--fov", "6", "-o", str(out),
                   "--instrument", "4x3", "--pa", "15", "--circle", "1", "--mirror", "--stretch", "sqrt", "--no-invert"])
    assert rc == 0 and out.stat().st_size > 10_000


def test_cli_downloads_from_skyview(monkeypatch, tmp_path):
    from astroquery.skyview import SkyView
    calls = {}

    def fake_get_images(position, survey, width, height, pixels):
        calls.update(survey=survey, width=width)
        data, wcs, hdr = core.load_fits(DSS)
        return [fits.HDUList([fits.PrimaryHDU(data, wcs.to_header())])]
    monkeypatch.setattr(SkyView, "get_images", fake_get_images)
    out = tmp_path / "s.png"
    assert cli.main(["05h40m59s -02d27m30s", "--fov", "7", "--survey", "dss2b", "-o", str(out)]) == 0
    assert calls["survey"] == ["DSS2 Blue"] and calls["width"].to(u.arcmin).value == 7
    assert out.exists()


def test_cli_gaia_option(monkeypatch, tmp_path):
    from astroquery.gaia import Gaia
    monkeypatch.setattr(Gaia, "launch_job", lambda q: FakeJob())
    csv_out = tmp_path / "comps.csv"
    rc = cli.main(["05h40m59s -02d27m30s", "--fits", str(DSS), "--fov", "6", "--gaia", "14", "--labels", "3",
                   "--stars-csv", str(csv_out), "-o", str(tmp_path / "g.png")])
    assert rc == 0 and len(list(csv.DictReader(open(csv_out)))) == 3


def test_cli_batch(tmp_path, capsys):
    batch = tmp_path / "targets.csv"
    batch.write_text("target,fov,title\n05h40m59s -02d27m30s,5,Horsehead\n12h00m00s +30d00m00s,5,Off image\n")
    rc = cli.main(["--batch", str(batch), "--fits", str(DSS), "-o", str(tmp_path / "charts")])
    assert rc == 1                                            # one target skipped, reported
    assert (tmp_path / "charts" / "05h40m59s_-02d27m30s.png").exists()
    assert "skipped 12h00m00s +30d00m00s" in capsys.readouterr().err


def test_cli_errors_are_friendly(capsys, tmp_path):
    assert cli.main([]) == 2
    assert cli.main(["05h40m59s -02d27m30s", "--fits", str(tmp_path / "missing.fits")]) == 1
    assert "FITS file not found" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        cli.main(["x", "--instrument", "banana"])


def test_backward_compatible_script_defaults(monkeypatch):
    import runpy
    seen = {}
    monkeypatch.setattr(cli, "main", lambda args: seen.setdefault("args", args) and 0)
    monkeypatch.setattr("sys.argv", ["finder_chart.py"])
    with pytest.raises(SystemExit):
        runpy.run_path(str(Path(__file__).resolve().parents[1] / "finder_chart.py"), run_name="__main__")
    assert "02h20m30s +57d07m00s" in seen["args"] and "45" in seen["args"]
