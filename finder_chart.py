"""Backward-compatible entry point: the original script made a chart of the
Double Cluster (NGC 869 and NGC 884). Run it with no arguments to make that chart,
or pass any finderchart arguments, e.g.  python finder_chart.py "M 57" --fov 10"""
import sys

from finderchart.cli import main

if __name__ == "__main__":
    args = sys.argv[1:] or ["02h20m30s +57d07m00s", "--fov", "45", "--title", "NGC 869 and NGC 884 (the Double Cluster)",
                            "-o", "ngc_884_869_image.png"]
    sys.exit(main(args))
