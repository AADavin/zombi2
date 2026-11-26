"""Export the approved four-level overview from its editable SVG source.

Run: python figures/scripts/fig_overview.py
Requires librsvg's rsvg-convert on PATH. This redraws the packaged illustration;
it does not rerun the historical simulations. Provenance: assets/overview.json.
"""
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main():
    converter = shutil.which("rsvg-convert")
    if converter is None:
        raise SystemExit("Install librsvg and put rsvg-convert on PATH.")
    source = ROOT / "assets/overview.svg"
    target = ROOT / "assets/overview.png"
    temporary = target.with_suffix(".tmp.png")
    subprocess.run([converter, "-w", "3200", "-o", str(temporary), str(source)], check=True)
    temporary.replace(target)
    print(target)


if __name__ == "__main__":
    main()
