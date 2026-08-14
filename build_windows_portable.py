from __future__ import annotations

import os
import shutil
import subprocess
import sys
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUILD_ENV = ROOT / ".build-venv"
PYTHON = BUILD_ENV / "Scripts" / "python.exe"


def run(command: list[str]) -> None:
    print("\n>", " ".join(map(str, command)))
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    if os.name != "nt":
        print("Build the Windows executable on Windows. PyInstaller bundles separately for each operating system.")
        return 2
    if not (3, 11) <= sys.version_info[:2] <= (3, 14):
        print("Python 3.11–3.14 is required to build the portable application.")
        return 2
    if not PYTHON.exists():
        venv.EnvBuilder(with_pip=True).create(BUILD_ENV)
    run([str(PYTHON), "-m", "pip", "install", "--upgrade", "pip"])
    run([str(PYTHON), "-m", "pip", "install", "--requirement", "requirements.txt", "--requirement", "build-requirements.txt"])
    dist = ROOT / "dist"
    work = ROOT / "build"
    if dist.exists():
        shutil.rmtree(dist)
    if work.exists():
        shutil.rmtree(work)
    run([
        str(PYTHON), "-m", "PyInstaller",
        "--noconfirm", "--clean", "--windowed", "--onedir",
        "--name", "LCMS_Compound_Curation",
        "--collect-all", "PySide6",
        "--hidden-import", "openpyxl.cell._writer",
        "app.py",
    ])
    folder = dist / "LCMS_Compound_Curation"
    archive_base = dist / "LCMS_Compound_Curation_Windows_Portable"
    archive = shutil.make_archive(str(archive_base), "zip", root_dir=folder)
    print(f"\nPortable Windows build created:\n{folder}\n{archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

