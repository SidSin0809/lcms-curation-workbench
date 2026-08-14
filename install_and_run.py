from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import struct
import subprocess
import sys
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv"
REQUIREMENTS = ROOT / "requirements.txt"
MARKER = VENV / ".lcms_dependencies"


def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def dependency_signature() -> str:
    payload = REQUIREMENTS.read_bytes() + f"{sys.version_info.major}.{sys.version_info.minor}".encode()
    return hashlib.sha256(payload).hexdigest()


def run(command: list[str]) -> None:
    print("\n>", " ".join(command))
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Install dependencies into an isolated local environment and launch the LC–MS GUI.")
    parser.add_argument("--repair", action="store_true", help="Recreate the isolated environment before launching.")
    args = parser.parse_args()
    if not (3, 11) <= sys.version_info[:2] <= (3, 14):
        print("Python 3.11–3.14 is required. Install a supported 64-bit Python from python.org and run this launcher again.")
        return 2
    if struct.calcsize("P") * 8 != 64:
        print("A 64-bit Python installation is required for the pinned Qt and data-analysis dependencies.")
        return 2
    print(f"Using Python {sys.version.split()[0]} from {sys.executable}")
    if not REQUIREMENTS.is_file():
        print(f"Dependency manifest is missing: {REQUIREMENTS}")
        return 2
    if args.repair and VENV.exists():
        shutil.rmtree(VENV)
    if not venv_python().exists():
        print("Creating an isolated Python environment inside the application folder…")
        venv.EnvBuilder(with_pip=True, clear=False).create(VENV)
    python = str(venv_python())
    signature = dependency_signature()
    installed_signature = MARKER.read_text(encoding="utf-8").strip() if MARKER.exists() else ""
    if installed_signature != signature:
        print("Installing the GUI, analysis, and Excel dependencies…")
        try:
            run([python, "-m", "pip", "install", "--upgrade", "pip"])
            run([python, "-m", "pip", "install", "--requirement", str(REQUIREMENTS)])
        except subprocess.CalledProcessError as exc:
            print(f"Dependency installation failed with exit code {exc.returncode}.")
            print("Check the internet connection, then rerun this launcher. Use --repair if the local environment is incomplete.")
            return exc.returncode or 1
        MARKER.write_text(signature, encoding="utf-8")
    print("Launching LC–MS Compound Curation Workbench…")
    return subprocess.call([python, "-m", "lcms_curation"], cwd=ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
