from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PortablePackageValidation(unittest.TestCase):
    def test_runtime_dependencies_are_pinned(self) -> None:
        requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
        names = {line.partition("==")[0]: line.partition("==")[2] for line in requirements if "==" in line}
        self.assertEqual(set(names), {"PySide6", "pandas", "numpy", "openpyxl"})
        for package, version in names.items():
            with self.subTest(package=package):
                self.assertRegex(version, r"^\d+\.\d+\.\d+$")

    def test_windows_launcher_bootstraps_the_local_environment(self) -> None:
        batch = (ROOT / "RUN_WINDOWS.bat").read_text(encoding="utf-8")
        installer = (ROOT / "install_and_run.py").read_text(encoding="utf-8")
        self.assertIn("install_and_run.py", batch)
        self.assertIn('VENV = ROOT / ".venv"', installer)
        self.assertIn("venv.EnvBuilder(with_pip=True", installer)
        self.assertIn('"--requirement", str(REQUIREMENTS)', installer)
        self.assertIn('subprocess.call([python, "-m", "lcms_curation"]', installer)

    def test_release_version_is_consistent(self) -> None:
        package_text = (ROOT / "lcms_curation" / "__init__.py").read_text(encoding="utf-8")
        gui_text = (ROOT / "lcms_curation" / "gui.py").read_text(encoding="utf-8")
        project_text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        package_version = re.search(r'__version__ = "([^"]+)"', package_text).group(1)
        gui_version = re.search(r'APP_VERSION = "([^"]+)"', gui_text).group(1)
        project_version = re.search(r'^version = "([^"]+)"', project_text, re.MULTILINE).group(1)
        self.assertEqual(package_version, gui_version)
        self.assertEqual(package_version, project_version)

    def test_v2_workflow_and_modules_are_packaged(self) -> None:
        gui = (ROOT / "lcms_curation" / "gui.py").read_text(encoding="utf-8")
        for label in (
            "Project + files",
            "Sample map",
            "Structure + MS²",
            "Identity audit",
            "Read resolution",
            "QC + filter",
            "Context model",
            "Export",
        ):
            self.assertIn(label, gui)
        for module in ("metadata.py", "qc.py", "chemistry.py", "provenance.py", "reporting.py"):
            self.assertTrue((ROOT / "lcms_curation" / module).is_file())
        self.assertIn('set "QT_STYLE_OVERRIDE=Fusion"', (ROOT / "RUN_WINDOWS.bat").read_text(encoding="utf-8"))

    def test_background_progress_notification_and_journal_report_contract(self) -> None:
        gui = (ROOT / "lcms_curation" / "gui.py").read_text(encoding="utf-8")
        exports = (ROOT / "lcms_curation" / "exports.py").read_text(encoding="utf-8")
        theme = (ROOT / "lcms_curation" / "theme.py").read_text(encoding="utf-8")
        for marker in (
            "BACKGROUND PIPELINE",
            "Background task and notification audit",
            "Generate / refresh report",
            "_run_background_task",
            "_show_task_notification",
            "journal_report_markdown",
        ):
            self.assertIn(marker, gui)
        for output in (
            "14_Task_History.csv",
            "JOURNAL_READY_METHODS_AND_ANALYSIS.md",
            '"14 Journal Methods"',
            '"15 Task History"',
        ):
            self.assertIn(output, exports)
        self.assertIn("QTextBrowser#journalView", theme)
        self.assertIn("QLabel#taskNotification", theme)


if __name__ == "__main__":
    unittest.main()
