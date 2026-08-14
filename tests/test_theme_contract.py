from __future__ import annotations

import unittest

from lcms_curation.theme import COLORS, build_stylesheet


def relative_luminance(color: str) -> float:
    channels = [int(color[index : index + 2], 16) / 255 for index in (1, 3, 5)]
    linear = [channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4 for channel in channels]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast_ratio(foreground: str, background: str) -> float:
    light, dark = sorted((relative_luminance(foreground), relative_luminance(background)), reverse=True)
    return (light + 0.05) / (dark + 0.05)


class ThemeContractValidation(unittest.TestCase):
    def test_normal_text_pairs_have_accessible_contrast(self) -> None:
        pairs = (
            (COLORS["ink"], COLORS["window"]),
            (COLORS["ink"], COLORS["base"]),
            (COLORS["muted"], COLORS["base"]),
            (COLORS["base"], COLORS["navy"]),
            (COLORS["blue"], COLORS["pale_blue"]),
        )
        for foreground, background in pairs:
            with self.subTest(foreground=foreground, background=background):
                self.assertGreaterEqual(contrast_ratio(foreground, background), 4.5)

    def test_every_problem_widget_has_an_explicit_foreground_contract(self) -> None:
        stylesheet = " ".join(build_stylesheet().split())
        expected_rules = (
            "QLabel { color: #182B47;",
            "QFrame#panel QLabel, QGroupBox QLabel { color: #182B47;",
            "QLineEdit, QComboBox, QDoubleSpinBox { background: #FFFFFF; color: #182B47;",
            "QComboBox QAbstractItemView { background: #FFFFFF; color: #182B47;",
            "QCheckBox { color: #182B47;",
            "QTableView { background: #FFFFFF; alternate-background-color: #F6FAFC; color: #182B47;",
            "QProgressBar { background: #E6EBF1; color: #13284A;",
            "QStatusBar { background: #F5F7FA; color: #52627A;",
        )
        for rule in expected_rules:
            with self.subTest(rule=rule):
                self.assertIn(rule, stylesheet)

    def test_header_and_page_ghost_buttons_use_different_readable_colors(self) -> None:
        stylesheet = build_stylesheet()
        self.assertIn("QFrame#appHeader QPushButton#ghostButton", stylesheet)
        self.assertIn("QWidget#page QPushButton#ghostButton", stylesheet)
        rules = [line.strip() for line in stylesheet.splitlines()]
        self.assertNotIn("QPushButton#ghostButton {", rules)


if __name__ == "__main__":
    unittest.main()
