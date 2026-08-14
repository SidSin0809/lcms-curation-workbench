from __future__ import annotations

from typing import Any


COLORS = {
    "window": "#F5F7FA",
    "base": "#FFFFFF",
    "alternate": "#F6FAFC",
    "ink": "#182B47",
    "navy": "#13284A",
    "muted": "#657188",
    "disabled": "#7B8799",
    "border": "#D9E0E9",
    "control_border": "#C9D3E0",
    "blue": "#1F5AE0",
    "blue_hover": "#194CC0",
    "pale_blue": "#DDE9FF",
    "aqua": "#2BC2B0",
    "success": "#19735F",
    "warning": "#A15B11",
    "danger": "#B42318",
}


def build_stylesheet() -> str:
    """Return the application stylesheet without importing Qt.

    Keeping this function Qt-free lets the Windows color contract be tested on
    build machines that do not have a graphical Qt runtime.
    """

    return """
        QWidget {
            color: #182B47;
            background-color: transparent;
        }
        QWidget:disabled { color: #7B8799; }
        QMainWindow, QDialog, QMessageBox, QWidget#page, QScrollArea,
        QScrollArea > QWidget > QWidget, QStackedWidget {
            background: #F5F7FA;
            color: #182B47;
        }
        QLabel {
            color: #182B47;
            background-color: transparent;
        }
        QFrame#appHeader { background: #13284A; }
        QLabel#brandBadge {
            background: #2BC2B0;
            color: #0C2444;
            font-size: 17px;
            font-weight: 800;
            border-radius: 7px;
            padding: 10px;
        }
        QLabel#brandTitle { color: #FFFFFF; font-size: 19px; font-weight: 700; }
        QLabel#brandSubtitle { color: #B9C9DF; font-size: 9px; letter-spacing: 1px; }

        QFrame#sidebar { background: #10213C; min-width: 250px; max-width: 320px; }
        QLabel#sidebarLabel { color: #8FA7C7; font-size: 10px; font-weight: 700; letter-spacing: 2px; }
        QLabel#sidebarNote {
            color: #B9C9DF;
            font-size: 10px;
            padding: 10px;
            background: #172E50;
            border-radius: 6px;
        }
        QListWidget#stageList {
            background: transparent;
            border: none;
            color: #CBD5E4;
            outline: none;
        }
        QListWidget#stageList::item { padding: 14px 10px; margin: 3px 0; border-radius: 6px; }
        QListWidget#stageList::item:selected {
            background: #1E416C;
            color: #FFFFFF;
            border-left: 3px solid #2BC2B0;
        }
        QListWidget#stageList::item:disabled { color: #7388A5; }

        QLabel#kicker { color: #1F5AE0; font-size: 10px; font-weight: 800; letter-spacing: 1.5px; }
        QLabel#sectionTitle { color: #13284A; font-size: 26px; font-weight: 750; }
        QLabel#bodyText { color: #657188; font-size: 12px; }
        QLabel#subheading { color: #13284A; font-size: 11px; font-weight: 800; letter-spacing: 1px; }

        QFrame#panel, QGroupBox {
            background: #FFFFFF;
            color: #182B47;
            border: 1px solid #D9E0E9;
            border-radius: 8px;
        }
        QFrame#panel QLabel, QGroupBox QLabel { color: #182B47; background: transparent; }
        QGroupBox { margin-top: 12px; padding-top: 12px; font-weight: 700; }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 12px;
            padding: 0 5px;
            color: #13284A;
            background: #FFFFFF;
        }
        QLabel#roleBadge {
            background: #DDE9FF;
            color: #1F5AE0;
            border-radius: 5px;
            font-weight: 800;
            padding: 6px 8px;
        }

        QLineEdit, QComboBox, QDoubleSpinBox {
            background: #FFFFFF;
            color: #182B47;
            border: 1px solid #C9D3E0;
            border-radius: 5px;
            padding: 7px 8px;
            min-height: 20px;
            selection-background-color: #DDE9FF;
            selection-color: #13284A;
        }
        QLineEdit { placeholder-text-color: #657188; }
        QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus { border: 1px solid #1F5AE0; }
        QLineEdit:disabled, QComboBox:disabled, QDoubleSpinBox:disabled {
            background: #F2F4F7;
            color: #7B8799;
            border-color: #D9E0E9;
        }
        QComboBox::drop-down, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {
            background: #EEF3FB;
            border: none;
            border-left: 1px solid #C9D3E0;
            width: 22px;
        }
        QComboBox QAbstractItemView {
            background: #FFFFFF;
            color: #182B47;
            border: 1px solid #9FB6E2;
            outline: none;
            selection-background-color: #DDE9FF;
            selection-color: #13284A;
        }

        QCheckBox { color: #182B47; background: transparent; spacing: 7px; }
        QCheckBox:disabled { color: #7B8799; }

        QPushButton {
            background: #FFFFFF;
            color: #182B47;
            border: 1px solid #C9D3E0;
            border-radius: 5px;
            padding: 8px 13px;
            font-weight: 650;
        }
        QPushButton:hover { background: #EEF3FB; }
        QPushButton:pressed { background: #DDE9FF; }
        QPushButton:disabled { background: #F2F4F7; color: #7B8799; border-color: #D9E0E9; }
        QPushButton#primaryButton { background: #1F5AE0; color: #FFFFFF; border: 1px solid #1F5AE0; }
        QPushButton#primaryButton:hover { background: #194CC0; }
        QPushButton#primaryButton:pressed { background: #133FA6; }
        QPushButton#secondaryButton, QPushButton#presetButton {
            background: #FFFFFF;
            color: #1F5AE0;
            border: 1px solid #9FB6E2;
        }
        QPushButton#secondaryButton:hover, QPushButton#presetButton:hover { background: #EAF1FF; }
        QFrame#appHeader QPushButton#ghostButton {
            background: transparent;
            color: #D7E3F3;
            border: 1px solid #627DA1;
        }
        QFrame#appHeader QPushButton#ghostButton:hover { background: #1E416C; color: #FFFFFF; }
        QWidget#page QPushButton#ghostButton {
            background: transparent;
            color: #52627A;
            border: 1px solid #AEBAC9;
        }
        QWidget#page QPushButton#ghostButton:hover { background: #EEF3FB; color: #13284A; }
        QPushButton#smallButton {
            background: #EAF1FF;
            color: #1F5AE0;
            border: 1px solid #C9D8F4;
            padding: 6px 10px;
        }

        QFrame#metricCard {
            background: #FFFFFF;
            border: 1px solid #D9E0E9;
            border-top: 3px solid #2BC2B0;
            border-radius: 8px;
        }
        QLabel#metricTitle { color: #657188; font-size: 9px; font-weight: 750; letter-spacing: 1px; }
        QLabel#metricValue { color: #13284A; font-size: 23px; font-weight: 800; }
        QLabel#metricNote { color: #657188; font-size: 9px; }
        QFrame#callout { background: #EAF8F7; border: 1px solid #B8E5DF; border-radius: 7px; }
        QFrame#warningBox { background: #FFF8E9; border: 1px solid #F0D79C; border-radius: 7px; }
        QFrame#callout QLabel, QFrame#warningBox QLabel { color: #182B47; background: transparent; }
        QLabel#calloutMarker {
            background: #2BC2B0;
            color: #FFFFFF;
            font-weight: 800;
            border-radius: 9px;
            min-width: 18px;
            max-width: 18px;
            min-height: 18px;
            max-height: 18px;
            qproperty-alignment: AlignCenter;
        }
        QLabel#methodBox { background: #EEF3FB; color: #182B47; border-left: 4px solid #1F5AE0; padding: 14px; }
        QLabel#resultLabel { color: #13284A; font-size: 15px; font-weight: 750; }
        QLabel#ciLabel { color: #19735F; font-weight: 700; }

        QTableView {
            background: #FFFFFF;
            alternate-background-color: #F6FAFC;
            color: #182B47;
            border: 1px solid #D9E0E9;
            border-radius: 5px;
            gridline-color: #E5EAF0;
            selection-background-color: #DDE9FF;
            selection-color: #13284A;
        }
        QTableView::item { color: #182B47; padding: 4px; }
        QTableView::item:selected { background: #DDE9FF; color: #13284A; }
        QHeaderView { background: #13284A; color: #FFFFFF; }
        QHeaderView::section {
            background: #13284A;
            color: #FFFFFF;
            border: none;
            border-right: 1px solid #314B70;
            padding: 7px;
            font-weight: 700;
        }
        QTableCornerButton::section { background: #13284A; border: 1px solid #314B70; }

        QProgressBar {
            background: #E6EBF1;
            color: #13284A;
            border: none;
            border-radius: 5px;
            text-align: center;
            min-height: 18px;
        }
        QProgressBar::chunk { background: #2BC2B0; border-radius: 5px; }
        QFrame#taskCenter {
            background: #FFFFFF;
            color: #182B47;
            border-top: 1px solid #C9D3E0;
        }
        QFrame#taskCenter QLabel { background: transparent; color: #182B47; }
        QLabel#taskStage {
            color: #1F5AE0;
            font-size: 9px;
            font-weight: 800;
            letter-spacing: 1px;
            min-width: 170px;
        }
        QLabel#taskMessage { color: #52627A; min-width: 360px; }
        QLabel#taskNotification {
            color: #146B5A;
            background: #DDF5F2;
            border: 1px solid #9ED7CF;
            border-radius: 5px;
            padding: 6px 10px;
            font-weight: 700;
        }
        QLabel#taskNotification[status="failed"] {
            color: #A61B13;
            background: #FDE8E7;
            border-color: #E8AAA5;
        }
        QTextBrowser#journalView {
            background: #FFFFFF;
            color: #182B47;
            border: 1px solid #C9D3E0;
            border-radius: 5px;
            padding: 12px;
            selection-background-color: #DDE9FF;
            selection-color: #13284A;
        }
        QTextBrowser#journalView a { color: #1F5AE0; }
        QStatusBar { background: #F5F7FA; color: #52627A; border-top: 1px solid #D9E0E9; }
        QStatusBar QLabel { color: #52627A; }
        QToolTip { background: #13284A; color: #FFFFFF; border: 1px solid #496486; padding: 5px; }
        QSplitter::handle { background: #D9E0E9; }
        QSplitter::handle:hover { background: #9FB6E2; }
        QScrollBar:vertical, QScrollBar:horizontal { background: #EEF1F5; border: none; }
        QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
            background: #9CAAC0;
            border-radius: 5px;
            min-height: 24px;
            min-width: 24px;
        }
        QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
    """


def apply_accessible_light_theme(application: Any) -> None:
    """Apply a platform-independent Fusion style, palette, and stylesheet."""

    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtWidgets import QStyleFactory

    fusion = QStyleFactory.create("Fusion")
    if fusion is not None:
        application.setStyle(fusion)
    palette = application.style().standardPalette()
    all_group_roles = {
        QPalette.ColorRole.Window: COLORS["window"],
        QPalette.ColorRole.WindowText: COLORS["ink"],
        QPalette.ColorRole.Base: COLORS["base"],
        QPalette.ColorRole.AlternateBase: COLORS["alternate"],
        QPalette.ColorRole.ToolTipBase: COLORS["navy"],
        QPalette.ColorRole.ToolTipText: COLORS["base"],
        QPalette.ColorRole.Text: COLORS["ink"],
        QPalette.ColorRole.Button: COLORS["base"],
        QPalette.ColorRole.ButtonText: COLORS["ink"],
        QPalette.ColorRole.BrightText: COLORS["base"],
        QPalette.ColorRole.Link: COLORS["blue"],
        QPalette.ColorRole.LinkVisited: "#5A3EB5",
        QPalette.ColorRole.Highlight: COLORS["pale_blue"],
        QPalette.ColorRole.HighlightedText: COLORS["navy"],
        QPalette.ColorRole.PlaceholderText: COLORS["muted"],
        QPalette.ColorRole.Light: COLORS["base"],
        QPalette.ColorRole.Midlight: "#EEF1F5",
        QPalette.ColorRole.Mid: "#B8C2CF",
        QPalette.ColorRole.Dark: "#657188",
        QPalette.ColorRole.Shadow: "#31425A",
    }
    for role, color in all_group_roles.items():
        palette.setColor(role, QColor(color))
    for role in (QPalette.ColorRole.WindowText, QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText):
        palette.setColor(QPalette.ColorGroup.Disabled, role, QColor(COLORS["disabled"]))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Base, QColor("#F2F4F7"))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Button, QColor("#F2F4F7"))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Highlight, QColor("#E1E6ED"))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.HighlightedText, QColor(COLORS["disabled"]))
    application.setPalette(palette)
    application.setStyleSheet(build_stylesheet())
