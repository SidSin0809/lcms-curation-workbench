from __future__ import annotations

import html
import json
import os
import re
import sys
import time
import traceback
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSettings, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QStyledItemDelegate,
    QTableView,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from .engine import (
    FILE_ROLES,
    REQUIRED_FIELDS,
    ROLE_LABELS,
    AnalysisBundle,
    analyze_files,
    apply_filters,
    evaluate_filters,
    rebuild_with_sample_metadata,
)
from .exports import export_complete_folder, export_csv, export_workbook, threshold_table
from .metadata import (
    SAMPLE_ROLES,
    build_sample_metadata,
    merge_sample_metadata,
    metadata_requirement_table,
    suggested_results_name,
)
from .provenance import score_provenance
from .qc import default_qc_settings
from .reporting import journal_report_markdown
from .theme import apply_accessible_light_theme

APP_NAME = "LC–MS Compound Curation Workbench"
APP_VERSION = "2.1.0"


EXTRACTION_METHODS = [
    ("Unknown / not declared", "unknown"),
    ("MTBE biphasic LLE", "mtbe"),
    ("Folch chloroform/methanol LLE", "folch"),
    ("Bligh–Dyer LLE", "bligh-dyer"),
    ("Generic liquid–liquid extraction", "generic-lle"),
    ("Ethyl acetate LLE", "ethyl-acetate"),
    ("Butanol LLE", "butanol"),
    ("Protein precipitation", "protein-precipitation"),
    ("SPE · HLB", "spe-hlb"),
    ("SPE · C18 / reversed phase", "spe-c18"),
    ("SPE · HILIC", "spe-hilic"),
    ("SPE · cation exchange", "spe-cation"),
    ("SPE · anion exchange", "spe-anion"),
    ("SPE · mixed mode", "spe-mixed-mode"),
    ("Monophasic polar extraction", "monophasic-polar"),
    ("Other / custom", "other"),
]

ANALYZED_PHASES = [
    ("Unknown / not declared", "unknown"),
    ("Upper organic", "upper-organic"),
    ("Lower organic", "lower-organic"),
    ("Aqueous-rich phase", "aqueous"),
    ("Upper aqueous", "upper-aqueous"),
    ("Lower aqueous", "lower-aqueous"),
    ("Supernatant", "supernatant"),
    ("SPE eluate / retained fraction", "spe-eluate"),
    ("SPE flow-through", "spe-flowthrough"),
    ("Interphase / pellet", "interphase-pellet"),
    ("Both organic and aqueous", "both"),
]

BIOLOGICAL_SYSTEMS = [
    ("Human biofluid / tissue", "human"),
    ("Human + microbiome", "human-microbiome"),
    ("Mammalian cell line", "cell-line"),
    ("Bacterial / microbial culture", "bacterial-culture"),
    ("Virus-exposure model", "virus-exposure"),
    ("Environmental / non-human", "environmental"),
    ("Other / mixed", "other"),
]

SAMPLE_MATRICES = [
    "Serum",
    "Plasma",
    "Saliva",
    "Urine",
    "Human tissue",
    "Placenta",
    "Cell line",
    "Bacterial culture",
    "Virus-exposure cell culture",
    "Feces / stool",
    "Environmental sample",
    "Other",
]


class DataFrameModel(QAbstractTableModel):
    def __init__(self, frame: pd.DataFrame | None = None, parent: QWidget | None = None):
        super().__init__(parent)
        self._frame = (frame if frame is not None else pd.DataFrame()).reset_index(drop=True)

    def set_frame(self, frame: pd.DataFrame) -> None:
        self.beginResetModel()
        self._frame = frame.reset_index(drop=True).copy()
        self.endResetModel()

    def frame(self) -> pd.DataFrame:
        return self._frame.copy()

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._frame)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._frame.columns)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if not index.isValid():
            return None
        value = self._frame.iat[index.row(), index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            if value is None or (isinstance(value, float) and np.isnan(value)):
                return ""
            if isinstance(value, (float, np.floating)):
                absolute = abs(float(value))
                if absolute and (absolute < 0.001 or absolute >= 1_000_000):
                    return f"{float(value):.4g}"
                return f"{float(value):.4f}".rstrip("0").rstrip(".")
            if isinstance(value, (bool, np.bool_)):
                return "Yes" if value else "No"
            return str(value)
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(value, (int, float, np.number)):
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.ToolTipRole:
            return "" if value is None else str(value)
        if role == Qt.ItemDataRole.ForegroundRole:
            column = str(self._frame.columns[index.column()])
            if column in {"Status", "Filter decision", "Analytical QC decision"}:
                return QColor("#19735f" if str(value) == "Pass" else "#a15b11" if str(value) in {"Review", "Not evaluable"} else "#b42318")
        return None

    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if orientation == Qt.Orientation.Horizontal and 0 <= section < len(self._frame.columns):
            if role in {Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.ToolTipRole}:
                return str(self._frame.columns[section])
        if orientation == Qt.Orientation.Vertical and role == Qt.ItemDataRole.DisplayRole:
            return str(section + 1)
        return None

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        if not 0 <= column < len(self._frame.columns):
            return
        self.layoutAboutToBeChanged.emit()
        name = self._frame.columns[column]
        self._frame = self._frame.sort_values(name, ascending=order == Qt.SortOrder.AscendingOrder, kind="mergesort", na_position="last").reset_index(drop=True)
        self.layoutChanged.emit()


class EditableDataFrameModel(DataFrameModel):
    def __init__(self, frame: pd.DataFrame | None = None, editable: set[str] | None = None, parent: QWidget | None = None):
        super().__init__(frame, parent)
        self.editable = editable or set()

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        flags = super().flags(index)
        if not index.isValid():
            return flags
        column = str(self._frame.columns[index.column()])
        if column in self.editable:
            flags |= Qt.ItemFlag.ItemIsEditable
            if column in {"Include", "Technical replicate"}:
                flags |= Qt.ItemFlag.ItemIsUserCheckable
        return flags

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        if index.isValid() and str(self._frame.columns[index.column()]) in {"Include", "Technical replicate"}:
            if role == Qt.ItemDataRole.CheckStateRole:
                return Qt.CheckState.Checked if bool(self._frame.iat[index.row(), index.column()]) else Qt.CheckState.Unchecked
            if role == Qt.ItemDataRole.DisplayRole:
                return ""
        return super().data(index, role)

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.EditRole) -> bool:
        if not index.isValid():
            return False
        column = str(self._frame.columns[index.column()])
        if column not in self.editable:
            return False
        if role == Qt.ItemDataRole.CheckStateRole and column in {"Include", "Technical replicate"}:
            self._frame.iat[index.row(), index.column()] = value == Qt.CheckState.Checked.value or value == Qt.CheckState.Checked
        elif role == Qt.ItemDataRole.EditRole:
            if column in {"Injection order", "Dilution factor"}:
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    return False
            self._frame.iat[index.row(), index.column()] = value
        else:
            return False
        self.dataChanged.emit(index, index, [role])
        return True


class ComboBoxDelegate(QStyledItemDelegate):
    def __init__(self, choices: list[str] | tuple[str, ...], parent: QWidget | None = None):
        super().__init__(parent)
        self.choices = list(choices)

    def createEditor(self, parent: QWidget, _option: Any, _index: QModelIndex) -> QWidget:
        editor = QComboBox(parent)
        editor.addItems(self.choices)
        return editor

    def setEditorData(self, editor: QWidget, index: QModelIndex) -> None:
        if isinstance(editor, QComboBox):
            editor.setCurrentText(str(index.data(Qt.ItemDataRole.DisplayRole) or ""))

    def setModelData(self, editor: QWidget, model: QAbstractTableModel, index: QModelIndex) -> None:
        if isinstance(editor, QComboBox):
            model.setData(index, editor.currentText(), Qt.ItemDataRole.EditRole)


class AnalysisWorker(QThread):
    progressed = Signal(int, str)
    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(self, files: dict[str, str], parameters: dict[str, Any], parent: QWidget | None = None):
        super().__init__(parent)
        self.files = files
        self.parameters = parameters

    def run(self) -> None:
        try:
            result = analyze_files(self.files, self.parameters, lambda percent, text: self.progressed.emit(percent, text))
            self.succeeded.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())


class FunctionWorker(QThread):
    progressed = Signal(int, str)
    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(self, function: Callable[[Callable[[int, str], None]], Any], parent: QWidget | None = None):
        super().__init__(parent)
        self.function = function

    def run(self) -> None:
        try:
            self.succeeded.emit(self.function(lambda percent, text: self.progressed.emit(percent, text)))
        except Exception:
            self.failed.emit(traceback.format_exc())


class MetricCard(QFrame):
    def __init__(self, title: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("metricCard")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        title_label = QLabel(title.upper())
        title_label.setObjectName("metricTitle")
        self.value = QLabel("—")
        self.value.setObjectName("metricValue")
        self.note = QLabel("")
        self.note.setObjectName("metricNote")
        self.note.setWordWrap(True)
        layout.addWidget(title_label)
        layout.addWidget(self.value)
        layout.addWidget(self.note)

    def set_metric(self, value: str, note: str = "") -> None:
        self.value.setText(value)
        self.note.setText(note)


def _section_heading(kicker: str, title: str, body: str = "") -> QWidget:
    widget = QWidget()
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 8)
    label = QLabel(kicker.upper())
    label.setObjectName("kicker")
    heading = QLabel(title)
    heading.setObjectName("sectionTitle")
    heading.setWordWrap(True)
    layout.addWidget(label)
    layout.addWidget(heading)
    if body:
        paragraph = QLabel(body)
        paragraph.setObjectName("bodyText")
        paragraph.setWordWrap(True)
        layout.addWidget(paragraph)
    return widget


def _callout(text: str, kind: str = "info") -> QFrame:
    frame = QFrame()
    frame.setObjectName("warningBox" if kind == "warning" else "callout")
    layout = QHBoxLayout(frame)
    layout.setContentsMargins(14, 12, 14, 12)
    marker = QLabel("!" if kind == "warning" else "i")
    marker.setObjectName("calloutMarker")
    label = QLabel(text)
    label.setWordWrap(True)
    label.setTextFormat(Qt.TextFormat.RichText)
    layout.addWidget(marker, 0, Qt.AlignmentFlag.AlignTop)
    layout.addWidget(label, 1)
    return frame


def _table(minimum_height: int = 300) -> tuple[QTableView, DataFrameModel]:
    table = QTableView()
    model = DataFrameModel()
    table.setModel(model)
    table.setAlternatingRowColors(True)
    table.setSortingEnabled(True)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
    table.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
    table.verticalHeader().setDefaultSectionSize(30)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    table.horizontalHeader().setSectionsMovable(True)
    table.horizontalHeader().setStretchLastSection(False)
    table.setWordWrap(False)
    table.setMinimumHeight(minimum_height)
    return table, model


def _fit_table(table: QTableView, maximum: int = 240) -> None:
    table.resizeColumnsToContents()
    for column in range(table.model().columnCount()):
        table.setColumnWidth(column, min(maximum, max(80, table.columnWidth(column))))


def _add_combo_items(combo: QComboBox, items: list[tuple[str, str]]) -> None:
    for label, value in items:
        combo.addItem(label, value)


def _set_combo_data(combo: QComboBox, value: Any) -> None:
    index = combo.findData(value)
    if index >= 0:
        combo.setCurrentIndex(index)


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"{APP_NAME} · {APP_VERSION}")
        self.resize(1540, 960)
        self.setMinimumSize(1160, 760)
        self.settings = QSettings("SiddharthSinghResearch", "LCMSCurationWorkbenchV2")
        self.analysis: AnalysisBundle | None = None
        self.shortlist = pd.DataFrame()
        self.filter_decisions = pd.DataFrame()
        self.provenance = pd.DataFrame()
        self.filter_settings: dict[str, Any] = {}
        self.provenance_settings: dict[str, Any] | None = None
        self.journal_report = ""
        self.pending_sample_metadata: pd.DataFrame | None = None
        self._workers: list[QThread] = []
        self._active_task: dict[str, Any] | None = None
        self.task_history: list[dict[str, Any]] = []
        self._max_stage = 0
        self._build_ui()

    def _build_ui(self) -> None:
        root = QWidget()
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)
        root_layout.addWidget(self._build_header())
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self._build_sidebar())
        self.stack = QStackedWidget()
        builders = (
            self._build_input_page,
            self._build_sample_page,
            self._build_structure_page,
            self._build_identity_page,
            self._build_selection_page,
            self._build_qc_filter_page,
            self._build_context_page,
            self._build_export_page,
        )
        for builder in builders:
            self.stack.addWidget(self._scroll_page(builder()))
        splitter.addWidget(self.stack)
        splitter.setSizes([300, 1240])
        root_layout.addWidget(splitter, 1)
        root_layout.addWidget(self._build_task_center())
        self.setCentralWidget(root)
        self.statusBar().showMessage("Create or load a project, then select all six files.")

    def _build_task_center(self) -> QWidget:
        center = QFrame()
        center.setObjectName("taskCenter")
        layout = QVBoxLayout(center)
        layout.setContentsMargins(18, 8, 18, 8)
        layout.setSpacing(5)
        row = QHBoxLayout()
        self.task_stage = QLabel("BACKGROUND PIPELINE")
        self.task_stage.setObjectName("taskStage")
        self.task_message = QLabel("Ready — completed tasks and failures will be reported here.")
        self.task_message.setObjectName("taskMessage")
        self.global_progress = QProgressBar()
        self.global_progress.setObjectName("globalProgress")
        self.global_progress.setRange(0, 100)
        self.global_progress.setValue(0)
        self.global_progress.setFormat("Ready")
        row.addWidget(self.task_stage)
        row.addWidget(self.global_progress, 1)
        row.addWidget(self.task_message, 2)
        layout.addLayout(row)
        self.task_notification = QLabel("")
        self.task_notification.setObjectName("taskNotification")
        self.task_notification.setWordWrap(True)
        self.task_notification.hide()
        layout.addWidget(self.task_notification)
        self.notification_timer = QTimer(self)
        self.notification_timer.setSingleShot(True)
        self.notification_timer.timeout.connect(self.task_notification.hide)
        return center

    def _build_header(self) -> QWidget:
        header = QFrame()
        header.setObjectName("appHeader")
        layout = QHBoxLayout(header)
        layout.setContentsMargins(24, 14, 24, 14)
        badge = QLabel("LC")
        badge.setObjectName("brandBadge")
        title_box = QVBoxLayout()
        title = QLabel("Compound Curation")
        title.setObjectName("brandTitle")
        subtitle = QLabel("PROJECT · SAMPLE MAP · SIX-FILE EVIDENCE · MS² · QC · CONTEXT · EXPORT")
        subtitle.setObjectName("brandSubtitle")
        title_box.addWidget(title)
        title_box.addWidget(subtitle)
        layout.addWidget(badge)
        layout.addLayout(title_box)
        layout.addStretch(1)
        load = QPushButton("Load project")
        load.setObjectName("ghostButton")
        load.clicked.connect(self._load_project)
        save = QPushButton("Save project")
        save.setObjectName("ghostButton")
        save.clicked.connect(self._save_project)
        about = QPushButton("Method & version")
        about.setObjectName("ghostButton")
        about.clicked.connect(self._show_about)
        layout.addWidget(load)
        layout.addWidget(save)
        layout.addWidget(about)
        return header

    def _build_sidebar(self) -> QWidget:
        sidebar = QFrame()
        sidebar.setObjectName("sidebar")
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(16, 22, 16, 18)
        label = QLabel("WORKFLOW")
        label.setObjectName("sidebarLabel")
        layout.addWidget(label)
        self.stage_list = QListWidget()
        self.stage_list.setObjectName("stageList")
        stages = [
            ("01  Project + files", "Six required exports and method context"),
            ("02  Sample map", "Roles, groups, batch, order, dilution"),
            ("03  Structure + MS²", "Join audit, spectra, fragment peaks"),
            ("04  Identity audit", "Accepted identities and repeated reads"),
            ("05  Read resolution", "Deterministic representative ranking"),
            ("06  QC + filter", "Analytical QC and adaptive thresholds"),
            ("07  Context model", "Chemistry, source, extraction phase"),
            ("08  Export", "CSV, XLSX, manifest, method notes"),
        ]
        for index, (title, subtitle) in enumerate(stages):
            item = QListWidgetItem(f"{title}\n     {subtitle}")
            item.setData(Qt.ItemDataRole.UserRole, index)
            if index > 0:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.stage_list.addItem(item)
        self.stage_list.setCurrentRow(0)
        self.stage_list.currentRowChanged.connect(self._navigate)
        layout.addWidget(self.stage_list, 1)
        privacy = QLabel("Runs locally. The application contains no upload or web-call code. Project files store paths and settings only.")
        privacy.setObjectName("sidebarNote")
        privacy.setWordWrap(True)
        layout.addWidget(privacy)
        return sidebar

    def _scroll_page(self, content: QWidget) -> QScrollArea:
        area = QScrollArea()
        area.setWidgetResizable(True)
        area.setFrameShape(QFrame.Shape.NoFrame)
        area.setWidget(content)
        return area

    def _page(self) -> tuple[QWidget, QVBoxLayout]:
        page = QWidget()
        page.setObjectName("page")
        layout = QVBoxLayout(page)
        layout.setContentsMargins(30, 26, 30, 34)
        layout.setSpacing(18)
        return page, layout

    def _begin_task(self, stage: str, task: str) -> bool:
        if self._active_task is not None:
            QMessageBox.information(
                self,
                "Background task in progress",
                f"Finish the current task first:\n{self._active_task['Task']}",
            )
            return False
        self._active_task = {
            "Stage": stage,
            "Task": task,
            "Status": "Running",
            "Started": datetime.now().astimezone().isoformat(timespec="seconds"),
            "Finished": "",
            "Duration (s)": np.nan,
            "Outcome": "In progress",
            "_started_perf": time.perf_counter(),
        }
        self.task_stage.setText(stage.upper())
        self.task_message.setText(task)
        self.global_progress.setValue(0)
        self.global_progress.setFormat("0%")
        self.task_notification.hide()
        self.notification_timer.stop()
        self.statusBar().showMessage(f"{stage} · {task}")
        self._refresh_task_history()
        return True

    def _task_progress(self, percent: int, message: str) -> None:
        value = max(0, min(100, int(percent)))
        self.global_progress.setValue(value)
        self.global_progress.setFormat(f"{value}%")
        self.task_message.setText(message)
        self.statusBar().showMessage(message)

    def _finish_task(self, status: str, outcome: str) -> None:
        if self._active_task is None:
            return
        entry = dict(self._active_task)
        entry["Status"] = status
        entry["Finished"] = datetime.now().astimezone().isoformat(timespec="seconds")
        entry["Duration (s)"] = round(time.perf_counter() - float(entry.pop("_started_perf")), 3)
        entry["Outcome"] = outcome
        self.task_history.append(entry)
        self._active_task = None
        self.global_progress.setValue(100 if status == "Completed" else 0)
        self.global_progress.setFormat(status)
        self.task_message.setText(outcome)
        self.statusBar().showMessage(outcome)
        self._show_task_notification(status, f"{entry['Task']}: {outcome}")
        self._refresh_task_history()

    def _complete_task(self, outcome: str) -> None:
        self._finish_task("Completed", outcome)

    def _fail_task(self, outcome: str) -> None:
        self._finish_task("Failed", outcome)

    def _show_task_notification(self, status: str, message: str) -> None:
        self.task_notification.setProperty("status", status.casefold())
        self.task_notification.setText(("✓ " if status == "Completed" else "! ") + message)
        self.task_notification.style().unpolish(self.task_notification)
        self.task_notification.style().polish(self.task_notification)
        self.task_notification.show()
        self.notification_timer.start(10_000)
        QApplication.beep()

    def _task_history_frame(self, include_active: bool = False) -> pd.DataFrame:
        rows = [dict(row) for row in self.task_history]
        if include_active and self._active_task is not None:
            active = dict(self._active_task)
            active.pop("_started_perf", None)
            rows.append(active)
        columns = ["Stage", "Task", "Status", "Started", "Finished", "Duration (s)", "Outcome"]
        return pd.DataFrame(rows, columns=columns)

    def _refresh_task_history(self) -> None:
        if not hasattr(self, "task_history_model"):
            return
        self.task_history_model.set_frame(self._task_history_frame(include_active=True))
        _fit_table(self.task_history_table, 260)

    def _run_background_task(
        self,
        stage: str,
        task: str,
        function: Callable[[Callable[[int, str], None]], Any],
        on_success: Callable[[Any], None] | None = None,
        outcome: str | Callable[[Any], str] = "Task completed successfully.",
        failure_title: str = "Task failed",
    ) -> None:
        if not self._begin_task(stage, task):
            return
        worker = FunctionWorker(function, self)
        self._workers.append(worker)
        worker.progressed.connect(self._task_progress)

        def succeeded(result: Any) -> None:
            message = outcome(result) if callable(outcome) else outcome
            self._complete_task(message)
            if on_success is not None:
                on_success(result)

        def failed(trace: str) -> None:
            final_line = trace.strip().splitlines()[-1] if trace.strip() else "Unknown error"
            self._fail_task(final_line)
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Icon.Critical)
            box.setWindowTitle(failure_title)
            box.setText(final_line)
            box.setDetailedText(trace)
            box.exec()

        worker.succeeded.connect(succeeded)
        worker.failed.connect(failed)
        worker.finished.connect(lambda: self._forget_worker(worker))
        worker.start()

    def _record_instant_task(self, stage: str, task: str, outcome: str) -> None:
        if not self._begin_task(stage, task):
            return
        self._task_progress(100, outcome)
        self._complete_task(outcome)

    def _build_input_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 01 / Project and input integrity",
            "Define the experiment, then completely scan all six exports",
            "Method context and file provenance are captured before identity curation. The scanner retains raw and normalized intensities, every CI candidate, isotope nodes, spectra, and fragment peaks.",
        ))
        project = QGroupBox("Project and analytical context")
        grid = QGridLayout(project)
        grid.addWidget(QLabel("Project name"), 0, 0)
        self.project_name = QLineEdit(str(self.settings.value("project/name", "LCMS Curation Project")))
        grid.addWidget(self.project_name, 0, 1)
        grid.addWidget(QLabel("Analyst / laboratory"), 0, 2)
        self.analyst = QLineEdit(str(self.settings.value("project/analyst", "")))
        grid.addWidget(self.analyst, 0, 3)
        grid.addWidget(QLabel("Assay type"), 1, 0)
        self.assay_type = QComboBox()
        self.assay_type.addItem("Metabolomics", "metabolomics")
        self.assay_type.addItem("Lipidomics", "lipidomics")
        self.assay_type.addItem("Mixed / other", "mixed")
        grid.addWidget(self.assay_type, 1, 1)
        grid.addWidget(QLabel("Biological system"), 1, 2)
        self.biological_system = QComboBox()
        _add_combo_items(self.biological_system, BIOLOGICAL_SYSTEMS)
        grid.addWidget(self.biological_system, 1, 3)
        grid.addWidget(QLabel("Sample matrix"), 2, 0)
        self.sample_matrix = QComboBox()
        self.sample_matrix.setEditable(True)
        self.sample_matrix.addItems(SAMPLE_MATRICES)
        grid.addWidget(self.sample_matrix, 2, 1)
        grid.addWidget(QLabel("Extraction method"), 2, 2)
        self.extraction_method = QComboBox()
        _add_combo_items(self.extraction_method, EXTRACTION_METHODS)
        grid.addWidget(self.extraction_method, 2, 3)
        grid.addWidget(QLabel("Analyzed phase / fraction"), 3, 0)
        self.analyzed_phase = QComboBox()
        _add_combo_items(self.analyzed_phase, ANALYZED_PHASES)
        grid.addWidget(self.analyzed_phase, 3, 1)
        grid.addWidget(QLabel("QC pooling point"), 3, 2)
        self.qc_pooling = QComboBox()
        self.qc_pooling.addItem("Unknown / not declared", "unknown")
        self.qc_pooling.addItem("Pooled before extraction", "before-extraction")
        self.qc_pooling.addItem("Pooled after extraction", "after-extraction")
        grid.addWidget(self.qc_pooling, 3, 3)
        grid.addWidget(QLabel("Ion mode"), 4, 0)
        self.ion_mode = QComboBox()
        self.ion_mode.addItem("Negative", "negative")
        self.ion_mode.addItem("Positive", "positive")
        self.ion_mode.addItem("Mixed", "mixed")
        grid.addWidget(self.ion_mode, 4, 1)
        grid.addWidget(QLabel("MS acquisition"), 4, 2)
        self.acquisition_mode = QComboBox()
        for label, value in (
            ("Unknown / not declared", "unknown"), ("DDA MS/MS", "dda"), ("DIA / MSE", "dia-mse"),
            ("HDMSE", "hdmse"), ("All-ion fragmentation", "aif"), ("Targeted MS/MS", "targeted-ms2"), ("MS1 only", "ms1-only"),
        ):
            self.acquisition_mode.addItem(label, value)
        grid.addWidget(self.acquisition_mode, 4, 3)
        grid.addWidget(QLabel("Search tolerance (ppm)"), 5, 0)
        self.mass_tolerance = QDoubleSpinBox()
        self.mass_tolerance.setRange(0.1, 100.0)
        self.mass_tolerance.setDecimals(1)
        self.mass_tolerance.setValue(float(self.settings.value("parameters/mass_tolerance", 10.0)))
        grid.addWidget(self.mass_tolerance, 5, 1)
        grid.addWidget(QLabel("Pooled-QC name regex"), 5, 2)
        self.qc_pattern = QLineEdit(str(self.settings.value("parameters/qc_pattern", "QC|Pool")))
        self.qc_pattern.setPlaceholderText("e.g. QC|Pool")
        grid.addWidget(self.qc_pattern, 5, 3)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        layout.addWidget(project)

        toolbar = QHBoxLayout()
        select_all = QPushButton("Select all six files")
        select_all.setObjectName("secondaryButton")
        select_all.clicked.connect(self._select_all_files)
        clear = QPushButton("Clear paths")
        clear.setObjectName("ghostButton")
        clear.clicked.connect(self._clear_paths)
        toolbar.addWidget(select_all)
        toolbar.addWidget(clear)
        toolbar.addStretch(1)
        layout.addLayout(toolbar)

        file_panel = QFrame()
        file_panel.setObjectName("panel")
        file_grid = QGridLayout(file_panel)
        file_grid.setContentsMargins(18, 18, 18, 18)
        file_grid.addWidget(QLabel("ROLE"), 0, 0)
        file_grid.addWidget(QLabel("REQUIRED FILE"), 0, 1)
        file_grid.addWidget(QLabel("PATH"), 0, 2)
        self.path_edits: dict[str, QLineEdit] = {}
        for row, role in enumerate(FILE_ROLES, 1):
            role_label = QLabel(role)
            role_label.setObjectName("roleBadge")
            description = QLabel(ROLE_LABELS[role])
            description.setWordWrap(True)
            edit = QLineEdit(str(self.settings.value(f"path/{role}", "")))
            edit.setPlaceholderText(f"Select {role} file…")
            browse = QPushButton("Browse")
            browse.setObjectName("smallButton")
            browse.clicked.connect(lambda _checked=False, selected_role=role: self._browse_role(selected_role))
            path_row = QWidget()
            path_layout = QHBoxLayout(path_row)
            path_layout.setContentsMargins(0, 0, 0, 0)
            path_layout.addWidget(edit, 1)
            path_layout.addWidget(browse)
            file_grid.addWidget(role_label, row, 0)
            file_grid.addWidget(description, row, 1)
            file_grid.addWidget(path_row, row, 2)
            self.path_edits[role] = edit
        file_grid.setColumnStretch(2, 1)
        layout.addWidget(file_panel)

        requirements = QGroupBox("Minimum data contract and downstream evidence")
        req_layout = QVBoxLayout(requirements)
        for role in FILE_ROLES:
            label = QLabel(f"<b>{role}</b> · {html.escape(', '.join(REQUIRED_FIELDS[role]))}")
            label.setTextFormat(Qt.TextFormat.RichText)
            label.setWordWrap(True)
            req_layout.addWidget(label)
        req_layout.addWidget(_callout(
            "<b>Required for analytical QC:</b> pooled-QC roles and replicate injections. <b>Required for contamination assessment:</b> process, extraction, and/or solvent blanks. "
            "<b>Required for stronger identification claims:</b> authentic standards, target-decoy results, RT/CCS references, and documented spectral-library search settings.",
        ))
        layout.addWidget(requirements)

        progress_panel = QFrame()
        progress_panel.setObjectName("panel")
        progress_layout = QHBoxLayout(progress_panel)
        self.scan_progress = QProgressBar()
        self.scan_progress.setRange(0, 100)
        self.scan_status = QLabel("Ready")
        self.scan_button = QPushButton("Scan all six files")
        self.scan_button.setObjectName("primaryButton")
        self.scan_button.clicked.connect(self._start_analysis)
        progress_layout.addWidget(self.scan_progress, 1)
        progress_layout.addWidget(self.scan_status)
        progress_layout.addWidget(self.scan_button)
        layout.addWidget(progress_panel)
        return page

    def _build_sample_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 02 / Sample metadata",
            "Map pooled QC, study samples, blanks, batches, and replicates",
            "Automatic roles use the CM condition/sample labels and your QC regex. Review the table: role errors propagate into CV, drift, blank, detection, and representative tie-break calculations.",
        ))
        cards = QHBoxLayout()
        self.sample_cards = [MetricCard("Sample columns"), MetricCard("Biological"), MetricCard("Pooled QC"), MetricCard("Blanks")]
        for card in self.sample_cards:
            cards.addWidget(card)
        layout.addLayout(cards)
        layout.addWidget(_callout(
            "A pooled QC prepared <b>before extraction</b> contains preparation plus analytical variation; pooling <b>after extraction</b> mainly reflects analytical variation. "
            "The tool records this distinction but does not infer it from filenames.",
        ))
        buttons = QHBoxLayout()
        import_button = QPushButton("Import metadata CSV")
        import_button.setObjectName("secondaryButton")
        import_button.clicked.connect(self._import_sample_metadata)
        export_button = QPushButton("Export metadata template")
        export_button.setObjectName("secondaryButton")
        export_button.clicked.connect(lambda: self._export_stage_csv("sample_metadata"))
        reset_button = QPushButton("Re-run automatic role rules")
        reset_button.setObjectName("ghostButton")
        reset_button.clicked.connect(self._reset_sample_roles)
        buttons.addWidget(import_button)
        buttons.addWidget(export_button)
        buttons.addWidget(reset_button)
        buttons.addStretch(1)
        layout.addLayout(buttons)
        editable = {
            "Condition / group", "Sample role", "Include", "Subject / biological unit", "Technical replicate",
            "Batch", "Injection order", "Dilution factor", "Sample matrix", "Biological system", "Extraction method", "Analyzed phase", "Notes",
        }
        self.sample_model = EditableDataFrameModel(pd.DataFrame(), editable)
        self.sample_table = QTableView()
        self.sample_table.setModel(self.sample_model)
        self.sample_table.setAlternatingRowColors(True)
        self.sample_table.setSortingEnabled(False)
        self.sample_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.sample_table.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.sample_table.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked | QAbstractItemView.EditTrigger.EditKeyPressed | QAbstractItemView.EditTrigger.SelectedClicked)
        self.sample_table.verticalHeader().setDefaultSectionSize(31)
        self.sample_table.horizontalHeader().setSectionsMovable(True)
        self.sample_table.setMinimumHeight(410)
        layout.addWidget(self.sample_table)
        self.metadata_requirements, self.metadata_requirements_model = _table(210)
        self.metadata_requirements_model.set_frame(metadata_requirement_table())
        _fit_table(self.metadata_requirements, 420)
        layout.addWidget(self.metadata_requirements)
        action = QHBoxLayout()
        self.sample_apply_status = QLabel("Review sample roles before applying.")
        apply_button = QPushButton("Apply sample map and rebuild evidence")
        apply_button.setObjectName("primaryButton")
        apply_button.clicked.connect(self._apply_sample_metadata)
        action.addWidget(self.sample_apply_status, 1)
        action.addWidget(apply_button)
        layout.addLayout(action)
        return page

    def _build_structure_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 03 / Structure and MS2 evidence",
            "Inspect the six-file join and every fragmentation spectrum",
            "The audit compares source keys with CM, records duplicate rows and SHA-256 hashes, and retains every MSP peak. Select a spectrum to inspect its peak list.",
        ))
        cards = QHBoxLayout()
        self.structure_cards = [MetricCard("CM features"), MetricCard("CI candidates"), MetricCard("MS2 spectra"), MetricCard("Fragment peaks")]
        for card in self.structure_cards:
            cards.addWidget(card)
        layout.addLayout(cards)
        self.audit_table, self.audit_model = _table(235)
        layout.addWidget(self.audit_table)
        self.structure_warning_layout = QVBoxLayout()
        layout.addLayout(self.structure_warning_layout)
        search = QHBoxLayout()
        search.addWidget(QLabel("Search spectra"))
        self.spectrum_search = QLineEdit()
        self.spectrum_search.setPlaceholderText("Feature ID, database ID, precursor type, name, or formula…")
        self.spectrum_search.textChanged.connect(self._filter_spectra)
        search.addWidget(self.spectrum_search, 1)
        spectrum_export = QPushButton("Export all MS2 peaks CSV")
        spectrum_export.setObjectName("secondaryButton")
        spectrum_export.clicked.connect(lambda: self._export_stage_csv("ms2_peaks"))
        search.addWidget(spectrum_export)
        layout.addLayout(search)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.spectrum_table, self.spectrum_model = _table(260)
        self.peak_table, self.peak_model = _table(220)
        splitter.addWidget(self.spectrum_table)
        splitter.addWidget(self.peak_table)
        splitter.setSizes([360, 270])
        layout.addWidget(splitter)
        self.spectrum_table.selectionModel().currentRowChanged.connect(self._spectrum_selected)
        return page

    def _build_identity_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 04 / Accepted identity audit",
            "Explore accepted compounds and repeated chromatographic reads",
            "Accepted Compound ID is the primary entity key. Description, Formula, Adducts, source links, chemistry, CI counts, and six-file evidence remain attached to every read.",
        ))
        cards = QHBoxLayout()
        self.identity_cards = [MetricCard("Accepted reads"), MetricCard("Unique entities"), MetricCard("Repeated-ID groups"), MetricCard("All CI candidates")]
        for card in self.identity_cards:
            cards.addWidget(card)
        layout.addLayout(cards)
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Search compounds"))
        self.identity_search = QLineEdit()
        self.identity_search.setPlaceholderText("Accepted ID, description, formula, adduct, family, or Feature ID…")
        self.identity_search.textChanged.connect(self._filter_identity_table)
        self.repeated_only = QCheckBox("Repeated identities only")
        self.repeated_only.toggled.connect(self._filter_identity_table)
        controls.addWidget(self.identity_search, 1)
        controls.addWidget(self.repeated_only)
        layout.addLayout(controls)
        self.identity_table, self.identity_model = _table(470)
        layout.addWidget(self.identity_table)
        buttons = QHBoxLayout()
        export_all = QPushButton("Export all accepted reads")
        export_all.setObjectName("secondaryButton")
        export_all.clicked.connect(lambda: self._export_stage_csv("accepted"))
        export_repeated = QPushButton("Export repeated reads")
        export_repeated.setObjectName("secondaryButton")
        export_repeated.clicked.connect(lambda: self._export_stage_csv("duplicates"))
        next_button = QPushButton("Review representative selection →")
        next_button.setObjectName("primaryButton")
        next_button.clicked.connect(lambda: self._go_to(4))
        buttons.addWidget(export_all)
        buttons.addWidget(export_repeated)
        buttons.addStretch(1)
        buttons.addWidget(next_button)
        layout.addLayout(buttons)
        return page

    def _build_selection_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 05 / Representative selection",
            "Choose one reproducible read per accepted identity",
            "Unique entities pass unchanged. Repeated reads use a deterministic lexicographic hierarchy, preserving every alternative and the complete comparison trace.",
        ))
        layout.addWidget(_callout(
            "<b>Primary hierarchy:</b> Score ↓ → Fragmentation Score ↓ → absolute Mass Error ↑ → Isotope Similarity ↓. "
            "<b>If all four tie:</b> evidence completeness ↓ → linked fragment peaks ↓ → pooled-QC detection ↓ → QC CV ↑ → biological detection/abundance ↓ → stable Feature ID. "
            "A tie across all available dimensions is flagged for manual review.",
        ))
        cards = QHBoxLayout()
        self.selection_cards = [MetricCard("Representatives"), MetricCard("Repeated-ID winners"), MetricCard("Reads collapsed"), MetricCard("Manual-review ties")]
        for card in self.selection_cards:
            cards.addWidget(card)
        layout.addLayout(cards)
        self.selection_table, self.selection_model = _table(470)
        layout.addWidget(self.selection_table)
        buttons = QHBoxLayout()
        export_selected = QPushButton("Export representatives")
        export_selected.setObjectName("secondaryButton")
        export_selected.clicked.connect(lambda: self._export_stage_csv("selected"))
        next_button = QPushButton("Review QC and thresholds →")
        next_button.setObjectName("primaryButton")
        next_button.clicked.connect(lambda: self._go_to(5))
        buttons.addWidget(export_selected)
        buttons.addStretch(1)
        buttons.addWidget(next_button)
        layout.addLayout(buttons)
        return page

    def _build_qc_filter_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 06 / Analytical QC and confidence filter",
            "Combine dataset-adaptive identity evidence with mapped-sample QC",
            "Recommendations are estimated only after identity collapse. Every rule remains editable; every entity receives individual pass/fail columns and reasons.",
        ))
        preset_row = QHBoxLayout()
        self.preset_buttons: dict[str, QPushButton] = {}
        for name in ("Inclusive", "Balanced", "Stringent"):
            button = QPushButton(name)
            button.setObjectName("secondaryButton")
            button.clicked.connect(lambda _checked=False, selected=name: self._set_preset(selected))
            preset_row.addWidget(button)
            self.preset_buttons[name] = button
        layout.addLayout(preset_row)
        split = QSplitter(Qt.Orientation.Horizontal)
        identity_group = QGroupBox("Identification-evidence thresholds")
        identity_grid = QGridLayout(identity_group)
        headers = ("Metric", "Manual cutoff", "Bootstrap 95% CI", "Method and rationale")
        for column, header in enumerate(headers):
            identity_grid.addWidget(QLabel(header.upper()), 0, column)
        self.threshold_spins: dict[str, QDoubleSpinBox] = {}
        self.threshold_meta: dict[str, tuple[QLabel, QLabel]] = {}
        specs = (
            ("score", "Score ≥", 0, 100, 1),
            ("fragmentation", "Fragmentation Score ≥", 0, 100, 1),
            ("abs_mass_error", "Absolute mass error ≤ ppm", 0, 100, 1),
            ("isotope", "Isotope Similarity ≥", 0, 100, 1),
        )
        for row, (key, label, low, high, decimals) in enumerate(specs, 1):
            identity_grid.addWidget(QLabel(label), row, 0)
            spin = QDoubleSpinBox()
            spin.setRange(low, high)
            spin.setDecimals(decimals)
            spin.valueChanged.connect(self._update_expected_passes)
            ci = QLabel("—")
            rationale = QLabel("—")
            rationale.setWordWrap(True)
            identity_grid.addWidget(spin, row, 1)
            identity_grid.addWidget(ci, row, 2)
            identity_grid.addWidget(rationale, row, 3)
            self.threshold_spins[key] = spin
            self.threshold_meta[key] = (ci, rationale)
        self.allow_zero = QCheckBox("Permit Fragmentation Score = 0 as an explicit unavailable/disabled/unmatched state")
        self.allow_zero.toggled.connect(self._update_expected_passes)
        identity_grid.addWidget(self.allow_zero, 5, 0, 1, 4)
        identity_grid.setColumnStretch(3, 1)

        qc_group = QGroupBox("Analytical-QC rules")
        qc_grid = QGridLayout(qc_group)
        self.apply_qc_rules = QCheckBox("Apply pooled-QC detection and CV")
        self.apply_qc_rules.setChecked(True)
        self.apply_blank_rule = QCheckBox("Apply sample/blank ratio when blanks exist")
        self.apply_blank_rule.setChecked(True)
        self.apply_drift_rule = QCheckBox("Apply QC drift limit")
        self.apply_dratio_rule = QCheckBox("Apply D-ratio limit")
        for checkbox in (self.apply_qc_rules, self.apply_blank_rule, self.apply_drift_rule, self.apply_dratio_rule):
            checkbox.toggled.connect(self._update_expected_passes)
        qc_grid.addWidget(self.apply_qc_rules, 0, 0, 1, 2)
        qc_grid.addWidget(self.apply_blank_rule, 1, 0, 1, 2)
        qc_grid.addWidget(self.apply_drift_rule, 2, 0, 1, 2)
        qc_grid.addWidget(self.apply_dratio_rule, 3, 0, 1, 2)
        self.qc_spins: dict[str, QDoubleSpinBox] = {}
        qc_specs = (
            ("min_qc_detection", "Minimum QC detection (%)", 0, 100, 80),
            ("max_qc_cv", "Maximum QC CV (%)", 0, 300, 30),
            ("min_biological_detection", "Minimum biological detection (%)", 0, 100, 20),
            ("min_biological_blank_ratio", "Minimum sample/blank median ratio", 0, 1000, 3),
            ("max_abs_qc_drift", "Maximum absolute QC drift (%)", 0, 1000, 40),
            ("max_d_ratio", "Maximum D-ratio (%)", 0, 100, 50),
        )
        for row, (key, label, low, high, value) in enumerate(qc_specs, 4):
            qc_grid.addWidget(QLabel(label), row, 0)
            spin = QDoubleSpinBox()
            spin.setRange(low, high)
            spin.setDecimals(1)
            spin.setValue(value)
            spin.valueChanged.connect(self._update_expected_passes)
            qc_grid.addWidget(spin, row, 1)
            self.qc_spins[key] = spin
        self.require_ms2 = QCheckBox("Require linked MS2 peaks")
        self.require_ms2.toggled.connect(self._update_expected_passes)
        self.min_fragment_peaks = QDoubleSpinBox()
        self.min_fragment_peaks.setRange(1, 500)
        self.min_fragment_peaks.setDecimals(0)
        self.min_fragment_peaks.setValue(3)
        self.min_fragment_peaks.valueChanged.connect(self._update_expected_passes)
        qc_grid.addWidget(self.require_ms2, 10, 0)
        qc_grid.addWidget(self.min_fragment_peaks, 10, 1)
        split.addWidget(identity_group)
        split.addWidget(qc_group)
        split.setSizes([850, 520])
        layout.addWidget(split)
        self.threshold_logic = QLabel("Recommendations appear after scanning.")
        self.threshold_logic.setObjectName("logicBox")
        self.threshold_logic.setWordWrap(True)
        self.threshold_logic.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(self.threshold_logic)
        self.profile_table, self.profile_model = _table(190)
        layout.addWidget(self.profile_table)
        self.qc_preview_table, self.qc_preview_model = _table(310)
        layout.addWidget(self.qc_preview_table)
        action = QHBoxLayout()
        self.expected_passes = QLabel("Scan files to estimate a result.")
        self.expected_passes.setObjectName("resultLabel")
        apply_button = QPushButton("Apply all filters and continue →")
        apply_button.setObjectName("primaryButton")
        apply_button.clicked.connect(self._apply_thresholds)
        action.addWidget(self.expected_passes, 1)
        action.addWidget(apply_button)
        layout.addLayout(action)
        layout.addWidget(_callout(
            "No universal Score/CV/blank threshold is inferred as a confirmed boundary. Presets are transparent triage settings. Without target-decoy results or authentic standards, do not report the shortlist as an identification FDR.",
            "warning",
        ))
        return page

    def _build_context_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 07 / Chemistry and contextual interpretation",
            "Score source evidence and extraction-location plausibility",
            "The model combines explicit description/ontology cues with declared assay context. Formula properties and evidence provenance are exported so low-confidence, prior-driven results can be reviewed.",
        ))
        context = QGroupBox("Context used for scoring")
        grid = QGridLayout(context)
        grid.addWidget(QLabel("Assay type"), 0, 0)
        self.context_assay = QComboBox()
        self.context_assay.addItem("Metabolomics", "metabolomics")
        self.context_assay.addItem("Lipidomics", "lipidomics")
        grid.addWidget(self.context_assay, 0, 1)
        grid.addWidget(QLabel("Biological system"), 0, 2)
        self.context_system = QComboBox()
        _add_combo_items(self.context_system, BIOLOGICAL_SYSTEMS)
        grid.addWidget(self.context_system, 0, 3)
        grid.addWidget(QLabel("Sample matrix"), 1, 0)
        self.context_matrix = QComboBox()
        self.context_matrix.setEditable(True)
        self.context_matrix.addItems(SAMPLE_MATRICES)
        grid.addWidget(self.context_matrix, 1, 1)
        grid.addWidget(QLabel("Exposure context"), 1, 2)
        self.exposure_context = QLineEdit()
        self.exposure_context.setPlaceholderText("e.g. none, antibiotic treatment, virus exposure, pesticide exposure")
        grid.addWidget(self.exposure_context, 1, 3)
        grid.addWidget(QLabel("Extraction method"), 2, 0)
        self.context_extraction = QComboBox()
        _add_combo_items(self.context_extraction, EXTRACTION_METHODS)
        grid.addWidget(self.context_extraction, 2, 1)
        grid.addWidget(QLabel("Analyzed phase"), 2, 2)
        self.context_phase = QComboBox()
        _add_combo_items(self.context_phase, ANALYZED_PHASES)
        grid.addWidget(self.context_phase, 2, 3)
        score_button = QPushButton("Score source and extraction-location evidence")
        score_button.setObjectName("primaryButton")
        score_button.clicked.connect(self._score_source)
        grid.addWidget(score_button, 3, 0, 1, 4)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        layout.addWidget(context)
        layout.addWidget(_callout(
            "<b>Source model:</b> six competing evidence scores plus entropy, margin, direct-evidence count, and rule provenance. "
            "<b>Phase model:</b> explicit behavior labels, parsed MolLogP/TPSA when present, formula-polarity proxy, lipid cues, and method-specific LLE/SPE logic. "
            "Outputs are review likelihoods—not causal origin or measured recovery.",
        ))
        cards = QHBoxLayout()
        self.source_cards = [MetricCard("Shortlisted"), MetricCard("Human endogenous"), MetricCard("Microbiome-derived"), MetricCard("Co-metabolic / other")]
        for card in self.source_cards:
            cards.addWidget(card)
        layout.addLayout(cards)
        self.source_table, self.source_model = _table(480)
        layout.addWidget(self.source_table)
        action = QHBoxLayout()
        export_short = QPushButton("Export shortlist")
        export_short.setObjectName("secondaryButton")
        export_short.clicked.connect(lambda: self._export_stage_csv("shortlist"))
        export_source = QPushButton("Export chemistry + source + phase")
        export_source.setObjectName("secondaryButton")
        export_source.clicked.connect(lambda: self._export_stage_csv("provenance"))
        next_button = QPushButton("Review exports →")
        next_button.setObjectName("primaryButton")
        next_button.clicked.connect(lambda: self._go_to(7))
        action.addWidget(export_short)
        action.addWidget(export_source)
        action.addStretch(1)
        action.addWidget(next_button)
        layout.addLayout(action)
        return page

    def _build_export_page(self) -> QWidget:
        page, layout = self._page()
        layout.addWidget(_section_heading(
            "Step 08 / Journal-ready reporting and reproducible export",
            "Review the dataset-specific methodology, analysis, and complete evidence package",
            "The manuscript text is generated from the selected options, mapped samples, audited pipeline state, active rules, and final shortlist. The export retains the report, task history, numbered CSVs, formatted XLSX, JSON manifest, and scientific limitations.",
        ))
        report_group = QGroupBox("Dataset-specific methodology and curation analysis")
        report_layout = QVBoxLayout(report_group)
        self.journal_status = QLabel("Apply filters and open this stage to generate the report.")
        self.journal_status.setObjectName("bodyText")
        self.journal_status.setWordWrap(True)
        report_toolbar = QHBoxLayout()
        refresh_report = QPushButton("Generate / refresh report")
        refresh_report.setObjectName("primaryButton")
        refresh_report.clicked.connect(self._generate_journal_report)
        copy_report = QPushButton("Copy Markdown")
        copy_report.setObjectName("secondaryButton")
        copy_report.clicked.connect(self._copy_journal_report)
        save_report = QPushButton("Save report .md")
        save_report.setObjectName("secondaryButton")
        save_report.clicked.connect(self._save_journal_report)
        report_toolbar.addWidget(refresh_report)
        report_toolbar.addWidget(copy_report)
        report_toolbar.addWidget(save_report)
        report_toolbar.addStretch(1)
        report_toolbar.addWidget(self.journal_status)
        report_layout.addLayout(report_toolbar)
        self.journal_view = QTextBrowser()
        self.journal_view.setObjectName("journalView")
        self.journal_view.setOpenExternalLinks(True)
        self.journal_view.setMinimumHeight(520)
        self.journal_view.setMarkdown(
            "# Manuscript-ready LC–MS curation methodology and analysis\n\n"
            "The report will appear here after filtering. Contextual source and extraction scoring can be run first so those results are included."
        )
        report_layout.addWidget(self.journal_view)
        layout.addWidget(report_group)

        history_group = QGroupBox("Background task and notification audit")
        history_layout = QVBoxLayout(history_group)
        history_layout.addWidget(QLabel("Every background computation records its stage, completion state, timestamps, duration, and outcome. This table is included in complete exports."))
        self.task_history_table, self.task_history_model = _table(220)
        history_layout.addWidget(self.task_history_table)
        layout.addWidget(history_group)

        self.export_table, self.export_model = _table(380)
        layout.addWidget(self.export_table)
        options = QFrame()
        options.setObjectName("panel")
        option_layout = QHBoxLayout(options)
        self.include_candidates = QCheckBox("Include every CI candidate in XLSX (large workbook)")
        option_layout.addWidget(self.include_candidates)
        option_layout.addStretch(1)
        workbook_button = QPushButton("Complete XLSX workbook")
        workbook_button.setObjectName("secondaryButton")
        workbook_button.clicked.connect(self._export_workbook)
        folder_button = QPushButton("Export complete results folder")
        folder_button.setObjectName("primaryButton")
        folder_button.clicked.connect(self._export_complete)
        option_layout.addWidget(workbook_button)
        option_layout.addWidget(folder_button)
        layout.addWidget(options)
        layout.addWidget(_callout(
            "The generated text is a dataset-specific draft, not a substitute for the experimental LC–MS method. Complete all checklist items that cannot be recovered from the six exports. "
            "Keep the manifest, sample metadata, task history, file audit, repeated reads, filter decisions, and report with downstream analysis.",
        ))
        return page

    def _navigate(self, index: int) -> None:
        if index < 0:
            return
        if index > self._max_stage:
            self.stage_list.blockSignals(True)
            self.stage_list.setCurrentRow(self.stack.currentIndex())
            self.stage_list.blockSignals(False)
            return
        self.stack.setCurrentIndex(index)

    def _unlock_stage(self, stage: int) -> None:
        self._max_stage = max(self._max_stage, stage)
        for index in range(self.stage_list.count()):
            item = self.stage_list.item(index)
            if index <= self._max_stage:
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEnabled)

    def _go_to(self, stage: int) -> None:
        self._unlock_stage(stage)
        self.stage_list.setCurrentRow(stage)
        if stage == 7 and self.filter_settings and not self.journal_report and self._active_task is None:
            QTimer.singleShot(0, self._generate_journal_report)

    def _browse_role(self, role: str) -> None:
        filters = {"II": "XML (*.xml)", "FD": "MSP (*.msp)", "CM": "CSV (*.csv)", "AM": "CSV (*.csv)", "CI": "CSV (*.csv)", "ACP": "CSV (*.csv)"}
        current = self.path_edits[role].text() or str(Path.home())
        path, _ = QFileDialog.getOpenFileName(self, f"Select {role}", current, filters[role] + ";;All files (*.*)")
        if path:
            self.path_edits[role].setText(path)
            self.settings.setValue(f"path/{role}", path)

    def _select_all_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(self, "Select CM, AM, CI, II, FD, and ACP files", str(Path.home()), "LC–MS exports (*.csv *.xml *.msp);;All files (*.*)")
        if not paths:
            return
        unresolved = []
        for path in paths:
            stem = Path(path).stem.upper()
            role = next((candidate for candidate in FILE_ROLES if stem == candidate or stem.startswith(candidate + "_")), None)
            if role:
                self.path_edits[role].setText(path)
                self.settings.setValue(f"path/{role}", path)
            else:
                unresolved.append(Path(path).name)
        if unresolved:
            QMessageBox.information(self, "Unassigned files", "These filenames could not be assigned automatically:\n" + "\n".join(unresolved) + "\n\nUse each Browse button to assign them.")

    def _clear_paths(self) -> None:
        for role, edit in self.path_edits.items():
            edit.clear()
            self.settings.remove(f"path/{role}")

    def _collect_parameters(self) -> dict[str, Any]:
        parameters = {
            "project_name": self.project_name.text().strip() or "LCMS Curation Project",
            "analyst": self.analyst.text().strip(),
            "assay_type": self.assay_type.currentData(),
            "study_type": self.assay_type.currentData(),
            "biological_system": self.biological_system.currentData(),
            "sample_matrix": self.sample_matrix.currentText().strip(),
            "extraction_method": self.extraction_method.currentData(),
            "analyzed_phase": self.analyzed_phase.currentData(),
            "qc_pooling_point": self.qc_pooling.currentData(),
            "ion_mode": self.ion_mode.currentData(),
            "acquisition_mode": self.acquisition_mode.currentData(),
            "mass_tolerance": self.mass_tolerance.value(),
            "qc_pattern": self.qc_pattern.text().strip() or "QC|Pool",
        }
        if self.pending_sample_metadata is not None:
            parameters["sample_metadata"] = self.pending_sample_metadata
        return parameters

    def _start_analysis(self) -> None:
        files = {role: self.path_edits[role].text().strip() for role in FILE_ROLES}
        missing = [role for role, path in files.items() if not path or not Path(path).is_file()]
        if missing:
            QMessageBox.warning(self, "Files required", "Select a valid file for: " + ", ".join(missing))
            return
        if not self._begin_task("Step 01 · Project + files", "Scan, validate, and reconcile all six LC–MS exports"):
            return
        parameters = self._collect_parameters()
        self.settings.setValue("project/name", parameters["project_name"])
        self.settings.setValue("project/analyst", parameters["analyst"])
        self.settings.setValue("parameters/mass_tolerance", parameters["mass_tolerance"])
        self.settings.setValue("parameters/qc_pattern", parameters["qc_pattern"])
        self.scan_button.setEnabled(False)
        self.scan_progress.setValue(0)
        self.scan_status.setText("Starting…")
        worker = AnalysisWorker(files, parameters, self)
        self._workers.append(worker)
        worker.progressed.connect(self._analysis_progress)
        worker.succeeded.connect(self._analysis_ready)
        worker.failed.connect(self._analysis_failed)
        worker.finished.connect(lambda: self._forget_worker(worker))
        worker.start()

    def _analysis_progress(self, percent: int, message: str) -> None:
        self.scan_progress.setValue(percent)
        self.scan_status.setText(message)
        self._task_progress(percent, message)

    def _analysis_ready(self, result: AnalysisBundle) -> None:
        self.analysis = result
        self.pending_sample_metadata = None
        self.shortlist = pd.DataFrame()
        self.filter_decisions = pd.DataFrame()
        self.provenance = pd.DataFrame()
        self.filter_settings = {}
        self.provenance_settings = None
        self.journal_report = ""
        self.scan_button.setEnabled(True)
        self.scan_progress.setValue(100)
        self.scan_status.setText("Complete")
        self._populate_sample_page()
        self._populate_structure_page()
        self._populate_identity_page()
        self._populate_selection_page()
        self._populate_threshold_page()
        self._populate_export_page()
        self._sync_context_controls()
        self._unlock_stage(5)
        self._go_to(1)
        self._complete_task(
            f"Analysis ready: {len(result.selected_reads):,} unique entities; {result.duplicate_groups:,} repeated IDs; {len(result.fd_peaks):,} MS² peaks."
        )

    def _analysis_failed(self, trace: str) -> None:
        self.scan_button.setEnabled(True)
        self.scan_status.setText("Failed")
        final_line = trace.strip().splitlines()[-1] if trace.strip() else "Unknown error"
        self._fail_task(final_line)
        details = QMessageBox(self)
        details.setIcon(QMessageBox.Icon.Critical)
        details.setWindowTitle("Analysis failed")
        details.setText(final_line)
        details.setDetailedText(trace)
        details.exec()

    def _populate_sample_page(self) -> None:
        assert self.analysis is not None
        frame = self.analysis.sample_metadata.copy()
        self.sample_model.set_frame(frame)
        role_column = frame.columns.get_loc("Sample role")
        self.sample_table.setItemDelegateForColumn(role_column, ComboBoxDelegate(SAMPLE_ROLES, self.sample_table))
        _fit_table(self.sample_table, 220)
        for hidden in ("Normalized abundance column", "Raw abundance column"):
            if hidden in frame.columns:
                self.sample_table.setColumnHidden(frame.columns.get_loc(hidden), True)
        counts = frame["Sample role"].value_counts()
        blanks = sum(int(counts.get(role, 0)) for role in ("Process blank", "Extraction blank", "Solvent blank"))
        self.sample_cards[0].set_metric(f"{len(frame):,}", "normalised sample columns")
        self.sample_cards[1].set_metric(f"{int(counts.get('Biological', 0)):,}", "study samples")
        self.sample_cards[2].set_metric(f"{int(counts.get('Pooled QC', 0)):,}", "repeatability/drift evidence")
        self.sample_cards[3].set_metric(f"{blanks:,}", "contamination evidence" if blanks else "not evaluable")

    def _import_sample_metadata(self) -> None:
        if self.analysis is None:
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import sample metadata", str(Path.home()), "CSV (*.csv);;All files (*.*)")
        if not path:
            return
        try:
            imported = pd.read_csv(path, dtype=str, keep_default_na=False)
            merged = merge_sample_metadata(self.sample_model.frame(), imported)
            self.sample_model.set_frame(merged)
            self.sample_apply_status.setText(f"Imported metadata for {len(imported):,} rows; click Apply to recompute evidence.")
            self._record_instant_task("Step 02 · Sample map", "Import sample metadata", f"Imported and reconciled {len(imported):,} metadata rows.")
        except Exception as exc:
            QMessageBox.critical(self, "Metadata import failed", str(exc))

    def _reset_sample_roles(self) -> None:
        if self.analysis is None:
            return
        frame = build_sample_metadata(
            self.analysis.normalized_columns,
            self.analysis.raw_columns,
            self.qc_pattern.text().strip() or "QC|Pool",
            self._collect_parameters(),
        )
        self.sample_model.set_frame(frame)
        self.sample_apply_status.setText("Automatic role rules re-applied; review and click Apply.")
        self._record_instant_task("Step 02 · Sample map", "Re-run automatic sample-role rules", f"Automatic roles regenerated for {len(frame):,} sample columns.")

    def _apply_sample_metadata(self) -> None:
        if self.analysis is None:
            return
        analysis = self.analysis
        metadata = self.sample_model.frame()
        self._run_background_task(
            "Step 02 · Sample map",
            "Recompute sample-role, abundance, QC, representative, and threshold evidence",
            lambda progress: rebuild_with_sample_metadata(analysis, metadata, progress),
            on_success=self._sample_metadata_ready,
            outcome=lambda result: (
                f"Sample map applied to {len(result.sample_metadata):,} samples; "
                f"{len(result.selected_reads):,} representatives and adaptive thresholds rebuilt."
            ),
            failure_title="Sample map failed",
        )

    def _sample_metadata_ready(self, result: AnalysisBundle) -> None:
        self.analysis = result
        self.shortlist = pd.DataFrame()
        self.filter_decisions = pd.DataFrame()
        self.provenance = pd.DataFrame()
        self.filter_settings = {}
        self.provenance_settings = None
        self.journal_report = ""
        self._populate_sample_page()
        self._populate_structure_page()
        self._populate_identity_page()
        self._populate_selection_page()
        self._populate_threshold_page()
        self._populate_export_page()
        self.sample_apply_status.setText("Sample map applied; representative tie-breaks, QC metrics, and thresholds rebuilt.")
        self._go_to(2)

    def _populate_structure_page(self) -> None:
        assert self.analysis is not None
        self.structure_cards[0].set_metric(f"{len(self.analysis.cm):,}", "CM feature rows")
        self.structure_cards[1].set_metric(f"{len(self.analysis.ci):,}", "accepted + rejected")
        self.structure_cards[2].set_metric(f"{len(self.analysis.fd):,}", "MSP records")
        self.structure_cards[3].set_metric(f"{len(self.analysis.fd_peaks):,}", "peak-level rows")
        self.audit_model.set_frame(self.analysis.audits)
        _fit_table(self.audit_table, 300)
        while self.structure_warning_layout.count():
            item = self.structure_warning_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        for warning in self.analysis.warnings:
            self.structure_warning_layout.addWidget(_callout(html.escape(warning), "warning"))
        self._filter_spectra()

    def _filter_spectra(self) -> None:
        if self.analysis is None:
            return
        frame = self.analysis.fd
        query = self.spectrum_search.text().strip()
        if query:
            columns = [column for column in ("Feature ID", "FD database ID", "FD precursor type", "FD name", "FD formula") if column in frame.columns]
            mask = pd.Series(False, index=frame.index)
            for column in columns:
                mask |= frame[column].astype(str).str.contains(query, case=False, regex=False, na=False)
            frame = frame.loc[mask]
        columns = [
            "FD spectrum index", "Feature ID", "FD database ID", "FD precursor type", "FD precursor m/z", "FD formula",
            "FD observed peaks", "FD peak-count integrity", "FD base peak m/z", "FD spectral entropy",
            "FD normalized spectral entropy", "FD effective peak count", "FD top-five intensity share (%)", "FD spectrum quality flag",
        ]
        self.spectrum_model.set_frame(frame[[column for column in columns if column in frame.columns]])
        _fit_table(self.spectrum_table, 210)
        self.peak_model.set_frame(pd.DataFrame())

    def _spectrum_selected(self, current: QModelIndex, _previous: QModelIndex) -> None:
        if self.analysis is None or not current.isValid() or self.spectrum_model.frame().empty:
            return
        frame = self.spectrum_model.frame()
        if current.row() >= len(frame):
            return
        spectrum_index = frame.iloc[current.row()]["FD spectrum index"]
        peaks = self.analysis.fd_peaks.loc[self.analysis.fd_peaks["FD spectrum index"].eq(spectrum_index)]
        columns = ["Peak index", "Fragment m/z", "Fragment intensity", "Relative intensity (%)", "Neutral loss from precursor (Da)", "Base peak"]
        self.peak_model.set_frame(peaks[columns])
        _fit_table(self.peak_table, 260)

    def _populate_identity_page(self) -> None:
        assert self.analysis is not None
        analysis = self.analysis
        self.identity_cards[0].set_metric(f"{len(analysis.accepted_reads):,}", "accepted feature reads")
        self.identity_cards[1].set_metric(f"{len(analysis.selected_reads):,}", "accepted entities")
        self.identity_cards[2].set_metric(f"{analysis.duplicate_groups:,}", f"{len(analysis.duplicate_reads):,} reads preserved")
        self.identity_cards[3].set_metric(f"{len(analysis.ci):,}", "all candidates scanned")
        self._filter_identity_table()

    def _filter_identity_table(self) -> None:
        if self.analysis is None:
            return
        frame = self.analysis.accepted_reads
        query = self.identity_search.text().strip()
        if self.repeated_only.isChecked():
            frame = frame.loc[frame["Read count for entity"].gt(1)]
        if query:
            columns = [column for column in ("Feature ID", "Accepted Compound ID", "Accepted Description", "Formula", "Adducts", "Chemical family (text rule)") if column in frame.columns]
            mask = pd.Series(False, index=frame.index)
            for column in columns:
                mask |= frame[column].astype(str).str.contains(query, case=False, regex=False, na=False)
            frame = frame.loc[mask]
        columns = [
            "Feature ID", "Accepted Compound ID", "Accepted Description", "Formula", "Adducts", "Chemical family (text rule)",
            "Score", "Fragmentation Score", "Absolute Mass Error (ppm)", "Isotope Similarity", "Read count for entity",
            "FD spectrum count", "FD observed peaks", "QC CV%", "Biological detection rate", "Ion-mode/adduct consistency",
        ]
        self.identity_model.set_frame(frame[[column for column in columns if column in frame.columns]])
        _fit_table(self.identity_table, 230)
        self.statusBar().showMessage(f"Identity explorer: {len(frame):,} matching accepted reads.")

    def _populate_selection_page(self) -> None:
        assert self.analysis is not None
        selected = self.analysis.selected_reads
        winners = selected.loc[selected["Read count for entity"].gt(1)]
        self.selection_cards[0].set_metric(f"{len(selected):,}", "one per accepted entity")
        self.selection_cards[1].set_metric(f"{len(winners):,}", "one winner per repeated ID")
        self.selection_cards[2].set_metric(f"{self.analysis.extra_reads_collapsed:,}", "alternatives retained")
        self.selection_cards[3].set_metric(f"{int(selected['Manual identity review'].sum()):,}", "all available tie-breaks equal")
        columns = [
            "Accepted Compound ID", "Accepted Description", "Formula", "Adducts", "Feature ID", "Read count for entity",
            "Score", "Fragmentation Score", "Absolute Mass Error (ppm)", "Isotope Similarity", "Primary-metric tie",
            "Evidence completeness", "FD observed peaks", "QC detection rate", "QC CV%", "Selection reason", "Manual identity review", "Selection trace",
        ]
        self.selection_model.set_frame(winners[[column for column in columns if column in winners.columns]])
        _fit_table(self.selection_table, 250)

    def _populate_threshold_page(self) -> None:
        assert self.analysis is not None
        recommendation = self.analysis.thresholds
        for key in ("score", "fragmentation", "abs_mass_error", "isotope"):
            field = recommendation[key]
            ci_label, method_label = self.threshold_meta[key]
            ci_label.setText(f"{field['ci_low']:g}–{field['ci_high']:g}")
            method_label.setText(f"{field['method']}\n{field['rationale']}")
        diagnostics = recommendation["diagnostics"]
        mass_text = (
            f"Stable mass-error boundary <b>{diagnostics['mass_intersection']:.2f} ppm</b>; ΔBIC <b>{diagnostics['mass_delta_bic']:.1f}</b>, Ashman D <b>{diagnostics['mass_ashman_d']:.2f}</b>."
            if diagnostics["mass_mixture_used"]
            else "No defensible two-cluster mass-error split was accepted; a conservative robust bound was used."
        )
        qc_count = int((self.analysis.sample_metadata["Sample role"] == "Pooled QC").sum())
        blank_count = int(self.analysis.sample_metadata["Sample role"].isin({"Process blank", "Extraction blank", "Solvent blank"}).sum())
        self.threshold_logic.setText(
            "<b>Recommendation logic.</b> Score and isotope: selected-entity lower quartile. Fragmentation: lower quartile among positive scores, with zero modeled separately. "
            + mass_text
            + " Five hundred deterministic entity-level bootstrap resamples quantify cutoff stability. "
            + f"Mapped evidence: <b>{qc_count} pooled-QC</b> and <b>{blank_count} blank</b> samples. Missing evidence is marked not evaluable, not silently passed as observed quality."
        )
        profile_rows = [{"Metric": metric, **values} for metric, values in recommendation["profiles"].items()]
        self.profile_model.set_frame(pd.DataFrame(profile_rows)[["Metric", "n", "min", "q10", "q25", "median", "q75", "q90", "max"]])
        _fit_table(self.profile_table, 220)
        for preset in recommendation["presets"]:
            button = self.preset_buttons[preset["name"]]
            button.setText(f"{preset['name']}\n{preset['estimated_passes']:,} identity-only passes")
            button.setToolTip(preset["description"])
        preview_columns = [
            "Accepted Compound ID", "Accepted Description", "Score", "Fragmentation Score", "Absolute Mass Error (ppm)", "Isotope Similarity",
            "QC detection rate", "QC CV%", "QC robust CV% (MAD)", "QC drift across run (%)", "D-ratio (%)",
            "Biological detection rate", "Biological/blank median ratio", "Analytical QC decision", "Analytical QC fail reasons",
        ]
        self.qc_preview_model.set_frame(self.analysis.analytical_qc[[column for column in preview_columns if column in self.analysis.analytical_qc.columns]])
        _fit_table(self.qc_preview_table, 230)
        self._set_preset("Balanced")

    def _set_preset(self, name: str) -> None:
        if self.analysis is None:
            return
        preset = next(item for item in self.analysis.thresholds["presets"] if item["name"] == name)
        for key in ("score", "fragmentation", "abs_mass_error", "isotope"):
            self.threshold_spins[key].blockSignals(True)
            self.threshold_spins[key].setValue(float(preset[key]))
            self.threshold_spins[key].blockSignals(False)
        self.allow_zero.blockSignals(True)
        self.allow_zero.setChecked(bool(preset["allow_zero_fragmentation"]))
        self.allow_zero.blockSignals(False)
        self._update_expected_passes()

    def _current_filters(self) -> dict[str, Any]:
        settings = {
            "score": self.threshold_spins["score"].value(),
            "fragmentation": self.threshold_spins["fragmentation"].value(),
            "abs_mass_error": self.threshold_spins["abs_mass_error"].value(),
            "isotope": self.threshold_spins["isotope"].value(),
            "allow_zero_fragmentation": self.allow_zero.isChecked(),
            **default_qc_settings(),
        }
        settings.update(
            {
                "min_qc_detection": self.qc_spins["min_qc_detection"].value() / 100,
                "max_qc_cv": self.qc_spins["max_qc_cv"].value(),
                "min_biological_detection": self.qc_spins["min_biological_detection"].value() / 100,
                "min_biological_blank_ratio": self.qc_spins["min_biological_blank_ratio"].value(),
                "max_abs_qc_drift": self.qc_spins["max_abs_qc_drift"].value(),
                "max_d_ratio": self.qc_spins["max_d_ratio"].value(),
                "apply_qc_filter": self.apply_qc_rules.isChecked(),
                "apply_blank_filter": self.apply_blank_rule.isChecked(),
                "apply_drift_filter": self.apply_drift_rule.isChecked(),
                "apply_d_ratio_filter": self.apply_dratio_rule.isChecked(),
                "require_ms2": self.require_ms2.isChecked(),
                "min_fragment_peaks": int(self.min_fragment_peaks.value()),
            }
        )
        return settings

    def _update_expected_passes(self) -> None:
        if self.analysis is None:
            return
        decisions = evaluate_filters(self.analysis.selected_reads, self._current_filters())
        count = int(decisions["Filter pass"].sum())
        share = count / len(decisions) * 100 if len(decisions) else 0
        id_pass = int(decisions["Identification evidence pass"].sum())
        qc_fail = int((~decisions["Analytical QC pass"]).sum())
        self.expected_passes.setText(
            f"Expected shortlist: {count:,} / {len(decisions):,} ({share:.1f}%). Identity-only pass: {id_pass:,}; analytical-QC failures: {qc_fail:,}."
        )

    def _apply_thresholds(self) -> None:
        if self.analysis is None:
            return
        settings = self._current_filters()
        selected = self.analysis.selected_reads.copy()

        def compute(progress: Callable[[int, str], None]) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
            progress(8, "Evaluating identification-evidence thresholds for every representative…")
            decisions = evaluate_filters(selected, settings)
            progress(52, "Applying enabled and evaluable analytical-QC rules…")
            shortlist = apply_filters(selected, settings)
            progress(86, "Assembling rule-level decisions and shortlisted compound details…")
            progress(100, "Identity and analytical-QC filtering complete")
            return settings, decisions, shortlist

        self._run_background_task(
            "Step 06 · QC + filter",
            "Apply adaptive identity thresholds and analytical-QC rules",
            compute,
            on_success=self._thresholds_ready,
            outcome=lambda result: f"Filters applied: {len(result[2]):,} of {len(selected):,} representatives shortlisted.",
            failure_title="Filtering failed",
        )

    def _thresholds_ready(self, result: tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]) -> None:
        self.filter_settings, self.filter_decisions, self.shortlist = result
        self.provenance = pd.DataFrame()
        self.provenance_settings = None
        self.journal_report = ""
        self.source_cards[0].set_metric(f"{len(self.shortlist):,}", "passed identity + enabled QC rules")
        for card in self.source_cards[1:]:
            card.set_metric("—", "run contextual scoring")
        preview_columns = [
            "Accepted Compound ID", "Accepted Description", "Formula", "Chemical family (text rule)", "Formula DBE",
            "Formula polarity proxy (0-100)", "Formula polarity class", "Adducts", "Ion-mode/adduct consistency",
            "Score", "Fragmentation Score", "Absolute Mass Error (ppm)", "Isotope Similarity", "QC CV%", "Filter fail reasons",
        ]
        self.source_model.set_frame(self.shortlist[[column for column in preview_columns if column in self.shortlist.columns]])
        _fit_table(self.source_table, 230)
        self._populate_export_page()
        self._unlock_stage(6)
        self._go_to(6)

    def _sync_context_controls(self) -> None:
        _set_combo_data(self.context_assay, self.assay_type.currentData())
        _set_combo_data(self.context_system, self.biological_system.currentData())
        self.context_matrix.setCurrentText(self.sample_matrix.currentText())
        _set_combo_data(self.context_extraction, self.extraction_method.currentData())
        _set_combo_data(self.context_phase, self.analyzed_phase.currentData())

    def _score_source(self) -> None:
        if self.shortlist.empty:
            QMessageBox.warning(self, "Apply filters first", "Apply the identity and analytical-QC filters before contextual scoring.")
            return
        matrix = self.context_matrix.currentText().strip()
        if not matrix:
            QMessageBox.warning(self, "Sample matrix required", "Enter the sample matrix before scoring source evidence.")
            return
        settings = {
            "study_type": self.context_assay.currentData(),
            "assay_type": self.context_assay.currentData(),
            "biological_system": self.context_system.currentData(),
            "sample_matrix": matrix,
            "exposure_context": self.exposure_context.text().strip() or "none",
            "extraction_method": self.context_extraction.currentData(),
            "analyzed_phase": self.context_phase.currentData(),
        }
        shortlist = self.shortlist.copy()
        self._run_background_task(
            "Step 07 · Context model",
            "Score chemistry, competing source classes, and extraction-location plausibility",
            lambda progress: (settings, score_provenance(shortlist, settings, progress)),
            on_success=self._source_scoring_ready,
            outcome=lambda result: f"Contextual source and extraction-location evidence scored for {len(result[1]):,} compounds.",
            failure_title="Contextual scoring failed",
        )

    def _source_scoring_ready(self, result: tuple[dict[str, Any], pd.DataFrame]) -> None:
        self.provenance_settings, self.provenance = result
        self.journal_report = ""
        if self.analysis is not None:
            self.analysis.parameters.update(self.provenance_settings)
        counts = self.provenance["Primary source class"].value_counts()
        human = int(counts.get("Human endogenous", 0))
        microbe = int(counts.get("Microbiome-derived", 0))
        co = int(counts.get("Host–microbe co-metabolic", 0))
        other = len(self.provenance) - human - microbe - co
        self.source_cards[0].set_metric(f"{len(self.provenance):,}", "contextually scored")
        self.source_cards[1].set_metric(f"{human:,}", "primary evidence class")
        self.source_cards[2].set_metric(f"{microbe:,}", "primary evidence class")
        self.source_cards[3].set_metric(f"{co:,} / {other:,}", "co-metabolic / food, drug, environment")
        columns = [
            "Accepted Compound ID", "Accepted Description", "Formula", "Chemical family (text rule)", "Formula polarity proxy (0-100)",
            "Primary source class", "Source confidence", "Source model entropy (0-1)", "Source class margin (percentage points)",
            "Human endogenous likelihood (%)", "Microbiome-derived likelihood (%)", "Host–microbe co-metabolic likelihood (%)",
            "Food-derived likelihood (%)", "Drug-derived likelihood (%)", "Environmental-derived likelihood (%)", "Source evidence",
            "Dominant predicted phase", "Analyzed-phase likelihood (%)", "Phase confidence", "Phase evidence", "Interpretation limitation",
        ]
        self.source_model.set_frame(self.provenance[[column for column in columns if column in self.provenance.columns]])
        _fit_table(self.source_table, 240)
        self._populate_export_page()
        self._unlock_stage(7)

    def _generate_journal_report(self) -> None:
        if self.analysis is None or not self.filter_settings:
            QMessageBox.warning(self, "Filtered results required", "Apply the Step 06 filters before generating the journal-ready section.")
            return
        analysis = self.analysis
        shortlist = self.shortlist.copy()
        provenance = self.provenance.copy()
        filters = dict(self.filter_settings)
        context = dict(self.provenance_settings) if self.provenance_settings else None
        self._run_background_task(
            "Step 08 · Journal report",
            "Draft dataset-specific methodology, curation analysis, limitations, and reporting checklist",
            lambda progress: journal_report_markdown(analysis, shortlist, provenance, filters, context, progress),
            on_success=self._journal_report_ready,
            outcome=lambda report: f"Journal-ready methodology and analysis generated ({len(report.split()):,} words).",
            failure_title="Journal report generation failed",
        )

    def _journal_report_ready(self, report: str) -> None:
        self.journal_report = report
        self.journal_view.setMarkdown(report)
        section_count = report.count("\n## ")
        self.journal_status.setText(f"Ready · {len(report.split()):,} words · {section_count:,} major sections · verify experimental checklist before submission")
        self._populate_export_page()

    def _copy_journal_report(self) -> None:
        if not self.journal_report:
            QMessageBox.information(self, "Generate the report", "Generate the journal-ready report first.")
            return
        QApplication.clipboard().setText(self.journal_report)
        self._record_instant_task("Step 08 · Journal report", "Copy report Markdown", "Journal-ready report copied to the clipboard.")

    def _save_journal_report(self) -> None:
        if not self.journal_report:
            QMessageBox.information(self, "Generate the report", "Generate the journal-ready report first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save journal-ready methodology and analysis",
            str(Path.home() / "JOURNAL_READY_METHODS_AND_ANALYSIS.md"),
            "Markdown (*.md);;Text (*.txt)",
        )
        if not path:
            return
        try:
            Path(path).write_text(self.journal_report, encoding="utf-8")
            self._record_instant_task("Step 08 · Journal report", "Save journal report", f"Report saved: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Report save failed", str(exc))

    def _stage_mapping(self) -> dict[str, tuple[pd.DataFrame, str]]:
        if self.analysis is None:
            return {}
        return {
            "project_metadata": (pd.DataFrame([(key, value) for key, value in self.analysis.parameters.items() if not isinstance(value, pd.DataFrame)], columns=["Parameter", "Value"]), "00_Project_Metadata.csv"),
            "sample_metadata": (self.sample_model.frame() if self.sample_model.rowCount() else self.analysis.sample_metadata, "01_Sample_Metadata.csv"),
            "audit": (self.analysis.audits, "02_File_Audit.csv"),
            "accepted": (self.analysis.accepted_reads, "03_All_Accepted_Feature_Reads.csv"),
            "duplicates": (self.analysis.duplicate_reads, "04_Repeated_Accepted_ID_Reads.csv"),
            "ms2_spectra": (self.analysis.fd, "05_MS2_Spectra.csv"),
            "ms2_peaks": (self.analysis.fd_peaks, "06_MS2_Fragment_Peaks.csv"),
            "selected": (self.analysis.selected_reads, "07_Selected_Compound_Representatives.csv"),
            "analytical_qc": (self.analysis.analytical_qc, "08_Analytical_QC.csv"),
            "filter_decisions": (self.filter_decisions, "09_Filter_Decisions_All_Entities.csv"),
            "shortlist": (self.shortlist, "10_Filtered_Compound_Shortlist.csv"),
            "provenance": (self.provenance, "11_Chemistry_Source_and_Phase.csv"),
            "thresholds": (threshold_table(self.analysis), "12_Threshold_Method.csv"),
            "task_history": (self._task_history_frame(include_active=True), "14_Task_History.csv"),
        }

    def _export_stage_csv(self, stage: str) -> None:
        mapping = self._stage_mapping()
        if stage not in mapping:
            QMessageBox.warning(self, "No results", "Complete the preceding step before exporting this file.")
            return
        frame, name = mapping[stage]
        if frame.empty and stage not in {"provenance", "filter_decisions", "shortlist"}:
            QMessageBox.warning(self, "No results", "This stage has no rows to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", str(Path.home() / name), "CSV (*.csv)")
        if path:
            export_frame = frame.copy()

            def write_csv(progress: Callable[[int, str], None]) -> Path:
                progress(15, f"Preparing {len(export_frame):,} rows for safe CSV output…")
                result = export_csv(export_frame, path)
                progress(100, f"CSV written: {result}")
                return result

            self._run_background_task(
                "Stage output · CSV",
                f"Export {name}",
                write_csv,
                outcome=lambda result: f"CSV export complete: {result}",
                failure_title="CSV export failed",
            )

    def _populate_export_page(self) -> None:
        if self.analysis is None:
            return
        descriptions = {
            "project_metadata": "Project and analytical settings",
            "sample_metadata": "Sample roles, groups, batches, order, dilution, matrix, extraction",
            "audit": "Hashes, structure, joins, coverage, duplicate keys",
            "accepted": "Every accepted feature read with six-file evidence and intensities",
            "duplicates": "Every read in repeated Accepted-ID groups",
            "ms2_spectra": "Spectrum-level integrity and entropy descriptors",
            "ms2_peaks": "Every fragment m/z, intensity, relative intensity, neutral loss",
            "selected": "One deterministic representative per accepted identity",
            "analytical_qc": "Pooled-QC, blank, detection, drift, D-ratio diagnostics",
            "filter_decisions": "Every representative with rule-level pass/fail reasons",
            "shortlist": "Final compounds passing enabled rules",
            "provenance": "Chemical properties, source evidence, extraction plausibility",
            "thresholds": "Recommended cutoffs, bootstrap intervals, method rationale",
            "task_history": "Background task status, timestamps, duration, and outcome audit",
        }
        rows = []
        for stage, (frame, name) in self._stage_mapping().items():
            ready = "Ready" if not frame.empty or stage in {"provenance", "filter_decisions", "shortlist"} and self.filter_settings else "Pending"
            rows.append({"Output file": name, "Stage": stage, "Status": ready, "Rows": len(frame), "Columns": len(frame.columns), "Contents": descriptions[stage]})
        rows.extend([
            {"Output file": "LCMS_Compound_Curation_Results.xlsx", "Stage": "workbook", "Status": "Ready" if self.filter_settings else "Pending", "Rows": "multi-sheet", "Columns": "multi-sheet", "Contents": "All stage tables plus raw source summaries"},
            {"Output file": "analysis_manifest.json", "Stage": "manifest", "Status": "Ready" if self.filter_settings else "Pending", "Rows": 1, "Columns": "JSON", "Contents": "Counts, settings, hashes, role counts, limitations"},
            {"Output file": "JOURNAL_READY_METHODS_AND_ANALYSIS.md", "Stage": "journal report", "Status": "Ready" if self.journal_report else "Pending", "Rows": "n/a", "Columns": "markdown", "Contents": "Dataset-specific manuscript methodology, curation results, limitations, checklist, and references"},
            {"Output file": "METHOD_AND_LIMITATIONS.txt", "Stage": "method", "Status": "Ready", "Rows": "n/a", "Columns": "text", "Contents": "Scientific logic, citations, and interpretation limits"},
            {"Output file": "SCIENTIFIC_METHOD.md", "Stage": "method", "Status": "Ready", "Rows": "n/a", "Columns": "markdown", "Contents": "Equations, ranking, threshold, QC, MS2, chemistry, source, and extraction logic"},
        ])
        self.export_model.set_frame(pd.DataFrame(rows))
        _fit_table(self.export_table, 380)

    def _export_workbook(self) -> None:
        if self.analysis is None or not self.filter_settings:
            QMessageBox.warning(self, "Results required", "Apply the filters before exporting the workbook.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export complete XLSX", str(Path.home() / "LCMS_Compound_Curation_Results.xlsx"), "Excel workbook (*.xlsx)")
        if not path:
            return
        self._run_export(
            lambda progress: export_workbook(
                self.analysis,
                self.shortlist,
                self.provenance,
                self.filter_settings,
                self.provenance_settings,
                path,
                self.include_candidates.isChecked(),
                task_history=self._task_history_frame(),
                progress=progress,
            ),
            "Building the XLSX workbook…",
        )

    def _export_complete(self) -> None:
        if self.analysis is None or not self.filter_settings:
            QMessageBox.warning(self, "Results required", "Apply the filters before exporting the complete results folder.")
            return
        parent = QFileDialog.getExistingDirectory(self, "Choose parent directory for results", str(Path.home()))
        if not parent:
            return
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        target = Path(parent) / f"{suggested_results_name(self.project_name.text())}_{timestamp}"
        self._run_export(
            lambda progress: export_complete_folder(
                self.analysis,
                self.shortlist,
                self.provenance,
                self.filter_settings,
                self.provenance_settings,
                target,
                self.include_candidates.isChecked(),
                task_history=self._task_history_frame(),
                progress=progress,
            ),
            "Writing the complete results folder…",
        )

    def _run_export(self, function: Callable[[Callable[[int, str], None]], Any], message: str) -> None:
        self._run_background_task(
            "Step 08 · Export",
            message.rstrip("…"),
            function,
            outcome=lambda result: f"Export complete: {result}",
            failure_title="Export failed",
        )

    def _export_succeeded(self, result: Any) -> None:
        self.statusBar().showMessage(f"Export complete: {result}")
        QMessageBox.information(self, "Export complete", f"Saved successfully:\n{result}")

    def _export_failed(self, trace: str) -> None:
        final_line = trace.strip().splitlines()[-1] if trace.strip() else "Unknown error"
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Icon.Critical)
        box.setWindowTitle("Export failed")
        box.setText(final_line)
        box.setDetailedText(trace)
        box.exec()

    def _forget_worker(self, worker: QThread) -> None:
        if worker in self._workers:
            self._workers.remove(worker)
        worker.deleteLater()

    def _project_payload(self) -> dict[str, Any]:
        metadata = None
        if self.sample_model.rowCount():
            metadata = self.sample_model.frame().replace({np.nan: None}).to_dict(orient="records")
        return {
            "application": APP_NAME,
            "version": APP_VERSION,
            "saved_at": datetime.now().isoformat(),
            "files": {role: self.path_edits[role].text().strip() for role in FILE_ROLES},
            "parameters": {key: value for key, value in self._collect_parameters().items() if not isinstance(value, pd.DataFrame)},
            "sample_metadata": metadata,
            "filters": self.filter_settings or None,
            "provenance_settings": self.provenance_settings,
        }

    def _save_project(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save project configuration", str(Path.home() / "LCMS_Curation_Project.json"), "JSON project (*.json)")
        if not path:
            return
        try:
            Path(path).write_text(
                json.dumps(self._project_payload(), indent=2, ensure_ascii=False, default=lambda value: value.item() if isinstance(value, np.generic) else str(value)),
                encoding="utf-8",
            )
            self._record_instant_task("Project", "Save project configuration", f"Project configuration saved: {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Project save failed", str(exc))

    def _load_project(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load project configuration", str(Path.home()), "JSON project (*.json);;All files (*.*)")
        if not path:
            return
        try:
            payload = json.loads(Path(path).read_text(encoding="utf-8"))
            for role, value in payload.get("files", {}).items():
                if role in self.path_edits:
                    self.path_edits[role].setText(str(value or ""))
            parameters = payload.get("parameters", {})
            self.project_name.setText(str(parameters.get("project_name", "LCMS Curation Project")))
            self.analyst.setText(str(parameters.get("analyst", "")))
            _set_combo_data(self.assay_type, parameters.get("assay_type", "metabolomics"))
            _set_combo_data(self.biological_system, parameters.get("biological_system", "human"))
            self.sample_matrix.setCurrentText(str(parameters.get("sample_matrix", "")))
            _set_combo_data(self.extraction_method, parameters.get("extraction_method", "unknown"))
            _set_combo_data(self.analyzed_phase, parameters.get("analyzed_phase", "unknown"))
            _set_combo_data(self.qc_pooling, parameters.get("qc_pooling_point", "unknown"))
            _set_combo_data(self.ion_mode, parameters.get("ion_mode", "negative"))
            _set_combo_data(self.acquisition_mode, parameters.get("acquisition_mode", "unknown"))
            self.mass_tolerance.setValue(float(parameters.get("mass_tolerance", 10.0)))
            self.qc_pattern.setText(str(parameters.get("qc_pattern", "QC|Pool")))
            metadata = payload.get("sample_metadata")
            self.pending_sample_metadata = pd.DataFrame(metadata) if metadata else None
            self._record_instant_task("Project", "Load project configuration", "Project loaded; verify paths and scan all six files.")
        except Exception as exc:
            QMessageBox.critical(self, "Project load failed", str(exc))

    def _show_about(self) -> None:
        QMessageBox.about(
            self,
            "Method and version",
            f"<b>{APP_NAME}</b><br>Version {APP_VERSION}<br><br>"
            "Portable PySide6 desktop workbench for complete six-file reconciliation, editable sample/QC metadata, peak-level MS² extraction, deterministic identity resolution, analytical-QC review, adaptive confidence triage, formula chemistry, contextual source/extraction interpretation, background task auditing, and dataset-specific journal reporting.<br><br>"
            "Primary references and URLs are written into every complete export. Thresholds and contextual likelihoods are review aids—not identification FDR, Level-1 confirmation, causal probabilities, or measured recovery.",
        )


def run() -> int:
    if sys.platform.startswith("win"):
        os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
        os.environ.setdefault("QT_SCALE_FACTOR_ROUNDING_POLICY", "PassThrough")
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setOrganizationName("SiddharthSinghResearch")
    app.setFont(QFont("Segoe UI" if sys.platform.startswith("win") else "DejaVu Sans", 10))
    apply_accessible_light_theme(app)
    window = MainWindow()
    window.show()
    return app.exec()
