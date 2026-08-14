@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
set "QT_STYLE_OVERRIDE=Fusion"
title LC-MS Compound Curation Workbench 2.0 Launcher
where py >nul 2>nul
if not errorlevel 1 (
    py -3 install_and_run.py
) else (
    where python >nul 2>nul
    if not errorlevel 1 (
        python install_and_run.py
    ) else (
        where python3 >nul 2>nul
        if not errorlevel 1 (
            python3 install_and_run.py
        ) else (
            echo Python was not found.
            echo Install 64-bit Python 3.11 through 3.14 from https://www.python.org/downloads/windows/
            echo During installation, enable "Add python.exe to PATH", then run this file again.
            pause
            exit /b 2
        )
    )
)
if errorlevel 1 (
    echo.
    echo The application did not start. Review the message above.
    echo If dependencies were interrupted, run: py -3 install_and_run.py --repair
    pause
)
endlocal
