@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (
    py -3 build_windows_portable.py
) else (
    python build_windows_portable.py
)
pause
endlocal

