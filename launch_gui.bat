@echo off
setlocal

:: Check Python
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python non trovato. Installalo da https://www.python.org/downloads/
    pause
    exit /b 1
)

:: Check / install PyQt5
python -c "import PyQt5" >nul 2>&1
if errorlevel 1 (
    echo Installazione dipendenze GUI...
    python -m pip install PyQt5>=5.15.0 --quiet
    if errorlevel 1 (
        echo ERROR: Installazione PyQt5 fallita.
        pause
        exit /b 1
    )
)

:: Launch GUI from its own directory
cd /d "%~dp0"
python pipeline_gui.py

endlocal
