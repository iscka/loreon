@echo off
setlocal
cd /d "%~dp0"

if not exist "loreon_app_icon_2.png" (
    echo ERROR: loreon_app_icon_2.png non trovato in questa cartella.
    pause
    exit /b 1
)

echo Installazione dipendenze build...
python -m pip install pyinstaller pillow --quiet
if errorlevel 1 (
    echo ERROR: Installazione dipendenze fallita.
    pause
    exit /b 1
)

echo Conversione PNG in ICO...
python -c "from PIL import Image; img = Image.open('loreon_app_icon_2.png'); img.save('loreon_app_icon_2.ico', sizes=[(16,16),(32,32),(48,48),(64,64),(128,128),(256,256)])"
if errorlevel 1 (
    echo ERROR: Conversione icona fallita.
    pause
    exit /b 1
)

echo Build in corso...
python -m PyInstaller ^
    --onedir ^
    --windowed ^
    --name LOREON ^
    --icon loreon_app_icon_2.ico ^
    --hidden-import=_socket ^
    --hidden-import=select ^
    --hidden-import=PyQt5.sip ^
    --hidden-import=encodings.utf_8 ^
    --hidden-import=encodings.ascii ^
    --collect-all=PyQt5 ^
    --clean ^
    -y ^
    pipeline_gui.py

if errorlevel 1 (
    echo ERROR: Build fallita.
    pause
    exit /b 1
)

echo Copia file nella cartella dist...
copy Dockerfile dist\LOREON\
copy requirements.txt dist\LOREON\
if exist .dockerignore copy .dockerignore dist\LOREON\
copy loreon_app_icon_2.png dist\LOREON\
copy loreon.jpeg dist\LOREON\
copy credits.md dist\LOREON\
copy metaGenomics_new.py dist\LOREON\
copy OtuUtils.py dist\LOREON\
copy ResultsReader.py dist\LOREON\
copy PipelineLogger.py dist\LOREON\
copy report_generator.py dist\LOREON\
copy pipeline_profiler.py dist\LOREON\
copy template.html dist\LOREON\

echo.
echo BUILD COMPLETATA!
echo Eseguibile: dist\LOREON\LOREON.exe
pause
