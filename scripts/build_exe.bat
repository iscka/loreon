@echo off
setlocal
cd /d "%~dp0.."

if not exist "sources\loreon_app_icon_2.png" (
    echo ERROR: sources\loreon_app_icon_2.png non trovato.
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
python -c "from PIL import Image; img = Image.open('sources/loreon_app_icon_2.png'); img.save('sources/loreon_app_icon_2.ico', sizes=[(16,16),(32,32),(48,48),(64,64),(128,128),(256,256)])"
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
    --icon sources\loreon_app_icon_2.ico ^
    --hidden-import=_socket ^
    --hidden-import=select ^
    --hidden-import=PyQt5.sip ^
    --hidden-import=encodings.utf_8 ^
    --hidden-import=encodings.ascii ^
    --collect-all=PyQt5 ^
    --paths sources ^
    --clean ^
    -y ^
    sources\pipeline_gui.py

if errorlevel 1 (
    echo ERROR: Build fallita.
    pause
    exit /b 1
)

echo Copia file nella cartella dist...
copy sources\Dockerfile dist\LOREON\
copy sources\requirements.txt dist\LOREON\
if exist sources\.dockerignore copy sources\.dockerignore dist\LOREON\
copy sources\loreon_app_icon_2.png dist\LOREON\
copy sources\loreon_app_icon_2.ico dist\LOREON\
copy sources\loreon.jpeg dist\LOREON\
copy sources\credits.md dist\LOREON\
copy LICENSE dist\LOREON\
copy sources\metaGenomics_new.py dist\LOREON\
copy sources\OtuUtils.py dist\LOREON\
copy sources\ResultsReader.py dist\LOREON\
copy sources\PipelineLogger.py dist\LOREON\
copy sources\report_generator.py dist\LOREON\
copy sources\pipeline_profiler.py dist\LOREON\
copy sources\template.html dist\LOREON\
xcopy /E /I /Y sources\example_data dist\LOREON\example_data

echo.
echo BUILD COMPLETATA!
echo Eseguibile: dist\LOREON\LOREON.exe
pause
