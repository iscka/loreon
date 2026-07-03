@echo off
setlocal enabledelayedexpansion

echo ============================================================
echo  LOREON Installer Builder
echo ============================================================
echo.

REM ── Ricerca ISCC.exe ────────────────────────────────────────
set ISCC=

REM 1) Già in PATH?
where ISCC.exe >nul 2>&1
if %ERRORLEVEL%==0 (
    for /f "delims=" %%i in ('where ISCC.exe') do set "ISCC=%%i"
    goto :found
)

REM 2) Percorsi standard Program Files
for %%P in (
    "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
    "%ProgramFiles%\Inno Setup 6\ISCC.exe"
    "%LocalAppData%\Programs\Inno Setup 6\ISCC.exe"
    "C:\InnoSetup\ISCC.exe"
) do (
    if exist %%P (
        set "ISCC=%%~P"
        goto :found
    )
)

REM 3) Registro di sistema (InstallLocation dal setup di Inno Setup)
for %%K in (
    "HKLM\SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1"
    "HKLM\SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1"
) do (
    for /f "tokens=2*" %%a in (
        'reg query %%K /v InstallLocation 2^>nul ^| findstr InstallLocation'
    ) do (
        if exist "%%b\ISCC.exe" (
            set "ISCC=%%b\ISCC.exe"
            goto :found
        )
    )
)

REM 4) Non trovato — chiedi il percorso manualmente
echo [ERRORE] Inno Setup 6 non trovato nei percorsi standard.
echo.
echo Percorsi cercati:
echo   - In PATH di sistema
echo   - %ProgramFiles(x86)%\Inno Setup 6\
echo   - %ProgramFiles%\Inno Setup 6\
echo   - %LocalAppData%\Programs\Inno Setup 6\
echo   - Registro di sistema
echo.
echo Inserisci il percorso completo di ISCC.exe (o premi INVIO per aprire il download):
set /p ISCC="> "

if "!ISCC!"=="" (
    start https://jrsoftware.org/isdl.php
    echo Dopo l'installazione, riavvia questo script.
    pause
    exit /b 1
)

if not exist "!ISCC!" (
    echo [ERRORE] File non trovato: !ISCC!
    pause
    exit /b 1
)

:found
echo [OK] Inno Setup: %ISCC%
echo.

REM ── Verifica dist\LOREON ────────────────────────────────────
if not exist "dist\LOREON\LOREON.exe" (
    echo [ERRORE] dist\LOREON\LOREON.exe non trovato.
    echo Esegui prima build_exe.bat per compilare l'eseguibile.
    pause
    exit /b 1
)
echo [OK] dist\LOREON\ trovata.
echo.

REM ── Verifica icona ──────────────────────────────────────────
if not exist "loreon_app_icon_2.ico" (
    echo [ATTENZIONE] loreon_app_icon_2.ico non trovata.
    echo             Per l'icona corretta, esegui prima build_exe.bat
    echo.
)

REM ── Directory di output ─────────────────────────────────────
if not exist "installer_output" mkdir windows_setup

REM ── Compilazione ────────────────────────────────────────────
echo Compilazione installer in corso...
echo.
"%ISCC%" installer.iss

if %ERRORLEVEL% neq 0 (
    echo.
    echo [ERRORE] Compilazione fallita (codice %ERRORLEVEL%).
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ============================================================
echo  Installer creato: installer_output\LOREON_Setup_1.0.exe
echo ============================================================
echo.
pause
