; ============================================================
;  LOREON v2.6 - Inno Setup 6 Installer Script
;  Requires: Inno Setup 6.3+  (testato su 6.7.1)
; ============================================================

#define AppName      "LOREON"
#define AppVersion   "2.6"
#define AppPublisher "University of Perugia"
#define AppURL       "https://github.com/iscka/loreon"
#define AppExeName   "LOREON.exe"

; Minimum Python version required
#define PyMinMajor 3
#define PyMinMinor 14

; ── Setup ────────────────────────────────────────────────────
[Setup]
AppId={{F3A2D1B0-C4E5-4F67-89AB-0123456789CD}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}
AppUpdatesURL={#AppURL}

DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
AllowNoIcons=yes

OutputDir=windows_setup
OutputBaseFilename=LOREON_Setup_{#AppVersion}
SetupIconFile=loreon_app_icon_2.ico
UninstallDisplayIcon={app}\loreon_app_icon_2.ico

Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
WizardSizePercent=120
MinVersion=10.0
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequired=admin
PrivilegesRequiredOverridesAllowed=dialog
DisableProgramGroupPage=yes

; ── Languages ────────────────────────────────────────────────
[Languages]
Name: "italian"; MessagesFile: "compiler:Languages\Italian.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"

; ── Tasks (optional shortcuts) ───────────────────────────────
[Tasks]
Name: "desktopicon";  Description: "Crea collegamento sul Desktop"; GroupDescription: "Icone aggiuntive:"
Name: "startupicon";  Description: "Avvia LOREON con Windows";               GroupDescription: "Avvio automatico:";  Flags: unchecked

; ── Files ────────────────────────────────────────────────────
[Files]
; Tutto il contenuto della directory dist\LOREON\
Source: "dist\LOREON\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

; ── Icons ────────────────────────────────────────────────────
[Icons]
Name: "{group}\{#AppName}";              Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\loreon_app_icon_2.ico"
Name: "{group}\Disinstalla {#AppName}";  Filename: "{uninstallexe}"
Name: "{commondesktop}\{#AppName}";      Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\loreon_app_icon_2.ico"; Tasks: desktopicon
Name: "{commonstartup}\{#AppName}";      Filename: "{app}\{#AppExeName}"; Tasks: startupicon

; ── Post-install actions ─────────────────────────────────────
[Run]
; Installa dipendenze Python (solo se Python è disponibile)
Filename: "{cmd}"; Parameters: "/C python -m pip install -r ""{app}\requirements.txt"" --quiet 2>&1"; WorkingDir: "{app}"; Flags: runhidden; StatusMsg: "Installazione dipendenze Python in corso..."; Check: IsPythonAvailable

; Avvia l'applicazione al termine dell'installazione
Filename: "{app}\{#AppExeName}"; Description: "Avvia {#AppName} {#AppVersion}"; Flags: nowait postinstall skipifsilent

; ── Uninstall cleanup ─────────────────────────────────────────
[UninstallDelete]
Type: filesandordirs; Name: "{app}\__pycache__"
Type: filesandordirs; Name: "{app}\_internal\__pycache__"

; ── Pascal code section ──────────────────────────────────────
[Code]

{ ----------------------------------------------------------- }
{  Global state                                               }
{ ----------------------------------------------------------- }
const
  NL = #13#10;   { CRLF — defined here so #13#10 never starts a line }

var
  GDockerOK  : Boolean;
  GPythonOK  : Boolean;
  GPythonVer : String;   { e.g. "3.14.1" }

{ ----------------------------------------------------------- }
{  Utility: run a command and capture stdout+stderr           }
{ ----------------------------------------------------------- }
function ExecCapture(const Cmd, Args: String; var Output: String): Boolean;
var
  TmpFile   : String;
  ResultCode: Integer;
  AnsiOut   : AnsiString;
begin
  TmpFile := ExpandConstant('{tmp}\loreon_check_output.txt');
  Result  := Exec(ExpandConstant('{cmd}'),
                  '/C ' + Cmd + ' ' + Args + ' > "' + TmpFile + '" 2>&1',
                  '', SW_HIDE, ewWaitUntilTerminated, ResultCode);
  Output := '';
  if FileExists(TmpFile) then begin
    LoadStringFromFile(TmpFile, AnsiOut);
    Output := String(AnsiOut);
  end;
end;

{ ----------------------------------------------------------- }
{  Docker: check registry then filesystem                     }
{ ----------------------------------------------------------- }
function CheckDockerInstalled: Boolean;
var
  Dummy: String;
begin
  { HKLM 64-bit }
  Result := RegQueryStringValue(HKLM,
    'SOFTWARE\Docker Inc.\Docker Desktop', 'InstallDir', Dummy);
  if Result then Exit;

  { HKLM 32-bit node (WOW6432Node) }
  Result := RegKeyExists(HKLM,
    'SOFTWARE\WOW6432Node\Docker Inc.\Docker Desktop');
  if Result then Exit;

  { Fallback: exe on disk }
  Result := FileExists('C:\Program Files\Docker\Docker\Docker Desktop.exe');
end;

{ ----------------------------------------------------------- }
{  Python: parse "Python X.Y.Z" from --version output        }
{ ----------------------------------------------------------- }
function ParsePythonVersion(const Output: String;
                             out Major, Minor: Integer;
                             out Full: String): Boolean;
var
  S: String;
  P: Integer;
begin
  Result := False;
  Major  := 0;
  Minor  := 0;
  Full   := '';

  P := Pos('Python ', Output);
  if P = 0 then Exit;

  S := Trim(Copy(Output, P + 7, Length(Output)));

  { strip anything after the version digits }
  Full := S;

  { major }
  P := Pos('.', S);
  if P = 0 then Exit;
  Major := StrToIntDef(Copy(S, 1, P - 1), -1);
  if Major < 0 then Exit;
  S := Copy(S, P + 1, Length(S));

  { minor }
  P := Pos('.', S);
  if P > 0 then
    S := Copy(S, 1, P - 1);   { drop patch }
  Minor := StrToIntDef(Trim(S), -1);
  if Minor < 0 then Exit;

  Result := True;
end;

function CheckPythonVersion: Boolean;
var
  Output     : String;
  Major, Minor: Integer;
  Full        : String;
begin
  Result := False;

  { Try "python" first, then "python3" }
  if not ExecCapture('python', '--version', Output) then
    if not ExecCapture('python3', '--version', Output) then
      Exit;

  if not ParsePythonVersion(Output, Major, Minor, Full) then Exit;

  GPythonVer := Full;
  Result := (Major > {#PyMinMajor}) or
            ((Major = {#PyMinMajor}) and (Minor >= {#PyMinMinor}));
end;

{ ----------------------------------------------------------- }
{  Exposed check used in [Run] section                        }
{ ----------------------------------------------------------- }
function IsPythonAvailable: Boolean;
begin
  Result := GPythonOK;
end;

{ ----------------------------------------------------------- }
{  Open a URL in the default browser                          }
{ ----------------------------------------------------------- }
procedure OpenURL(const URL: String);
var
  ResultCode: Integer;
begin
  ShellExec('open', URL, '', '', SW_SHOWNORMAL, ewNoWait, ResultCode);
end;

{ ----------------------------------------------------------- }
{  Build the "missing components" message and offer downloads }
{ ----------------------------------------------------------- }
procedure CheckPrerequisites;
var
  Msg       : String;
  MsgResult : Integer;
  OpenDocker: Boolean;
  OpenPython: Boolean;
begin
  GDockerOK := CheckDockerInstalled;
  GPythonOK := CheckPythonVersion;

  if GDockerOK and GPythonOK then Exit;   { all good, skip dialog }

  { ---- Build message ---- }
  Msg := 'LOREON richiede i seguenti componenti, che non sono stati rilevati sul sistema:' + NL + NL;

  if not GDockerOK then
    Msg := Msg + '  [MANCANTE]  Docker Desktop' + NL;

  if not GPythonOK then begin
    if GPythonVer <> '' then
      Msg := Msg + '  [AGGIORNARE]  Python ' + GPythonVer +
                   '  (richiesta >= {#PyMinMajor}.{#PyMinMinor})' + NL
    else
      Msg := Msg + '  [MANCANTE]  Python >= {#PyMinMajor}.{#PyMinMinor}' + NL;
  end;

  Msg := Msg + NL +
         'L''installazione di LOREON può proseguire, ma l''applicazione' + NL +
         'non funzionerà correttamente senza i componenti indicati.' + NL + NL +
         'Vuoi aprire le pagine di download dei componenti mancanti?';

  MsgResult := MsgBox(Msg, mbConfirmation, MB_YESNO or MB_DEFBUTTON1);

  if MsgResult = IDYES then begin
    OpenDocker := not GDockerOK;
    OpenPython := not GPythonOK;

    if OpenDocker then
      OpenURL('https://www.docker.com/products/docker-desktop/');

    if OpenPython then begin
      Sleep(400);   { evita che il browser sovrapponga le schede }
      OpenURL('https://www.python.org/downloads/');
    end;

    { Chiedi se l'utente ha completato le installazioni e vuole ri-verificare }
    Sleep(1000);
    MsgResult := MsgBox(
      'Dopo aver installato i componenti mancanti, clicca OK per' + NL +
      'ri-verificare oppure Annulla per proseguire comunque.',
      mbConfirmation, MB_OKCANCEL);

    if MsgResult = IDOK then begin
      { re-check }
      GDockerOK := CheckDockerInstalled;
      GPythonOK := CheckPythonVersion;

      if GDockerOK and GPythonOK then
        MsgBox('Tutti i componenti sono stati rilevati correttamente. Ottimo!',
               mbInformation, MB_OK);
    end;
  end;
end;

{ ----------------------------------------------------------- }
{  Custom status line on wizard overview page                 }
{ ----------------------------------------------------------- }
function GetDockerStatus(Param: String): String;
begin
  if GDockerOK then Result := 'Docker Desktop: OK'
               else Result := 'Docker Desktop: NON TROVATO';
end;

function GetPythonStatus(Param: String): String;
begin
  if GPythonOK then
    Result := 'Python ' + GPythonVer + ': OK'
  else if GPythonVer <> '' then
    Result := 'Python ' + GPythonVer + ': VERSIONE INSUFFICIENTE'
  else
    Result := 'Python >= {#PyMinMajor}.{#PyMinMinor}: NON TROVATO';
end;

{ ----------------------------------------------------------- }
{  Wizard events                                              }
{ ----------------------------------------------------------- }
procedure InitializeWizard;
begin
  GDockerOK  := False;
  GPythonOK  := False;
  GPythonVer := '';
  CheckPrerequisites;
end;

function NextButtonClick(CurPageID: Integer): Boolean;
begin
  Result := True;
  { On the first page, remind if still missing components }
  if CurPageID = wpWelcome then begin
    if not GDockerOK or not GPythonOK then begin
      if MsgBox('Alcuni componenti sono ancora mancanti.' + NL +
                'Vuoi procedere comunque con l''installazione?',
                mbConfirmation, MB_YESNO) = IDNO then
        Result := False;
    end;
  end;
end;
