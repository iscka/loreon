; ============================================================
;  LOREON - Inno Setup 6 Installer Script
;  Requires: Inno Setup 6.3+  (testato su 6.7.1)
; ============================================================

#define AppName      "LOREON"
#define AppVersion   "1.1.1"
#define AppPublisher "University of Perugia"
#define AppURL       "https://github.com/iscka/loreon"
#define AppExeName   "LOREON.exe"

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

OutputDir=..\windows_setup
OutputBaseFilename=LOREON_Setup_{#AppVersion}
SetupIconFile=..\sources\loreon_app_icon_2.ico
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
; Tutto il contenuto della directory dist\LOREON\, nella radice del repository
Source: "..\dist\LOREON\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

; ── Icons ────────────────────────────────────────────────────
[Icons]
Name: "{group}\{#AppName}";              Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\loreon_app_icon_2.ico"
Name: "{group}\Disinstalla {#AppName}";  Filename: "{uninstallexe}"
Name: "{commondesktop}\{#AppName}";      Filename: "{app}\{#AppExeName}"; IconFilename: "{app}\loreon_app_icon_2.ico"; Tasks: desktopicon
Name: "{commonstartup}\{#AppName}";      Filename: "{app}\{#AppExeName}"; Tasks: startupicon

; ── Post-install actions ─────────────────────────────────────
[Run]
; Nessuna dipendenza Python da installare: LOREON.exe include la GUI e la pipeline gira in Docker
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
  GDockerOK : Boolean;

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
{  Open a URL in the default browser                          }
{ ----------------------------------------------------------- }
procedure OpenURL(const URL: String);
var
  ResultCode: Integer;
begin
  ShellExec('open', URL, '', '', SW_SHOWNORMAL, ewNoWait, ResultCode);
end;

{ ----------------------------------------------------------- }
{  Docker missing: explain and offer the download             }
{ ----------------------------------------------------------- }
procedure CheckPrerequisites;
var
  Msg       : String;
  MsgResult : Integer;
begin
  GDockerOK := CheckDockerInstalled;

  if GDockerOK then Exit;   { all good, skip dialog }

  { ---- Build message ---- }
  Msg := 'LOREON richiede Docker Desktop, che non è stato rilevato sul sistema.' + NL + NL +
         'L''installazione di LOREON può proseguire, ma l''applicazione' + NL +
         'non funzionerà correttamente senza Docker Desktop.' + NL + NL +
         'Vuoi aprire la pagina di download?';

  MsgResult := MsgBox(Msg, mbConfirmation, MB_YESNO or MB_DEFBUTTON1);

  if MsgResult = IDYES then begin
    OpenURL('https://www.docker.com/products/docker-desktop/');

    { Chiedi se l'utente ha completato l'installazione e vuole ri-verificare }
    Sleep(1000);
    MsgResult := MsgBox(
      'Dopo aver installato Docker Desktop, clicca OK per' + NL +
      'ri-verificare oppure Annulla per proseguire comunque.',
      mbConfirmation, MB_OKCANCEL);

    if MsgResult = IDOK then begin
      { re-check }
      GDockerOK := CheckDockerInstalled;

      if GDockerOK then
        MsgBox('Docker Desktop è stato rilevato correttamente. Ottimo!',
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

{ ----------------------------------------------------------- }
{  Wizard events                                              }
{ ----------------------------------------------------------- }
procedure InitializeWizard;
begin
  GDockerOK := False;
  CheckPrerequisites;
end;

function NextButtonClick(CurPageID: Integer): Boolean;
begin
  Result := True;
  { On the first page, remind if Docker is still missing }
  if CurPageID = wpWelcome then begin
    if not GDockerOK then begin
      if MsgBox('Docker Desktop è ancora mancante.' + NL +
                'Vuoi procedere comunque con l''installazione?',
                mbConfirmation, MB_YESNO) = IDNO then
        Result := False;
    end;
  end;
end;
