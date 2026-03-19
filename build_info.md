# LOREON — Build & Distribution Guide

**Version:** 2.6
**Last updated:** 2026-03-18

---

## Architettura di distribuzione

LOREON utilizza un'architettura ibrida per la distribuzione:

1. **PyInstaller** impacchetta solo la GUI (`pipeline_gui.py` + PyQt5) in un eseguibile nativo
2. **Al primo avvio**, l'applicazione crea automaticamente un **virtualenv** con tutte le dipendenze Python necessarie per la pipeline (pandas, duckdb, plotly, etc.)
3. Il virtualenv viene creato una sola volta e riutilizzato per tutti gli avvii successivi

### Perche' questo approccio

- **L'utente deve installare solo Python 3.8+** (e Docker su Windows)
- PyInstaller non riesce a impacchettare correttamente librerie complesse come duckdb/plotly/biopython
- Il virtualenv garantisce compatibilita' con il Python di sistema dell'utente
- Gli aggiornamenti delle dipendenze non richiedono un rebuild del pacchetto

### Flusso al primo avvio

```
Utente apre LOREON
  → pipeline_gui.py: _ensure_venv() verifica se il venv esiste
  → Se NO: mostra progress dialog, chiama venv_manager.create_venv()
    → Trova python3 di sistema
    → Crea venv in ~/Library/Application Support/LOREON/venv/ (macOS)
                    ~/.local/share/loreon/venv/ (Linux)
    → pip install -r requirements.txt
  → Se SI': avvio immediato
```

### Dove viene salvato il virtualenv

| Piattaforma | Percorso |
|---|---|
| macOS | `~/Library/Application Support/LOREON/venv/` |
| Linux | `~/.local/share/loreon/venv/` |
| Windows | `%LOCALAPPDATA%\LOREON\venv\` |

### File chiave

| File | Ruolo |
|---|---|
| `venv_manager.py` | Creazione/gestione del virtualenv, discovery di Python di sistema |
| `pipeline_worker.py` | Usa `get_venv_python()` per lanciare gli script con il Python del venv |
| `pipeline_gui.py` | Chiama `_ensure_venv()` prima di mostrare la finestra principale |

---

## Prerequisiti utente finale

| Piattaforma | Requisiti |
|---|---|
| **macOS** | Python 3.8+ (`brew install python@3` o python.org) |
| **Linux** | Python 3.8+, python3-venv, python3-pip (dichiarati come dipendenze .deb) |
| **Windows** | Python 3.8+ (python.org) + Docker Desktop |

Connessione internet richiesta al primo avvio per scaricare le dipendenze pip.

---

## Build Scripts

| Script | Platform | Output |
|---|---|---|
| `build_dmg.sh` | macOS | `dist/LOREON-2.6-macOS.dmg` |
| `build_deb.sh` | Ubuntu/Debian | `dist/loreon_2.6_amd64.deb` |

---

## macOS — DMG Installer (`build_dmg.sh`)

### Requisiti build machine

- Python 3.8+ con pip
- PyQt5 (`pip install PyQt5`)
- PyInstaller (`pip install pyinstaller`)
- Xcode Command Line Tools (per `sips`, `iconutil`, `hdiutil`)

### Build

```bash
chmod +x build_dmg.sh
./build_dmg.sh
```

### Output

- `dist/LOREON.app` — App bundle
- `dist/LOREON-2.6-macOS.dmg` — Disk image distribuibile

### Contenuto del DMG

- `LOREON.app` — Applicazione (PyInstaller --onedir + venv_manager)
- `Applications` — Link simbolico a `/Applications` per drag-and-drop install
- `README.txt` — Istruzioni per l'utente

### Caratteristiche

- **Bundle identifier:** `com.loreon.pipeline`
- **Icone Retina** complete (16px → 1024px, incluse varianti @2x)
- **Info.plist** aggiornato con versione e bundle ID
- **Compressione UDZO** ottimizzata (zlib-level=9)
- L'app **non e' code-signed**
- **Primo avvio**: crea automaticamente il virtualenv (~1-2 min)

### Installazione (utente finale)

1. Installare Python 3.8+ se non presente
2. Aprire il file `.dmg`
3. Trascinare `LOREON.app` nella cartella `Applications`
4. Al primo avvio: click destro → Apri → Apri
5. Oppure da Terminale:
   ```bash
   xattr -d com.apple.quarantine /Applications/LOREON.app
   ```

### Dipendenze runtime

```bash
brew install minimap2 samtools
# oppure
conda install -c bioconda minimap2 samtools
```

---

## Ubuntu/Debian — Pacchetto .deb (`build_deb.sh`)

### Requisiti build machine

- Ubuntu 20.04+ / Debian 11+
- Python 3.8+ con pip
- `dpkg-deb` (di solito preinstallato)
- `fakeroot` (`sudo apt install fakeroot`)
- PyQt5: `pip install PyQt5` oppure `sudo apt install python3-pyqt5`
- PyInstaller: `pip install pyinstaller`

### Build

```bash
chmod +x build_deb.sh
./build_deb.sh
```

### Output

- `dist/loreon_2.6_amd64.deb` — Pacchetto Debian installabile

### Struttura del pacchetto

| Percorso | Contenuto |
|---|---|
| `/opt/loreon/` | Eseguibile PyInstaller + moduli Python + dati |
| `/usr/local/bin/loreon` | Symlink → GUI (`/opt/loreon/LOREON`) |
| `/usr/local/bin/loreon-cli` | Symlink → CLI (`/opt/loreon/metaGenomics_new.py`) |
| `/usr/share/applications/loreon.desktop` | Integrazione menu applicazioni |
| `/usr/share/icons/hicolor/*/apps/loreon.png` | Icone (48, 128, 256px) |
| `/usr/share/pixmaps/loreon.png` | Icona fallback (legacy) |
| `/usr/share/man/man1/loreon.1.gz` | Man page |
| `/usr/share/doc/loreon/` | Crediti e licenza MIT |

### Dipendenze dichiarate nel pacchetto

| Tipo | Pacchetti |
|---|---|
| **Depends** (obbligatori) | `python3 (>= 3.8)`, `python3-venv`, `python3-pip`, `libxcb-xinerama0`, `libxcb-cursor0`, `libgl1`, `libglib2.0-0`, `libfontconfig1`, `libxkbcommon0` |
| **Recommends** | `samtools` |
| **Suggests** | `minimap2` |

### Script di manutenzione

- **postinst:** crea symlink in `/usr/local/bin/`, aggiorna desktop database e icon cache
- **prerm:** rimuove i symlink prima della disinstallazione
- **postrm:** aggiorna desktop/icon cache; su `purge` rimuove `/opt/loreon/`

### Installazione (utente finale)

```bash
# Metodo 1: apt (risolve dipendenze automaticamente)
sudo apt install ./loreon_2.6_amd64.deb

# Metodo 2: dpkg + apt
sudo dpkg -i loreon_2.6_amd64.deb
sudo apt-get install -f    # risolve dipendenze mancanti
```

### Utilizzo dopo installazione

```bash
loreon              # Avvia GUI
loreon-cli --help   # Modalita' CLI
man loreon          # Documentazione
```

### Disinstallazione

```bash
sudo dpkg -r loreon        # rimuove il pacchetto
sudo dpkg -P loreon        # rimuove + purge (elimina /opt/loreon/ e venv)
```

### Dipendenze runtime

```bash
# minimap2 e samtools devono essere nel PATH
conda install -c bioconda minimap2 samtools
# oppure (solo samtools):
sudo apt install samtools
```

---

## File inclusi in entrambi i pacchetti

| File | Ruolo |
|---|---|
| `pipeline_gui.py` | Entry point GUI (PyQt5) — impacchettato da PyInstaller |
| `pipeline_worker.py` | Worker thread per esecuzione pipeline |
| `venv_manager.py` | Gestione virtualenv: creazione, verifica, discovery Python |
| `metaGenomics_new.py` | Entry point CLI (argparse) |
| `OtuUtils.py` | Utilita' core: allineamento, filtro, tassonomia |
| `ResultsReader.py` | Generazione OTU table con DuckDB |
| `report_generator.py` | Generatore report HTML (Plotly + Jinja2) |
| `pipeline_profiler.py` | Profilazione performance |
| `PipelineLogger.py` | Logging esecuzione |
| `fileorganizer.py` | Utilita' organizzazione file |
| `template.html` | Template Jinja2 per report HTML |
| `Dockerfile` | Container per modalita' Docker (Windows) |
| `requirements.txt` | Dipendenze Python (installate nel venv al primo avvio) |
| `loreon_app_icon_2.png` | Icona applicazione |
| `loreon.jpeg` | Immagine splash screen |
| `credits.md` | Crediti |

---

## Note

- Il **rilevamento architettura** nel build .deb e' automatico (`dpkg --print-architecture`)
- Se **ImageMagick** (`convert`) e' disponibile sulla build machine, le icone vengono ridimensionate correttamente; altrimenti viene copiata l'immagine originale
- Il pacchetto .deb viene verificato con **lintian** se disponibile
- Entrambi gli script puliscono i file temporanei al termine del build
- La versione viene letta da `CITATION.cff` (attualmente `2.6`); aggiornare la variabile `APP_VERSION` negli script quando si rilascia una nuova versione
- **PyInstaller impacchetta solo PyQt5** — tutte le altre dipendenze (pandas, duckdb, plotly, biopython, etc.) vengono installate via pip nel virtualenv al primo avvio
- Il virtualenv viene creato nella home dell'utente, non in `/opt/loreon/`, per evitare problemi di permessi
- Se il virtualenv viene corrotto o cancellato, viene ricreato automaticamente al prossimo avvio
