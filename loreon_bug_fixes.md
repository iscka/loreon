# LOREON — Bug Fixes & Analisi del Codice
**Sessione di revisione: 2026-02-26**
**Analista: Claude (Anthropic) — revisione statica, nessuna esecuzione runtime**

---

## Indice

1. [Contesto dell'analisi](#1-contesto-dellanalisi)
2. [Mappa del codice sorgente](#2-mappa-del-codice-sorgente)
3. [Flusso dati end-to-end](#3-flusso-dati-end-to-end)
4. [BUG CRITICI (🔴) — da correggere prima di qualsiasi run](#4-bug-critici--da-correggere-prima-di-qualsiasi-run)
   - [BUG-1 · `import multiprocessing` fuori posto](#bug-1--import-multiprocessing-fuori-posto)
   - [BUG-2 · Nome file OTU sovrascritto](#bug-2--nome-file-otu-sovrascritto-con-stringa-hardcoded)
   - [BUG-3 · Directory dei risultati errata](#bug-3--directory-dei-risultati-errata--tabella-otu-sempre-vuota)
   - [BUG-4 · `base_name` non definita per file non-barcode](#bug-4--base_name-non-definita-per-file-non-barcode)
   - [BUG-5 · `process.wait()` senza timeout](#bug-5--processwait-senza-timeout--gui-si-blocca)
   - [BUG-6 · SQL DuckDB costruito con sintassi Python](#bug-6--sql-duckdb-costruito-con-sintassi-python-non-valida)
   - [BUG-7 · `fname` usata prima dell'assegnazione](#bug-7--fname-usata-prima-dellassegnazione-nellexception-handler)
5. [AVVERTENZE (⚠️) — da correggere prima della pubblicazione](#5-avvertenze--da-correggere-prima-della-pubblicazione)
   - [WARN-1 · Checkbox `force_tax_map` ignorata](#warn-1--checkbox-force_tax_map-ignorata-silenziosamente)
   - [WARN-2 · k-mer e window size non passati a minimap2](#warn-2--k-mer-e-window-size-non-passati-a-minimap2)
   - [WARN-3 · `shell=True` con path non quotati](#warn-3--shelltrue-con-path-non-quotati)
   - [WARN-4 · Divisione per zero nel report](#warn-4--divisione-per-zero-nel-report)
   - [WARN-5 · BAM intermedi non rimossi in caso di errore](#warn-5--bam-intermedi-non-rimossi-in-caso-di-errore)
   - [WARN-6 · Nomi colonne JS in italiano, Python in inglese](#warn-6--nomi-colonne-js-in-italiano-python-in-inglese--grafico-composizione-vuoto)
6. [Problemi di qualità del codice](#6-problemi-di-qualità-del-codice)
7. [Tabella di riepilogo](#7-tabella-di-riepilogo)
8. [Ordine di applicazione raccomandato](#8-ordine-di-applicazione-raccomandato)

---

## 1. Contesto dell'analisi

### Cosa fa LOREON

LOREON (Long-Read ONT Metagenomic Pipeline) è una **pipeline desktop** scritta in Python che elabora dati di sequenziamento Oxford Nanopore (ONT) per il metabarcoding fungino. Riceve in input le cartelle FASTQ demultiplexate per barcode, esegue filtraggio per lunghezza, allineamento contro un database di riferimento (UNITE, SILVA, CBS, Eukariome) e produce:

- `Final_OTU_Table.xlsx` — tabella OTU (taxa × campioni) con conteggi di letture
- `report_filtering_{min}_{max}.xlsx` — statistiche di qualità pre/post filtro
- `Report_{analisi}_{db}.html` — report interattivo con grafici Plotly

### Stack tecnologico

| Componente | Tecnologia |
|------------|------------|
| GUI | PyQt5 (QMainWindow + QThread) |
| Allineamento | `minimap2 -ax map-ont` |
| Elaborazione BAM | `samtools` (view, sort, index, idxstats) |
| Parsing sequenze | Biopython `SeqIO` |
| Aggregazione OTU | DuckDB SQL PIVOT |
| Report | Plotly + Jinja2 + Bootstrap 5 |
| Output tabelle | pandas + openpyxl |

### File analizzati (percorso base: `Applicazioni/loreon/`)

```
pipeline_gui.py        ← entry point, GUI PyQt5
pipeline_worker.py     ← bridge GUI↔CLI (QThread)
metaGenomics_new.py    ← orchestratore CLI (v2.6)
OtuUtils.py            ← libreria bioinformatica core
ResultsReader.py       ← aggregazione OTU con DuckDB
report_generator.py    ← generatore HTML (Plotly+Jinja2)
template.html          ← template Jinja2 con JS/Bootstrap
PipelineLogger.py      ← profiler DuckDB (non integrato)
fileorganizer.py       ← utility pre-processing .gz
requirements.txt       ← dipendenze Python
```

---

## 2. Mappa del codice sorgente

### Architettura dei processi

```
MainWindow (thread principale Qt)
    │
    │  start_pipeline_signal(dict settings)
    ▼
PipelineWorker (QThread)
    │
    │  subprocess.Popen(python metaGenomics_new.py ...)
    ▼
metaGenomics_new.py
    ├── Step 0: OtuUtils._create_tax_map()
    ├── Step 1: OtuUtils.filter_and_merge_directory() × N (multiprocessing.Pool)
    ├── Step 1.5: pandas → report_filtering.xlsx
    ├── Step 2: OtuUtils.mapping_improved() + tabeling_improved() × N (Pool)
    ├── Step 2.5: OtuUtils.reformat_tmp_indices()
    └── Step 3: ResultsReader.makeOtu_duckdb()
    │
    │  subprocess.Popen(python report_generator.py ...)
    ▼
report_generator.py → Report HTML
```

### Segnali Qt (thread-safe)

| Segnale | Direzione | Scopo |
|---------|-----------|-------|
| `start_pipeline_signal(dict)` | GUI → Worker | Avvia la pipeline con i parametri |
| `stop_pipeline_signal()` | GUI → Worker | Richiede interruzione |
| `log_signal(str)` | Worker → GUI | Ogni riga di stdout → log widget |
| `progress_signal(int, str)` | Worker → GUI | Aggiorna la progress bar |
| `finished_signal(bool, str)` | Worker → GUI | Successo/fallimento + path report |

### Aggiornamento progress bar (parsing sentinel da stdout)

| Stringa cercata in stdout | Progress % | Label |
|---------------------------|-----------|-------|
| `"--- Step 0:"` | 5% | Parsing Taxonomy |
| `"--- Step 1:"` | 10% | Filter |
| `"--- Step 1.5:"` | 40% | Filter Report |
| `"--- Step 2:"` | 50% | Mapping/Tabeling |
| `"--- Step 2.5:"` | 85% | Reformatting |
| `"--- Step 3:"` | 90% | OTU Aggregation |
| avvio report | 95% | HTML Report generation |

---

## 3. Flusso dati end-to-end

```
INPUT
  input_dir/
    barcode01/*.fastq.gz
    barcode02/*.fastq.gz
    unclassified/*.fastq.gz
  reference.fasta  (UNITE/SILVA/CBS/Eukariome)

  ▼ STEP 0 — _create_tax_map()
  Legge headers FASTA → estrae OTU_ID + stringa tassonomica
  Output: {db_stem}_taxonomy_map.tsv
  Formato: OTU_ID <TAB> Fungi;Ascomycota;...;Species

  ▼ STEP 1 — filter_and_merge_directory() [parallelo per barcode]
  Per ogni barcode: filtra per lunghezza, unisce tutti i file in uno
  Output: good_seq/barcode01.fastq   (letture OK)
          bad_seq/barcode01.fastq    (letture scartate)

  ▼ STEP 1.5 — pandas DataFrame
  Output: report_filtering_{min}_{max}.xlsx
  Colonne: barcode | mean_len | median_len | bad_count | good_count | total_count

  ▼ STEP 2 — mapping_improved() + tabeling_improved() [parallelo]
  minimap2 -ax map-ont → mapping/barcode01.sam
  samtools (8 step) → tabeling/barcode01_sorted_filtered.bam
                      tabeling/results/barcode01.tmp.txt
  Formato .tmp.txt: ref_id <TAB> ref_len <TAB> mapped <TAB> unmapped

  ▼ STEP 2.5 — reformat_tmp_indices()
  Riscrive colonna 1 dei .tmp.txt per allineare OTU_ID al formato del taxonomy map

  ▼ STEP 3 — makeOtu_duckdb()
  DuckDB PIVOT: tutti i .tmp.txt → matrice OTU_ID × barcode (conteggi)
  INNER JOIN con taxonomy_map.tsv
  Natural sort delle colonne (barcode1 < barcode2 < ... < unclassified)
  Output: Final_OTU_Table.xlsx
  Formato: OTU_ID (index) | Taxonomy | barcode01 | barcode02 | ...

  ▼ STEP 4 — report_generator.py
  Legge: report_filtering.xlsx + Final_OTU_Table.xlsx
  Genera: KPI, grafici Plotly (filtro, composizione tassonomica, heatmap)
  Output: Report_{analisi}_{db}.html (self-contained)

OUTPUT
  output_dir/
    mapping/barcode01.sam, ...
    tabeling/barcode01_sorted_filtered.bam, ...
    tabeling/results/barcode01.tmp.txt, ...
    good_seq/barcode01.fastq, ...
    bad_seq/barcode01.fastq, ...
    {db_stem}_taxonomy_map.tsv
    report_filtering_{min}_{max}.xlsx
    Final_OTU_Table.xlsx             ← BUG-2/BUG-3 impattano qui
    Report_{analisi}_{db}.html
```

---

## 4. Bug critici (🔴) — da correggere prima di qualsiasi run

---

### BUG-1 · `import multiprocessing` fuori posto

**File:** `pipeline_gui.py`
**Righe:** 101, 104 (uso) vs. 319 (import)
**Severità:** 🔴 Critico — crash al lancio su sistemi che importano il modulo come libreria

#### Causa

`multiprocessing.cpu_count()` viene chiamato dentro `init_ui()` (eseguita all'istanziazione di `MainWindow`), ma `import multiprocessing` è posizionato dentro il blocco `if __name__ == "__main__":` a riga 319.

Quando il file viene eseguito direttamente come script (`python pipeline_gui.py`), il blocco `__main__` viene eseguito *prima* che `MainWindow()` venga istanziata, quindi l'import avviene prima dell'uso — motivo per cui il bug non si manifesta nell'uso ordinario. Ma se il modulo venisse importato da un altro script, o se l'ordine di esecuzione cambiasse, si ottiene un `NameError` istantaneo.

#### Codice attuale

```python
# riga 101
self.total_threads_input.setRange(1, multiprocessing.cpu_count())  # ← usa multiprocessing
self.job_threads_input.setRange(1, multiprocessing.cpu_count())    # ← usa multiprocessing

# riga 317-319 — troppo tardi
if __name__ == "__main__":
    import multiprocessing          # ← import qui non è sicuro
    multiprocessing.freeze_support()
```

#### Fix proposto

```python
# --- INIZIO FILE (dopo gli altri import) ---
import sys
import multiprocessing              # ← spostare qui
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, ...
)
from PyQt5.QtCore import QThread, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices

# ...

if __name__ == "__main__":
    # import multiprocessing        ← eliminare questa riga
    multiprocessing.freeze_support()
    app = QApplication(sys.argv)
    window = MainWindow()
```

**Modifica:** spostare `import multiprocessing` dal blocco `__main__` (riga 319) alla sezione import in cima al file (dopo riga 11).

---

### BUG-2 · Nome file OTU sovrascritto con stringa hardcoded

**File:** `metaGenomics_new.py`
**Righe:** 157–161
**Severità:** 🔴 Critico — il worker cerca un file che non esiste mai con quel nome

#### Causa

Alla riga 157 viene costruito un nome file dinamico e descrittivo (`OTU_Table_{analisi}_{db}.xlsx`). Alla riga 158 viene assegnato a `final_output_file`. Alla riga 161, la stessa variabile viene **immediatamente sovrascritta** con `"Final_OTU_Table.xlsx"`.

Il file viene quindi salvato come `Final_OTU_Table.xlsx`, ma `pipeline_worker.py` (righe 84–85) costruisce il percorso `OTU_Table_{analisi}_{db}.xlsx` e lo passa a `report_generator.py`, che non trova il file e fallisce.

#### Codice attuale

```python
# riga 157-161
otu_filename = f"OTU_Table_{analysis_name}_{db_name}.xlsx"    # costruito correttamente
final_output_file = output_dir / otu_filename                  # assegnato
print(f"OTU Table output file name set on : {otu_filename}")   # stampato
final_output_file = output_dir / "Final_OTU_Table.xlsx"        # ← SOVRASCRIVE (BUG)
```

#### Fix proposto

```python
# riga 157-160 (riga 161 da eliminare)
otu_filename = f"OTU_Table_{analysis_name}_{db_name}.xlsx"
final_output_file = output_dir / otu_filename
print(f"OTU Table output file name set on : {otu_filename}")
# ← eliminare la riga 161
```

**Modifica:** eliminare la sola riga 161.

---

### BUG-3 · Directory dei risultati errata → tabella OTU sempre vuota

**File:** `metaGenomics_new.py`
**Riga:** 147
**Severità:** 🔴 Critico — causa più impattante: la tabella OTU è sempre vuota

#### Causa

`tabeling_improved()` in `OtuUtils.py` (riga 540) scrive i file `.tmp.txt` in:
```
output_dir / 'tabeling' / 'results' / barcode01.tmp.txt
```

Ma `metaGenomics_new.py` (riga 147) punta `makeOtu_duckdb()` a:
```
output_dir / 'results'   ← cartella separata, sempre vuota
```

`makeOtu_duckdb()` non trova nessun file `.tmp.txt`, la query DuckDB PIVOT restituisce zero righe, e la tabella OTU viene salvata vuota. Nessun messaggio di errore viene emesso.

#### Codice attuale

```python
# metaGenomics_new.py riga 145-150
mapping_dir = output_dir / 'mapping'
tabeling_dir = output_dir / 'tabeling'
results_dir = output_dir / 'results'       # ← percorso sbagliato

mapping_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)  # ← crea la cartella sbagliata
```

#### Referenza in OtuUtils.py (dove i file vengono realmente scritti)

```python
# OtuUtils.py riga 540-541 — tabeling_improved()
results_dir = base_dir / 'tabeling' / 'results'   # ← qui vengono scritti i .tmp.txt
tab_dir = base_dir / 'tabeling'
results_dir.mkdir(parents=True, exist_ok=True)
```

#### Fix proposto

```python
# metaGenomics_new.py riga 145-150
mapping_dir = output_dir / 'mapping'
tabeling_dir = output_dir / 'tabeling'
results_dir = output_dir / 'tabeling' / 'results'   # ← corretto

mapping_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)
```

**Modifica:** riga 147, cambiare `output_dir / 'results'` in `output_dir / 'tabeling' / 'results'`.

---

### BUG-4 · `base_name` non definita per file non-barcode

**File:** `OtuUtils.py`
**Riga:** 219 (`filterByLen_improved`)
**Severità:** 🔴 Critico — `NameError` su qualsiasi file non nominato `barcodeXX`

#### Causa

Nel ramo `else` (file che non corrispondono al pattern `barcodeXX`), viene assegnato `barcode_name = "unclassified"` ma `base_name` non viene mai definita. La riga 221 (`good_file_path = good_dir / base_name`) solleva `NameError`, che viene catturato dal blocco `except` e restituisce `None`, silenziando l'errore ma scartando tutti i file `unclassified`.

#### Codice attuale

```python
# OtuUtils.py riga 209-221
barcode_name = None
match = re.search(r'(barcode\d+)', input_file.name)

if match:
    barcode_name = match.group(1)
    base_name = barcode_name + '.fastq'   # ← assegnata solo nel ramo if
else:
    barcode_name = "unclassified"
    # base_name = ???                     ← MAI assegnata nel ramo else

good_file_path = good_dir / base_name    # ← NameError se si è nel ramo else
bad_file_path = bad_dir / base_name
```

#### Fix proposto

```python
if match:
    barcode_name = match.group(1)
    base_name = barcode_name + '.fastq'
else:
    barcode_name = "unclassified"
    base_name = "unclassified.fastq"     # ← aggiungere questa riga

good_file_path = good_dir / base_name
bad_file_path = bad_dir / base_name
```

**Modifica:** aggiungere `base_name = "unclassified.fastq"` alla riga 220 (dopo `barcode_name = "unclassified"`).

---

### BUG-5 · `process.wait()` senza timeout → GUI si blocca

**File:** `pipeline_worker.py`
**Riga:** 141
**Severità:** 🔴 Critico — la GUI può bloccarsi permanentemente dopo Stop

#### Causa

Quando l'utente clicca "Stop", `worker.stop()` chiama `process.terminate()` (SIGTERM). Il loop di lettura stdout esce (`if not self.is_running: break`), poi arriva alla riga 141: `self.process.wait()` — che blocca il thread Qt **senza timeout**. Se il sottoprocesso (es. `minimap2` durante un allineamento pesante) ignora SIGTERM, il worker thread non termina mai, rendendo la GUI non responsiva.

#### Codice attuale

```python
# pipeline_worker.py riga 133-142
with self.process.stdout:
    for line in iter(self.process.stdout.readline, ''):
        line = line.strip()
        if line:
            self.log_signal.emit(line)
            self.update_progress_from_log(line)
        if not self.is_running:
            break                           # ← esce dal loop correttamente
self.process.wait()                         # ← blocca senza timeout
local_return_code = self.process.returncode
```

#### Fix proposto

```python
with self.process.stdout:
    for line in iter(self.process.stdout.readline, ''):
        line = line.strip()
        if line:
            self.log_signal.emit(line)
            self.update_progress_from_log(line)
        if not self.is_running:
            break
try:
    self.process.wait(timeout=10)           # ← attendi max 10 secondi
except subprocess.TimeoutExpired:
    self.log_signal.emit("[WARN] Process did not terminate, sending SIGKILL...")
    self.process.kill()                     # ← forza la terminazione
    self.process.wait()
local_return_code = self.process.returncode
```

**Modifica:** sostituire `self.process.wait()` (riga 141) con il blocco try/except mostrato sopra.

---

### BUG-6 · SQL DuckDB costruito con sintassi Python non valida

**File:** `ResultsReader.py`
**Righe:** 347–358
**Severità:** 🔴 Critico — la query può fallire silenziosamente in alcune versioni di DuckDB

#### Causa

`str(col_types)` produce la rappresentazione Python di un dizionario:
```
{'otu_id': 'VARCHAR', 'length': 'BIGINT', 'mapped_count': 'BIGINT', 'unmapped_count': 'BIGINT'}
```
`str(clean_files)` produce la rappresentazione Python di una lista:
```
['/path/a/barcode01.tmp.txt', '/path/b/barcode02.tmp.txt']
```
Entrambe vengono interpolate direttamente nella stringa SQL. Alcune versioni di DuckDB le parsano in modo tollerante, ma è comportamento non documentato e version-dipendente.

#### Codice attuale

```python
# ResultsReader.py riga 342-366
col_types = {
    'otu_id': 'VARCHAR', 'length': 'BIGINT',
    'mapped_count': 'BIGINT', 'unmapped_count': 'BIGINT'
}

col_types_str = str(col_types)            # ← rappresentazione Python, non SQL
clean_files_str_list = str(clean_files)   # ← rappresentazione Python, non SQL

sql_query = f"""
WITH all_data AS (
    SELECT ...
    FROM read_csv(
        {clean_files_str_list}, delim='\t', header=False,
        columns={col_types_str}, filename=True
    )
    WHERE otu_id != '*'
)
PIVOT ...
"""
```

#### Fix proposto

```python
# Costruisce lista SQL valida: ['/path/a', '/path/b']
files_sql = "[" + ", ".join(f"'{f}'" for f in clean_files) + "]"

# Stringa literal SQL per i tipi colonna
col_types_sql = "{'otu_id': 'VARCHAR', 'length': 'BIGINT', 'mapped_count': 'BIGINT', 'unmapped_count': 'BIGINT'}"

sql_query = f"""
WITH all_data AS (
    SELECT
        replace(split_part(filename, '/', -1), '.tmp.txt', '') AS sample_name,
        otu_id, mapped_count
    FROM read_csv(
        {files_sql}, delim='\t', header=False,
        columns={col_types_sql}, filename=True
    )
    WHERE otu_id != '*'
)
PIVOT ( SELECT * FROM all_data )
ON sample_name
USING sum(mapped_count)
GROUP BY otu_id;
"""
```

**Modifica:** righe 347–358, sostituire `str(col_types)` e `str(clean_files)` con i literal SQL corretti.

---

### BUG-7 · `fname` usata prima dell'assegnazione nell'exception handler

**File:** `OtuUtils.py`
**Righe:** 528–557 (`tabeling_improved`)
**Severità:** 🔴 Critico — maschera l'eccezione originale con un secondo `NameError`

#### Causa

`fname = sam_path.stem` è assegnata alla riga 534, **dentro** il blocco `try`. Se un'eccezione avviene prima di quella riga (es. `Path(samfile)` fallisce per un path malformato), il blocco `except` alla riga 557 referenzia `fname` che non è ancora stata definita, sollevando un `NameError` che nasconde completamente l'errore originale.

#### Codice attuale

```python
# OtuUtils.py riga 528-558
try:
    sam_path = Path(samfile)       # se fallisce qui...
    fname = sam_path.stem          # ← fname non è ancora definita
    base_dir = sam_path.parent.parent
    ...
except Exception as e:
    print(f"[{fname}] ERRORE in setup percorsi: {e}")  # ← NameError! fname non esiste
    return
```

#### Fix proposto

```python
fname = Path(samfile).stem if samfile else "unknown"   # ← inizializzazione sicura FUORI dal try

try:
    sam_path = Path(samfile)
    fname = sam_path.stem          # continua ad aggiornare normalmente
    base_dir = sam_path.parent.parent
    ...
except Exception as e:
    print(f"[{fname}] ERRORE in setup percorsi: {e}")  # ← ora fname è sempre definita
    return
```

**Modifica:** aggiungere `fname = Path(samfile).stem if samfile else "unknown"` alla riga 527 (prima del `try`).

---

## 5. Avvertenze (⚠️) — da correggere prima della pubblicazione

---

### WARN-1 · Checkbox `force_tax_map` ignorata silenziosamente

**File:** `pipeline_gui.py`
**Riga:** 221–236 (`run_pipeline`)

La checkbox `force_tax_map_check` è presente nella GUI (riga 289 usa `set_running_state`), ma il suo valore **non viene mai aggiunto al dizionario `settings`**. Il worker legge sempre `False` per questa chiave.

#### Codice attuale

```python
settings = {
    "input_dir": self.input_dir_label.text(),
    "output_dir": self.output_dir_label.text(),
    "db_path": self.db_label.text(),
    "format": self.format_combo.currentText(),
    "debug": self.debug_check.isChecked(),
    "enable_filter": self.filter_check.isChecked(),
    "min_len": self.min_len_input.value(),
    "max_len": self.max_len_input.value(),
    "total_threads": self.total_threads_input.value(),
    "threads_per_job": self.job_threads_input.value(),
    "kmer_size": self.kmer_input.value(),
    "window_size": self.window_input.value(),
    "open_report": self.open_report_check.isChecked()
    # "force_tax_map" ← MANCANTE
}
```

#### Fix proposto

```python
settings = {
    ...
    "open_report": self.open_report_check.isChecked(),
    "force_tax_map": self.force_tax_map_check.isChecked()   # ← aggiungere
}
```

---

### WARN-2 · k-mer e window size non passati a minimap2

**File:** `metaGenomics_new.py` (righe 302–308) + `OtuUtils.py` (riga 410)

I parametri `-k` e `-w` vengono esposti nella GUI e parsati come argomenti CLI, ma non vengono mai inclusi in `map_worker_kwargs` né nella signature di `mapping_improved()`.

#### Codice attuale

```python
# metaGenomics_new.py riga 302-308
map_worker_kwargs = {
    'db_path': str(db_path),
    'threads': threads_per_job_map,
    'output_dir': str(mapping_dir),
    'debug': options.debug
    # kmer_size e window_size: assenti
}

# OtuUtils.py riga 410 — firma della funzione
def mapping_improved(fastq_file: str, db_path: str, threads: int,
                     output_dir: str, debug: bool = False):
    # kmer_size e window_size: non accettati

# OtuUtils.py riga 448-454 — comando minimap2
cmd_list = [
    'minimap2', '-ax', 'map-ont',
    '-t', thread_str,
    str(db_path_obj), str(fastq_path)
    # -k e -w: mai aggiunti
]
```

#### Fix proposto (in tre punti)

```python
# 1. metaGenomics_new.py riga 302-308
map_worker_kwargs = {
    'db_path': str(db_path),
    'threads': threads_per_job_map,
    'output_dir': str(mapping_dir),
    'debug': options.debug,
    'kmer_size': options.kmer_size,       # ← aggiungere
    'window_size': options.window_size    # ← aggiungere
}

# 2. OtuUtils.py riga 410 — aggiungere parametri alla firma
def mapping_improved(fastq_file: str, db_path: str, threads: int,
                     output_dir: str, debug: bool = False,
                     kmer_size: int = None, window_size: int = None):  # ← aggiungere

# 3. OtuUtils.py riga 448-454 — aggiungere -k e -w al comando
cmd_list = ['minimap2', '-ax', 'map-ont', '-t', thread_str]
if kmer_size:
    cmd_list.extend(['-k', str(kmer_size)])     # ← aggiungere
if window_size:
    cmd_list.extend(['-w', str(window_size)])   # ← aggiungere
cmd_list.extend([str(db_path_obj), str(fastq_path)])
```

---

### WARN-3 · `shell=True` con path non quotati

**File:** `OtuUtils.py`
**Righe:** 566–595 (`tabeling_improved`)

Tutti i comandi `samtools` sono costruiti con f-string e passati a `subprocess.run(..., shell=True)`. Path contenenti spazi o metacaratteri shell (comuni su macOS con nomi di cartelle con spazi) rompono il comando silenziosamente.

#### Esempio del problema

```python
# OtuUtils.py riga 566
cmd_unfiltered = (
    f"samtools view -@{thread_str} -bS {samfile} | "   # ← samfile non quotato
    f"samtools sort -@{thread_str} -o {unfiltered_sorted_bam} -"
)
subprocess.run(cmd_unfiltered, shell=True, ...)
```

#### Fix proposto (approccio consigliato)

```python
import shlex

# Quotare ogni path con shlex.quote()
cmd_unfiltered = (
    f"samtools view -@{thread_str} -bS {shlex.quote(samfile)} | "
    f"samtools sort -@{thread_str} -o {shlex.quote(str(unfiltered_sorted_bam))} -"
)
```

---

### WARN-4 · Divisione per zero nel report

**File:** `report_generator.py`
**Riga:** 66

`total_seq` è la somma della colonna `'total count'`. Se il filtro scarta tutte le letture (o il file di filtro è vuoto), `total_seq == 0` e la divisione solleva `ZeroDivisionError`, facendo crashare il generatore del report.

#### Codice attuale

```python
# report_generator.py riga 63-67
total_seq = df_filter['total count'].sum()
good_seq = df_filter['good count'].sum()
kpis = {
    ...
    'rejection_rate': f"{(total_seq - good_seq) / total_seq:.1%}",  # ← ZeroDivisionError
    ...
}
```

#### Fix proposto

```python
'rejection_rate': (
    f"{(total_seq - good_seq) / total_seq:.1%}"
    if total_seq > 0 else "N/A"
),
```

---

### WARN-5 · BAM intermedi non rimossi in caso di errore

**File:** `OtuUtils.py`
**Righe:** 620–640 (`tabeling_improved`)

I file BAM intermedi (`.sorted.bam`, `.sorted.bam.bai`) vengono rimossi alla fine del blocco `try`. Se si verifica un'eccezione a metà elaborazione, i file temporanei rimangono su disco, potenzialmente occupando decine di GB.

#### Fix proposto

```python
# Aggiungere un blocco finally per il cleanup
try:
    # ... elaborazione BAM ...
except Exception as e:
    print(f"[{fname}] ERRORE: {e}")
finally:
    # Cleanup garantito anche in caso di errore
    for tmp_file in [unfiltered_sorted_bam,
                     Path(str(unfiltered_sorted_bam) + ".bai")]:
        if tmp_file.exists():
            tmp_file.unlink()
```

---

### WARN-6 · Nomi colonne JS in italiano, Python in inglese → grafico composizione vuoto

**File:** `template.html` (righe 140–146) vs `report_generator.py` (riga 168)
**Severità:** ⚠️ Alta — il grafico interattivo principale del report è sempre vuoto

#### Causa

`report_generator.py` `aggregate_by_rank()` produce un DataFrame long-format con colonne:
```python
# report_generator.py riga 163-168
df_long = df_rel.reset_index().melt(
    id_vars=rank_name,
    var_name='Sample',                   # ← inglese
    value_name='Relative Abundance (%)'  # ← inglese
)
df_long = df_long.rename(columns={rank_name: 'Taxonomy'})  # ← inglese
```

Il template JavaScript accede con nomi **italiani**:
```javascript
// template.html riga 140-146
const groups = [...new Set(filteredData.map(d => d.Tassonomia))];       // ← italiano
const groupData = filteredData.filter(d => d.Tassonomia === groupName);  // ← italiano
return {
    x: groupData.map(d => d.Campione),                    // ← italiano
    y: groupData.map(d => d['Abbondanza Relativa (%)']),  // ← italiano
```

`d.Tassonomia` è `undefined` → `groups` è un array vuoto → nessun trace → grafico vuoto.

#### Fix proposto (modificare il JavaScript per usare i nomi inglesi)

```javascript
// template.html riga 140-146 — DOPO la correzione
const groups = [...new Set(filteredData.map(d => d.Taxonomy))];
const groupData = filteredData.filter(d => d.Taxonomy === groupName);
return {
    x: groupData.map(d => d.Sample),
    y: groupData.map(d => d['Relative Abundance (%)']),
    name: groupName,
    type: 'bar'
};
```

**Alternativa:** rinominare le colonne in `report_generator.py` in italiano (ma l'approccio in inglese è preferibile per coerenza con il resto del codice).

---

## 6. Problemi di qualità del codice

| Problema | File | Dettaglio |
|----------|------|-----------|
| Codice duplicato | `OtuUtils.py` | `_run_fasta_parsing()` e `_create_tax_map()` hanno logica quasi identica |
| Codice morto | `OtuUtils.py` | 3 generazioni di funzioni filter/tabeling, solo l'ultima è usata |
| Classe inutilizzabile | `ResultsReader.py` | La classe `ResultsReader` importa `mappedBarcode` (non presente nel repo) |
| Colonne invertite | `ResultsReader.py` righe 62–63 | I valori di Family e Order sono assegnati alle colonne sbagliate |
| Typo nei nomi colonna | `ResultsReader.py` | `'kigdom'` → `'kingdom'`, `'philum'` → `'phylum'` |
| Lingua mista | Tutti i file | Commenti, print, nomi variabile mescolano italiano e inglese |
| Attributi di classe mutabili | `ResultsReader.py` righe 13–21 | `files_dict`, `barcodes`, `log` sono condivisi tra istanze |
| DataFrame mutato in-place | `report_generator.py` riga 136 | `aggregate_by_rank()` modifica il DataFrame del chiamante |
| `PipelineLogger` non integrato | — | Esiste ma non viene mai chiamato dalla pipeline attiva |

---

## 7. Tabella di riepilogo

| ID | Tipo | File | Riga/i | Descrizione | Impatto runtime |
|----|------|------|--------|-------------|-----------------|
| BUG-1 | 🔴 Critico | `pipeline_gui.py` | 101, 319 | `import multiprocessing` nel posto sbagliato | Crash potenziale al lancio |
| BUG-2 | 🔴 Critico | `metaGenomics_new.py` | 161 | Nome file OTU sovrascritto con hardcoded | Report non trova OTU table |
| BUG-3 | 🔴 Critico | `metaGenomics_new.py` | 147 | Directory risultati sbagliata | OTU table sempre vuota |
| BUG-4 | 🔴 Critico | `OtuUtils.py` | 219 | `base_name` non definita per file non-barcode | Crash silenzioso su unclassified |
| BUG-5 | 🔴 Critico | `pipeline_worker.py` | 141 | `process.wait()` senza timeout | GUI si blocca su Stop |
| BUG-6 | 🔴 Critico | `ResultsReader.py` | 347–358 | SQL DuckDB con sintassi Python | Query invalida in alcune versioni |
| BUG-7 | 🔴 Critico | `OtuUtils.py` | 528–557 | `fname` non inizializzata pre-try | NameError maschera eccezione reale |
| WARN-1 | ⚠️ Alta | `pipeline_gui.py` | 235 | `force_tax_map` mancante nel dict settings | Checkbox senza effetto |
| WARN-2 | ⚠️ Alta | `metaGenomics_new.py` + `OtuUtils.py` | 303–308, 410, 448 | k-mer e window size non passati | Parametri GUI ignorati |
| WARN-3 | ⚠️ Media | `OtuUtils.py` | 566–595 | `shell=True` con path non quotati | Crash su path con spazi |
| WARN-4 | ⚠️ Media | `report_generator.py` | 66 | Divisione per zero | Crash se 0 reads passano il filtro |
| WARN-5 | ⚠️ Bassa | `OtuUtils.py` | 620–640 | Cleanup BAM solo in try, non in finally | File temporanei su disco |
| WARN-6 | ⚠️ Alta | `template.html` | 140–146 | Nomi colonne JS italiani ≠ Python inglese | Grafico composizione sempre vuoto |

---

## 8. Ordine di applicazione raccomandato

### Fase 1 — Correzioni minimali per far funzionare la pipeline (< 30 min)

```
1. metaGenomics_new.py riga 147:  'results' → 'tabeling/results'         (BUG-3)
2. metaGenomics_new.py riga 161:  eliminare la riga di overwrite           (BUG-2)
3. template.html     righe 140-146: nomi JS italiano → inglese             (WARN-6)
```
> Queste tre correzioni insieme sono sufficienti per ottenere una prima run funzionante con output corretto.

### Fase 2 — Robustezza e correttezza (< 1 ora)

```
4. pipeline_gui.py   riga 319:   spostare import multiprocessing in cima   (BUG-1)
5. OtuUtils.py       riga 219:   aggiungere base_name = "unclassified.fastq" (BUG-4)
6. pipeline_worker.py riga 141:  aggiungere timeout a process.wait()       (BUG-5)
7. ResultsReader.py  riga 347:   SQL DuckDB con literal corretti            (BUG-6)
8. OtuUtils.py       riga 527:   inizializzare fname prima del try          (BUG-7)
```

### Fase 3 — Funzionalità mancanti (< 2 ore)

```
9.  pipeline_gui.py  riga 235:   aggiungere force_tax_map al dict settings  (WARN-1)
10. metaGenomics_new.py + OtuUtils.py: passare kmer_size/window_size        (WARN-2)
11. report_generator.py riga 66: guard divisione per zero                    (WARN-4)
12. OtuUtils.py tabeling_improved: aggiungere finally per cleanup BAM        (WARN-5)
```

### Fase 4 — Qualità e manutenibilità (refactoring, non urgente)

```
13. OtuUtils.py:      unificare _run_fasta_parsing e _create_tax_map
14. OtuUtils.py:      rimuovere funzioni legacy (filterByLenOld, tabeling vecchio)
15. ResultsReader.py: rimuovere classe ResultsReader (dipendenza mappedBarcode mancante)
16. ResultsReader.py: correggere typo 'kigdom'/'philum' e colonne Family/Order invertite
17. OtuUtils.py:      aggiungere shlex.quote() ai comandi samtools (WARN-3)
18. Tutti i file:     standardizzare la lingua (solo inglese)
19. Integrare PipelineLogger nella pipeline attiva per il profiling
```

---

*Fine documento — generato da revisione statica del codice sorgente in data 2026-02-26*
*Nessuna modifica è stata applicata al codice. Tutte le proposte sono da validare prima dell'applicazione.*
