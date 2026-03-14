# LOREON Pipeline — Optimization Report

**Data analisi:** 2026-03-13

---

## Priorità CRITICA

### OPT-1: Parser FASTQ raw al posto di BioPython SeqIO nel filtro
- **File:** `OtuUtils.py` — `filter_and_merge_directory`
- **Problema:** `SeqIO.parse` costruisce un oggetto `SeqRecord` per ogni read. Per il filtro per lunghezza servono solo 4 righe raw e `len(seq)`.
- **Soluzione:** Parser FASTQ nativo a 4 righe, skip completo di BioPython.
- **Speedup stimato:** 5-10x sul throughput del filtro.

### OPT-2: Riduzione chiamate samtools da 9 a 4 per sample
- **File:** `OtuUtils.py` — `tabeling_improved`
- **Problema:** 9 subprocess sequenziali per sample con ridondanze: `idxstats` chiamato 2 volte sullo stesso BAM, BAM unfiltered creato solo per ottenere conteggi ricavabili da `flagstat`.
- **Soluzione:** Pipeline diretta `minimap2 | samtools view -e | samtools sort`, eliminare subprocess ridondanti.
- **Speedup stimato:** 2-4x sulla fase di tabeling.

### OPT-3: `xlsxwriter` al posto di `openpyxl` per scrittura OTU
- **File:** `ResultsReader.py`
- **Problema:** `openpyxl` costruisce albero XML in memoria per ogni cella. Per 50k OTU × 24 campioni impiega 30-120s.
- **Soluzione:** `engine='xlsxwriter'` (drop-in replacement).
- **Speedup stimato:** 5-15x sulla scrittura della tabella OTU.

### OPT-4: Parser raw per header FASTA + unificazione funzioni duplicate
- **File:** `OtuUtils.py` — `_create_tax_map` / `_run_fasta_parsing`
- **Problema:** `SeqIO.parse` legge header + sequenza intera, ma servono solo gli header. Le due funzioni sono codice duplicato.
- **Soluzione:** `for line in f: if line.startswith('>')` + funzione unica.
- **Speedup stimato:** 3-5x sulla creazione della taxonomy map.

---

## Priorità ALTA

### OPT-5: Join tassonomico dentro DuckDB
- **File:** `ResultsReader.py`
- **Problema:** Carica intero TSV tassonomico (500k+ righe SILVA) in pandas, poi inner join. La maggior parte delle righe viene scartata.
- **Soluzione:** Spostare il join nella query DuckDB prima di materializzare il DataFrame.
- **Impatto:** Riduzione significativa del picco di memoria.

### OPT-6: Sovrapposizione filtering ↔ mapping (producer-consumer)
- **File:** `metaGenomics_new.py`
- **Problema:** Tutti i barcode finiscono il filtro prima che inizi il mapping. Con `imap_unordered` + callback i due step si sovrappongono.
- **Speedup stimato:** ~20-30s di wall time in meno.

### OPT-7: Formati intermedi XLSX → TSV
- **File:** `metaGenomics_new.py`, `report_generator.py`
- **Problema:** Il filtering report fa round-trip XLSX (write openpyxl → read openpyxl) per un DataFrame da ~24 righe.
- **Soluzione:** Usare `.tsv` come formato intermedio.
- **Speedup stimato:** ~100x su I/O intermedio.

---

## Priorità MEDIA

### OPT-8: Reservoir sampler per la mediana
- **File:** `OtuUtils.py` — `filter_and_merge_directory`
- **Problema:** `all_lengths` accumula tutti i length in RAM (1M reads/barcode × 24 worker = ~192 MB).
- **Soluzione:** Reservoir sampler da 10.000 elementi (~80 KB/worker).

### OPT-9: Vettorizzare `aggregate_by_rank`
- **File:** `report_generator.py`
- **Problema:** `.apply(lambda x: ...)` esegue 50k chiamate Python per rank (×3 = 150k).
- **Soluzione:** `.str.split(';', expand=True)` una volta sola per tutti i rank.
- **Speedup stimato:** 3-5x sulla generazione report.

### OPT-10: Limitare tabella OTU in HTML a top-N righe
- **File:** `report_generator.py`
- **Problema:** `to_html()` con 50k righe genera centinaia di MB di HTML.
- **Soluzione:** Mostrare top-500 OTU, nota che la tabella completa è nel file XLSX.

### OPT-11: Parallelizzare `reformat_tmp_indices`
- **File:** `OtuUtils.py` — Step 2.5 / `reformat_tmp_indices`
- **Problema:** File `.tmp.txt` indipendenti riformattati in sequenza.
- **Soluzione:** Pool parallelo come per filter e mapping.

---

## Priorità BASSA

### OPT-12: Cache backend memoria nel profiler
- **File:** `pipeline_profiler.py` — `_get_mem_mb`
- **Problema:** Try/import/except chain eseguita ad ogni `end_step`.
- **Soluzione:** Detect backend una volta in `__init__`, salvare callable.

### OPT-13: Rimuovere auto-scroll ridondante nel log GUI
- **File:** `pipeline_gui.py` — `append_log`
- **Problema:** `scrollBar.setValue(maximum())` ridondante, `appendPlainText` fa già auto-scroll.
- **Soluzione:** Rimuovere la chiamata esplicita a `setValue`.

### OPT-14: Opzione GUI "Delete BAM" per pulizia file intermedi
- **File:** `pipeline_gui.py`, `pipeline_worker.py`, `OtuUtils.py`
- **Problema:** SAM e BAM filtrati restano su disco (multi-GB per 24 barcode).
- **Soluzione:** Checkbox in GUI → flag passato al worker → pulizia a fine tabeling.

### OPT-15: Timeout `process.wait()` configurabile
- **File:** `pipeline_worker.py`
- **Problema:** 600s hardcoded, insufficiente per dataset grandi.
- **Soluzione:** Rimuovere timeout (già protetto dal loop stdout).

---

## Rinominazione Step

| Vecchio nome | Nuovo nome |
|---|---|
| Step 0: Taxonomy | Step 1: Taxonomy |
| Step 1: Filter | Step 2: Filter |
| Step 1.5: Filter Report | Step 3: Filter Report |
| Step 2: Mapping | Step 4: Mapping |
| Step 2.5: Reformat | Step 5: Reformat |
| Step 3: OTU Aggregation | Step 6: OTU Aggregation |
| Step 4: HTML Report (GUI) | Step 7: HTML Report (GUI) |
