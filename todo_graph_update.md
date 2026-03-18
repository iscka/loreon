# TODO: Graph Update — Integrazione Grafici EPI2ME in LOREON

## Analisi dei Report EPI2ME di Riferimento

### Report 1: `wf-16s-report.html` (Report Analisi 16S)

| # | Sezione | Tipo Grafico | Libreria | Dati Visualizzati |
|---|---------|-------------|----------|-------------------|
| 1 | **Read Summary** | 3 Bar Chart | Bokeh (JSON embedded) | Statistiche qualità reads per barcode (lunghezza, qualità, conteggi) |
| 2 | **Number of Reads** | Tabella HTML | Bootstrap | Conteggi reads per campione (totali, non classificati, %) |
| 3 | **Lineages** | Sankey Diagram interattivo | D3.js custom | Flusso tassonomico gerarchico (superkingdom → genus) con filtri abbondanza |
| 4 | Sunburst | Sunburst Chart | ECharts | Composizione tassonomica gerarchica |
| 5 | Taxonomy | 5 Bar Chart (tabs) | Bokeh | Top 9 taxa per rango (phylum→genus), abbondanza relativa |
| 6 | Abundances | 10 Tabelle (tabs) | SimpleDataTables | Conteggi per rango + versione rarefatta |
| 7 | Alpha Diversity | Tabella + 2 Bokeh | Bokeh + HTML | Indici diversità (Shannon, Simpson, ecc.) + curve rarefazione |

### Report 2: `execution/report.html` (Report Esecuzione Nextflow)

| # | Sezione | Tipo Grafico | Libreria | Dati Visualizzati |
|---|---------|-------------|----------|-------------------|
| 1 | **CPU** | Box Plot (2 tab: raw + % allocated) | Plotly.js 1.34 | Uso CPU per processo (% single core + % allocato) |
| 2 | **Memory** | Box Plot (3 tab: RAM + virtual + % allocated) | Plotly.js 1.34 | Memoria fisica, virtuale, % allocata per processo |
| 3 | **Job Duration** | Box Plot (2 tab: raw + % allocated) | Plotly.js 1.34 | Durata esecuzione per processo (minuti + % allocato) |
| 4 | **I/O** | Box Plot (2 tab: read + write) | Plotly.js 1.34 | Bytes letti/scritti per processo |
| 5 | Tasks Table | DataTable interattiva | jQuery DataTables | Dettaglio singoli task con 30+ campi |

---

## Piano di Integrazione

### PARTE A: Primi 3 grafici wf-16s → Report Analisi Finale (`report_generator.py` + `template.html`)

#### A1. Read Summary — 3 Bar Chart ✅ FATTIBILE (complessità: MEDIA)

**Stato attuale LOREON:** Il report ha già grafici simili ("Filter Count Report" con good/bad counts e "Mean Sequence Length"). La sezione Read Summary di EPI2ME mostra 3 grafici Bokeh con statistiche qualità per barcode.

**Azione proposta:** Riadattare/arricchire i grafici di filtro già esistenti in LOREON per includere:
- Distribuzione lunghezza reads (istogramma) per campione
- Qualità media reads per campione
- Conteggi reads (buoni/scartati/totali) — *già presente*

**Implementazione:**
- Libreria: **Plotly.js** (già in uso, nessuna nuova dipendenza)
- Dati: Già disponibili in `df_filter` (colonne: `barcode`, `total count`, `good count`, `bad count`, `mean length`, `median length`)
- Manca: distribuzione qualità reads (EPI2ME usa dati da FastQ quality scores non disponibili nel pipeline LOREON attuale)
- File da modificare: `report_generator.py` (funzione `generate_filter_graphs`), `template.html`

**Limitazioni:**
- I 3 grafici Bokeh di EPI2ME usano dati di qualità base-per-base che LOREON non raccoglie
- Si può replicare il concetto con i dati disponibili (conteggi + lunghezze)

---

#### A2. Number of Reads — Tabella ✅ FATTIBILE (complessità: BASSA)

**Stato attuale LOREON:** Esiste già la tabella `table_filter_html` con dati filtro per barcode.

**Azione proposta:** Aggiungere una riga riepilogativa con:
- Totale reads
- Reads non classificati (unmapped) — *disponibile se `mapping_stats.json` è presente*
- Percentuale non classificati

**Implementazione:**
- Dati: combinazione di `df_filter` + `mapping_data` (già caricati)
- File da modificare: `report_generator.py` (funzione `get_kpis` o nuovo riepilogo), `template.html`
- Nessuna nuova dipendenza

---

#### A3. Lineages — Sankey Diagram ⚠️ FATTIBILE (complessità: ALTA)

**Stato attuale LOREON:** Il report ha grafici di composizione tassonomica (stacked bar per Family/Genus/Species) e heatmap. Non ha un Sankey.

**Azione proposta:** Aggiungere un diagramma Sankey interattivo che mostri il flusso tassonomico gerarchico dai livelli superiori (Superkingdom) fino a Genus/Species.

**Implementazione — Opzione 1 (D3.js, come EPI2ME):**
- Pro: Risultato identico a EPI2ME, massima interattività (filtri abbondanza, ricerca, download SVG)
- Contro: Aggiunge D3.js come dipendenza (~250KB CDN), codice JavaScript complesso (~300 righe custom), data transformation significativa
- Dati necessari: La colonna `Taxonomy` dell'OTU table contiene già la gerarchia completa (`Superkingdom;Phylum;Class;Order;Family;Genus;Species`)

**Implementazione — Opzione 2 (Plotly Sankey, consigliata):**
- Pro: Usa la stessa libreria già presente (Plotly), API Sankey nativa (`go.Sankey`), meno codice
- Contro: Meno personalizzabile della versione D3 di EPI2ME, ma comunque interattivo
- Codice stimato: ~80 righe Python in `report_generator.py`
- Dati: Parsare `Taxonomy` per costruire nodi e link (edges) tra livelli tassonomici

**Struttura dati richiesta (da OTU table):**
```python
# Da: "Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus"
# A:  nodes = ["Bacteria", "Bacillota", "Bacilli", ...]
#     links = [{"source": 0, "target": 1, "value": 282054}, ...]
```

**File da modificare:** `report_generator.py` (nuova funzione `generate_sankey`), `template.html` (nuova sezione)

---

### PARTE B: Grafici execution/report → Performance Report (`pipeline_profiler.py`)

#### B1. CPU Box Plot ✅ FATTIBILE (complessità: MEDIA)

**Stato attuale LOREON:** Il `performance_report.html` mostra CPU time come valore singolo per step (bar chart). Non ha distribuzione statistica.

**Problema:** LOREON non è un workflow Nextflow — non ha processi paralleli con varianza CPU. Ogni step ha UN solo valore `cpu_time_s`. Un box plot non avrebbe senso con un singolo data point per step.

**Azione proposta — Alternativa adattata:**
- Se lo step esegue task paralleli (es. Step 4 Mapping con N file), raccogliere CPU time *per task* e mostrare box plot
- Per step sequenziali: mantenere il bar chart attuale
- Aggiungere grafico **CPU efficiency** (rapporto `cpu_time / wall_time` per step) come bar chart
- File da modificare: `pipeline_profiler.py` (sezione `_save_html`), `metaGenomics_new.py` (raccogliere metriche per-task)

**Prerequisito:** Modificare `PipelineProfiler` per accettare metriche aggregate (lista di valori) oltre ai singoli valori.

---

#### B2. Memory Box Plot ✅ FATTIBILE (complessità: MEDIA)

**Stato attuale LOREON:** Bar chart con RSS per step + linea peak. Singolo valore per step.

**Stesso problema del CPU:** Un box plot necessita di distribuzioni. Con singoli valori, non ha senso.

**Azione proposta — Alternativa adattata:**
- Grafico **Memory progression** (line chart): mostrare l'andamento RSS nel tempo (inizio vs fine di ogni step)
- Aggiungere un secondo asse con `wall_time` per vedere la correlazione tempo-memoria
- Richiederebbe campionamento memoria a intervalli regolari durante l'esecuzione (feature futura)

---

#### B3. Job Duration Box Plot ✅ FATTIBILE (complessità: BASSA)

**Stato attuale LOREON:** Già presente come bar chart + Gantt timeline. Funzionalmente equivalente.

**Azione proposta:**
- Mantenere i grafici attuali (bar chart + Gantt) che sono più informativi per un pipeline sequenziale
- Aggiungere **waterfall chart** (variazione incrementale): mostra quanto ogni step aggiunge al tempo totale
- Implementazione Plotly: `type: 'waterfall'` — ~30 righe JS

---

#### B4. I/O Chart ⚠️ FATTIBILE (complessità: ALTA)

**Stato attuale LOREON:** Non presente.

**Problema:** LOREON non raccoglie metriche I/O per step. Necessario:
1. Aggiungere tracking I/O in `PipelineProfiler` (via `resource.getrusage` su Linux/macOS)
2. `resource.RUSAGE_SELF` fornisce `ru_inblock` e `ru_oublock` (blocchi I/O)
3. Alternativa: usare `psutil.Process().io_counters()` (cross-platform)

**Azione proposta:**
- Fase 1: Aggiungere raccolta dati I/O in `PipelineProfiler.end_step()`
- Fase 2: Aggiungere grafico bar I/O nel report HTML
- Dipendenza opzionale: `psutil` (già tentato come fallback in `_detect_mem_backend`)

---

## Riepilogo Priorità e Roadmap

| ID | Grafico | Target Report | Complessità | Priorità | Dipendenze Nuove |
|----|---------|--------------|-------------|----------|-----------------|
| A1 | Read Summary (bar charts arricchiti) | Analisi | MEDIA | 🟢 ALTA | Nessuna |
| A2 | Number of Reads (tabella riepilogativa) | Analisi | BASSA | 🟢 ALTA | Nessuna |
| A3 | Lineages Sankey | Analisi | ALTA | 🟡 MEDIA | Nessuna (Plotly Sankey) |
| B1 | CPU Efficiency chart | Performance | MEDIA | 🟡 MEDIA | Modifiche al profiler |
| B2 | Memory Progression | Performance | MEDIA | 🟡 MEDIA | Campionamento memoria |
| B3 | Waterfall Duration | Performance | BASSA | 🟢 ALTA | Nessuna |
| B4 | I/O Chart | Performance | ALTA | 🔴 BASSA | `psutil` o `resource` |

### Ordine di Implementazione Consigliato

```
Sprint 1 (Rapido):  A2 → B3 → A1
Sprint 2 (Medio):   A3 (Sankey Plotly) → B1 (CPU efficiency)
Sprint 3 (Futuro):  B2 (memory sampling) → B4 (I/O tracking)
```

---

## Note Tecniche

### Compatibilità Librerie
- **LOREON usa già:** Plotly.js (latest), Bootstrap 5.3, DataTables, Jinja2
- **EPI2ME usa:** Bokeh (embedded), D3.js, ECharts, Bootstrap 4, Plotly 1.34
- **Strategia:** Reimplementare tutti i grafici in Plotly.js per uniformità (nessuna nuova dipendenza frontend)

### Formato Dati Disponibili in LOREON
- `df_filter`: barcode, total count, good count, bad count, mean length, median length
- `df_otu`: OTU_ID (index), Taxonomy, sample columns con conteggi
- `mapping_data`: total_reads, mapped_unfiltered, mapped_filtered, mapping_rate per campione
- `PipelineProfiler._steps`: wall_time_s, cpu_time_s, mem_mb + metriche custom per step

### Vincolo Principale
LOREON è un pipeline Python sequenziale, non un workflow Nextflow distribuito. I grafici statistici (box plot) di Nextflow presuppongono molteplici esecuzioni dello stesso processo — non applicabile direttamente. La strategia di adattamento usa grafici alternativi più informativi per il contesto LOREON.
