# LOREON example dataset

A small dataset to check that an installation works and to see what LOREON produces.

## Content

| Path | Description |
|------|-------------|
| `reads/MockA.fastq.gz`, `MockB.fastq.gz`, `MockC.fastq.gz` | 500 reads from each of three yeast mock communities. The amplicon spans the ITS region and the D1/D2 domains of the LSU (about 1.0–1.6 kb) and was sequenced on an Oxford Nanopore MinION. One file is one sample. |
| `database/unite_subset.fasta` | 1,041 sequences of the UNITE general FASTA release for Fungi: the 41 sequences hit by the example reads when they are mapped on the whole release, plus 1,000 sequences drawn at random with a fixed seed. Identifiers, taxonomy and sequences are unchanged. |
| `expected/OTU_Table_expected.tsv` | The OTU table that LOREON writes for these files with the default settings. |

## Running it

Command line:

```bash
python metaGenomics_new.py --example
```

Graphical interface: **RUN EXAMPLE DATASET**.

The run takes less than a minute. It writes its results to `./loreon_example_output` (command line; change it with `-o`) or to `LOREON_example_output` in the home folder (interface), compares the read counts of the OTU table with the expected table and reports whether they agree. On the command line the last line is `The installation works.` and the exit status is `0`; if the tables differ, the differing cells are listed and the exit status is `1`.

The expected table was produced with minimap2 2.28 and samtools 1.21; the same counts were obtained with minimap2 2.24 and samtools 1.13, on macOS and Linux, on x86-64 and ARM processors.

The counts are not those of an analysis on the whole UNITE release, because the database is a subset.

The same comparison is available for any analysis with `--expected-table <table.tsv>`.

## Sources and licences

**Reads.** NCBI BioProject [PRJNA862129](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA862129), *MinION rDNA sequences of yeast mock communities*, University of Perugia. The files hold the first 500 reads of mock communities A, B and C (biological replicate 2, technical replicate 1); the read headers were reduced to the read identifier.

**Database.** Abarenkov K., Zirk A., Piirmann T., Pöhönen R., Ivanov F., Nilsson R.H., Kõljalg U. (2025). *UNITE general FASTA release for Fungi*. Version 19.02.2025. UNITE Community. https://doi.org/10.15156/BIO/3301229. The release is distributed under the Creative Commons Attribution 4.0 licence (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/). The file in this folder is a selection of 1,041 of its 102,137 sequences and is distributed under the same licence.
