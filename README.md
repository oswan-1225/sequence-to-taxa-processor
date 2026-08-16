## Sequence to Taxa Processor

A from-scratch metagenomic taxonomic classifier: given raw DNA sequencing
reads (FASTA/FASTQ) and a set of reference genomes, classify each read to
its most likely species of origin using k-mer matching all in Python.

The tool is general-purpose and bring-your-own-data: point it at any folder
of reference genomes and any read file via CLI arguments. It does not
hardcode a specific species or dataset.

### Installation

```bash
pip install -r requirements.txt
```

### Quickstart

Build a k-mer reference index from a folder of reference genomes:

```bash
python src/build_reference.py \
  --genome_dir data/reference/Genomes \
  --output data/reference/kmer_index.pkl \
  --k 21
```

Classify a read file against that index:

```bash
python src/classify_reads.py \
  --index data/reference/kmer_index.pkl \
  --reads path/to/reads.fastq \
  --output results.csv \
  --db results.db \
  --source SRA
```

`--db`/`--source` are optional; when given, results are also written to a
SQLite database (see `src/database.py`) for cross-sample abundance queries.

Summarize per-sample diversity (species richness, Shannon index) from that
database:

```bash
python src/diversity_report.py --db results.db --output diversity.csv
```

### How it works

- **k-mer indexing**: every reference genome is broken into overlapping
  k-mers (default k=21) and indexed as `{kmer: set of species}`.
- **Classification**: each read is broken into the same k-mers, and the
  species with the most matching k-mers wins ("top hit"), with a confidence
  score based on vote share.
- **Why k=21**: short k-mers (k=3/4) are statistically near-guaranteed to
  appear in many genomes by chance (DNA has only 4 letters), making them
  non-discriminating. Real tools like Kraken2 use k in the 21-35 range for
  the same reason.

### Validated accuracy

Classified 390,381 real Illumina reads (SRA accession SRR10391187, from the
ZymoBIOMICS D6300 mock community — 8 bacterial species at ~12% each) against
an 8-species reference index (k=21, 30.5M k-mers):

- All 8 expected species correctly detected, no phantom species
- 95.9% of reads classified (4.1% unclassified)
- Per-species abundance in the 4.9%–18.2% range (theoretical truth: 12%
  each) — deviations are explainable by genome size, GC content/sequencing
  bias, and ambiguous k-mer vote-sharing, not unexplained noise

### Performance baseline

End-to-end timing (index load + classification, not just classification in
isolation) on the 390,381-read / 30.5M-k-mer benchmark above:

- Index load: ~19–22s
- Classification: ~56–60s (~6,450–6,940 reads/sec)
- Total: ~80s

See `src/benchmark.py`.

### Repo structure

```
sequence-to-taxa-processor/
├── src/
│   ├── fasta_utils.py          # FASTA/FASTQ parsing, k-mer extraction
│   ├── classifier_functions.py # index building, read classification
│   ├── build_reference.py      # CLI: build a k-mer index from a genome folder
│   ├── classify_reads.py       # CLI: classify reads against a saved index
│   ├── diversity_report.py     # CLI: per-sample diversity report from the DB
│   ├── benchmark.py            # end-to-end performance benchmark
│   ├── database.py             # SQLite storage + abundance queries
│   └── diversity.py            # species richness, Shannon diversity
├── data/
│   ├── reference/Genomes/      # example reference genomes (ZymoBIOMICS)
│   └── raw/sra_reads/          # example real read data (gitignored)
├── tests/                      # pytest suite (fasta_utils, classifier,
│                                #   database, diversity, build_reference)
└── requirements.txt
```

### Project status

- **Done**: core k-mer classifier, real FASTA/FASTQ parsing, SQLite storage
  with per-sample abundance queries, diversity metrics (species richness,
  Shannon index), CI running the test suite on every push.
- **Not yet started**: Nextflow workflow orchestration, Docker
  containerization, cloud deployment, and polished demo notebooks.

### Validation dataset

Validation (not core to the tool) uses the **ZymoBIOMICS Microbial Community
Standard (D6300)**, a commercially defined mock community with known
composition — BioProject PRJNA587452, SRA accession SRR10391187. Reference
genomes for its 8 bacterial species are committed under
`data/reference/Genomes/`.
