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

`src/pipeline.py` is the recommended entry point. It chains
build -> classify -> diversity -> plot in one command, and stops at the
first stage that fails instead of continuing on bad output:

```bash
python src/pipeline.py \
  --genome-dir data/reference/Genomes \
  --reads path/to/reads.fastq \
  --output-dir results/ \
  --k 21 \
  --source SRA
```

Everything it produces (index, classifications CSV, redistributed-abundance
CSV, SQLite database, diversity report, species-abundance CSV, abundance
plot) lands in `--output-dir`. A few flags worth knowing about:

- `--index path/to/kmer_index.pkl` reuses an existing index instead of
  rebuilding it. It's independent of `--genome-dir`. If you pass both,
  the genomes still get loaded for the GC-outlier check below, they just
  aren't used to rebuild the index.
- `--min-quality 20` drops FASTQ reads whose mean Phred quality (Phred+33)
  is below the threshold before classifying. FASTQ only, it'll raise if
  you use it on a FASTA file.
- `--top-n 10` controls how many species get plotted individually before
  the rest collapse into "Other" (default: 10).
- `--skip-plot` skips the plotting stage.

If `--genome-dir` is given, those reference genomes also get checked for
GC-content outliers (`src/qc.py`). You'll see a warning if a species'
genome GC% is a statistical outlier relative to the rest of the reference
set, since that's a plausible contributor to abundance skew (see
"Validated accuracy" below).

Each stage still works standalone if you want more control over an
individual step: `src/build_reference.py`, `src/classify_reads.py`
(also takes `--min-quality`), and `src/diversity_report.py`. Run any of
them with `--help` for their options.

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
- **Splitting votes for better abundance estimates**: alongside the
  per-read winner-take-all CSV, `classify_reads.py` also splits each
  read's k-mer votes proportionally across every species it hit, instead
  of giving the whole read to one winner. See "Validated accuracy" below
  for the measured effect.

### Validated accuracy

Classified 390,381 real Illumina reads (SRA accession SRR10391187) against a
10-species ZymoBIOMICS D6300 reference index (8 bacteria at 12% each, 2
yeasts at 2% each, k=21, 68.66M k-mers).

**Redistributed bacterial abundance lands within 12.3% mean relative
deviation of ground truth, better than [Zymo's own published <15% baseline](https://files.zymoresearch.com/protocols/_d6300_zymobiomics_microbial_community_standard.pdf)
from their in-house shotgun sequencing.** That's despite SRR10391187 turning
out to be 16S rRNA amplicon sequencing, not whole-genome shotgun as
originally assumed. That's more challenging for a whole-genome k-mer
classifier, since each read only covers a ~588bp PCR-amplified slice of the
genome rather than a uniform sample of it.

The two yeasts (Saccharomyces cerevisiae, Cryptococcus neoformans) landed at
roughly 1/10th to 1/30th of their expected 2% each, and vote-splitting
barely moved them. Confirmed via SRA metadata: bacteria-targeted 16S
primers generally don't amplify fungal rRNA, so most yeast DNA was never in
the sequenced pool to begin with.

All 10 species were still correctly detected (no phantoms), 96.1% of reads
classified.

### Validation graphic

![Redistributed vs. expected abundance for each species in the ZymoBIOMICS mock community](docs/abundance_validation.png)

Same validation run as above, plotted against the known ground truth. This
simply validates the tool's accuracy, it isn't a core component of the
pipeline. It only works because ZymoBIOMICS publishes an expected
composition to compare against, which your own data won't have.
Regenerate it after a fresh classification run with
`python docs/plot_zymobiomics_validation.py` (dataset-specific).

### Performance baseline

End-to-end timing (index load + classification, not just classification in
isolation) on the 390,381-read / 30.5M-k-mer benchmark above:

- Index load: ~19–22s
- Classification: ~56–60s (~6,450–6,940 reads/sec)
- Total: ~80s

See `src/benchmark.py`.

### Docker

A `Dockerfile` is included (`python:3.13-slim`, entrypoint wraps
`pipeline.py`, reference genomes baked into the image). Verified: a
container build reproduces the pipeline's output identically to a
bare-metal run, on the 8-species pre-redistribution baseline (see
"Validated accuracy" above for current numbers):

```bash
docker build -t sequence-to-taxa-processor .
docker run --rm -v ${PWD}/results:/app/results sequence-to-taxa-processor \
  --genome-dir data/reference/Genomes \
  --reads /app/reads.fastq \
  --output-dir /app/results
```

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
│   ├── diversity.py            # species richness, Shannon diversity
│   ├── qc.py                   # GC-outlier check, read quality filtering
│   ├── visualization.py        # species-abundance-by-sample chart
│   └── pipeline.py             # CLI: unified build->classify->diversity->plot
├── data/
│   ├── reference/Genomes/      # example reference genomes (ZymoBIOMICS)
│   └── raw/sra_reads/          # example real read data (gitignored)
├── tests/                      # pytest suite (fasta_utils, classifier,
│                                #   database, diversity, qc, build_reference,
│                                #   classify_reads, visualization, pipeline)
├── docs/                       # portfolio assets referenced by this README,
│                                #   + the script that regenerates them
├── Dockerfile                  # containerized pipeline.py entrypoint
└── requirements.txt
```

### Project status

- **Done**: core k-mer classifier, proportional vote-splitting for
  abundance estimates, real FASTA/FASTQ parsing, SQLite storage with
  per-sample abundance queries, diversity metrics (species richness,
  Shannon index), GC-outlier and read-quality QC (`qc.py`), a
  species-abundance-by-sample visualization, a unified `pipeline.py` CLI,
  CI running the test suite on every push, Docker containerization
  (see above).
- **Not yet started**: Nextflow workflow orchestration, cloud deployment,
  polished demo notebooks.

### Validation dataset

Validation (not core to the tool) uses the **ZymoBIOMICS Microbial Community
Standard (D6300)**, a commercially defined mock community with known
composition, BioProject PRJNA587452, SRA accession SRR10391187. Reference
genomes for its 8 bacterial species and 2 yeast species are committed under
`data/reference/Genomes/`.
