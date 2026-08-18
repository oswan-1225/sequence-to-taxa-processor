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
- **Vote redistribution (experimental)**: alongside the per-read
  winner-take-all CSV, `classify_reads.py` also splits each read's k-mer
  votes proportionally across every species it hit, instead of giving the
  whole read to one winner. See "Validated accuracy" below, it measures
  worse than winner-take-all on this dataset.

### Validated accuracy

Classified 390,381 real Illumina reads (SRA accession SRR10391187) against a
10-species ZymoBIOMICS D6300 reference index (k=21, 68.66M k-mers).

SRR10391187 is 16S rRNA amplicon sequencing (confirmed via SRA's own run
metadata), so accuracy is measured against Zymo's own 16S-adjusted
theoretical composition, avoiding a flat per-species split, since rRNA copy
number varies a lot between species:

| Species | Genomic DNA | 16S Only | 16S & 18S |
|---|---|---|---|
| Pseudomonas aeruginosa | 12% | 4.2% | 3.6% |
| Escherichia coli | 12% | 10.1% | 8.9% |
| Salmonella enterica | 12% | 10.4% | 9.1% |
| Lactobacillus fermentum | 12% | 18.4% | 16.1% |
| Enterococcus faecalis | 12% | 9.9% | 8.7% |
| Staphylococcus aureus | 12% | 15.5% | 13.6% |
| Listeria monocytogenes | 12% | 14.1% | 12.4% |
| Bacillus subtilis | 12% | 17.4% | 15.3% |
| Saccharomyces cerevisiae | 2% | NA | 9.3% |
| Cryptococcus neoformans | 2% | NA | 3.3% |

(Source: [Zymo's D6300 datasheet](https://files.zymoresearch.com/protocols/_d6300_zymobiomics_microbial_community_standard.pdf),
Table 1. "16S Only" is copy-number-adjusted for standard bacterial 16S
primers, which is why the yeasts show `NA`, not `0%`; those primers don't
amplify fungal rRNA at all.)

Winner-take-all classification (each read's full weight goes to its best
match) lands within 14.9% mean relative deviation of that target across
the 8 bacteria. The two yeasts are basically invisible in this data, which is intentional.  
Bacterial 16S primers don't pick up fungal
rRNA, and Zymo's own 16S numbers list the yeasts as not applicable. Our
observed 0.06% and 0.23% line up with that.

Two alternatives were tested and rejected in favor of winner-take-all.
Proportional vote-share redistribution (splitting a read's credit across
every species it hits, instead of giving it all to one winner) comes in
worse, at 27.3% mean deviation. Discarding ambiguous multi-genome reads
entirely, Zymo's own suggested fix for this kind of ambiguity, throws away
93.7% of the reads here and still lands at 138.7%, since with 16S amplicon
data almost every read overlaps a region shared across species.
Winner-take-all wins clearly, so that's the reported method.

### Validation graphic

![Observed vs. Zymo's 16S-adjusted expected abundance for each bacterial species in the ZymoBIOMICS mock community](docs/abundance_validation.png)

Same run as above, plotted against Zymo's 16S-adjusted composition. Yeasts
aren't included, no valid 16S target to compare them against. Not a core
part of the pipeline either; this only works because Zymo happens to
publish real reference data to check against. Regenerate it after a fresh
classification run with
`python docs/plot_zymobiomics_validation.py` (dataset-specific).

### Performance baseline

End-to-end timing on the 390,381-read / 68.66M-k-mer benchmark above:

- Index load: ~56s
- Classification: ~68s (~5,700 reads/sec)
- Total: ~134s

See `src/benchmark.py`.

### Docker

A `Dockerfile` is included (`python:3.13-slim`, entrypoint wraps
`pipeline.py`, reference genomes baked into the image). A container build
reproduces the pipeline's output identically to a bare-metal run, checked
against the 8-species benchmark (see "Validated accuracy" above for the
full 10-species numbers):

```bash
docker build -t sequence-to-taxa-processor
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
