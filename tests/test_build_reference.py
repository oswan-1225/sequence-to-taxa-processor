from build_reference import load_all_genomes


def write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def test_load_all_genomes_merges_multi_record_genome(tmp_path):
    write_fasta(tmp_path / "species_a.fna", [
        ("chromosome", "ACGT"),
        ("plasmid", "TTTT"),
    ])

    genomes = load_all_genomes(str(tmp_path))

    assert genomes == {"species_a": "ACGTTTTT"}


def test_load_all_genomes_keys_by_filename_stem_per_species(tmp_path):
    write_fasta(tmp_path / "species_a.fna", [("chr", "ACGT")])
    write_fasta(tmp_path / "species_b.fasta", [("chr", "GGCC")])

    genomes = load_all_genomes(str(tmp_path))

    assert set(genomes) == {"species_a", "species_b"}
    assert genomes["species_a"] == "ACGT"
    assert genomes["species_b"] == "GGCC"


def test_load_all_genomes_ignores_unrecognized_extensions(tmp_path):
    write_fasta(tmp_path / "species_a.fna", [("chr", "ACGT")])
    (tmp_path / "notes.txt").write_text("not a sequence file")

    genomes = load_all_genomes(str(tmp_path))

    assert set(genomes) == {"species_a"}
