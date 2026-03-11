"""Integration tests for Kraken2Runner using a real mini database.

These tests build a tiny kraken2 database containing one E. coli
(taxid 562) and one Homo sapiens (taxid 9606) reference sequence,
then classify mock reads and verify that bacterial, human, and
unclassified reads are correctly tallied.

Requirements:
    * ``kraken2`` and ``kraken2-build`` must be on PATH.
    * No network access required — all reference sequences and
      taxonomy files are generated locally.

The database is built once per session via a ``session``-scoped
pytest fixture and reused across all tests.
"""

import os
import random
import shutil
import subprocess
import tempfile

import pytest

from kmer_denovo_filter.kmer_utils import Kraken2Runner

_KRAKEN2_AVAILABLE = shutil.which("kraken2") is not None
_KRAKEN2_BUILD_AVAILABLE = shutil.which("kraken2-build") is not None

pytestmark = pytest.mark.skipif(
    not (_KRAKEN2_AVAILABLE and _KRAKEN2_BUILD_AVAILABLE),
    reason="kraken2 and kraken2-build must be on PATH",
)

# ── Taxonomy constants ────────────────────────────────────────────
_ECOLI_TAXID = 562
_HUMAN_TAXID = 9606


def _gen_seq(length, seed_val):
    """Generate a deterministic pseudo-random DNA sequence."""
    rng = random.Random(seed_val)
    return "".join(rng.choice("ACGT") for _ in range(length))


# ── Fixtures ──────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def mini_kraken2_db(tmp_path_factory):
    """Build a tiny kraken2 database with E. coli and Human sequences.

    The database is created once per test session and shared across all
    integration tests.  It uses a minimal taxonomy (root → Bacteria →
    E. coli, root → Eukaryota → Homo sapiens) and two short reference
    sequences (2000 bp each).
    """
    db = str(tmp_path_factory.mktemp("mini_kraken2_db"))

    # ── Taxonomy ──────────────────────────────────────────────────
    tax_dir = os.path.join(db, "taxonomy")
    os.makedirs(tax_dir, exist_ok=True)

    # nodes.dmp — full lineages so kraken2-build can resolve ancestors.
    nodes = [
        # (taxid, parent, rank)
        (1, 1, "no rank"),
        (131567, 1, "no rank"),          # cellular organisms
        # Bacteria lineage → E. coli
        (2, 131567, "superkingdom"),
        (1224, 2, "phylum"),
        (1236, 1224, "class"),
        (91347, 1236, "order"),
        (543, 91347, "family"),
        (561, 543, "genus"),
        (562, 561, "species"),
        # Eukaryota lineage → Homo sapiens
        (2759, 131567, "superkingdom"),
        (33154, 2759, "no rank"),
        (33208, 33154, "kingdom"),
        (6072, 33208, "no rank"),
        (33213, 6072, "no rank"),
        (33511, 33213, "no rank"),
        (7711, 33511, "phylum"),
        (89593, 7711, "subphylum"),
        (7742, 89593, "no rank"),
        (7776, 7742, "superclass"),
        (117570, 7776, "no rank"),
        (117571, 117570, "no rank"),
        (8287, 117571, "no rank"),
        (32523, 8287, "no rank"),
        (32524, 32523, "no rank"),
        (40674, 32524, "class"),
        (32525, 40674, "no rank"),
        (9347, 32525, "no rank"),
        (1437010, 9347, "no rank"),
        (314146, 1437010, "superorder"),
        (9443, 314146, "order"),
        (376913, 9443, "suborder"),
        (314293, 376913, "infraorder"),
        (9526, 314293, "parvorder"),
        (314295, 9526, "superfamily"),
        (9604, 314295, "family"),
        (207598, 9604, "subfamily"),
        (9605, 207598, "genus"),
        (9606, 9605, "species"),
    ]
    with open(os.path.join(tax_dir, "nodes.dmp"), "w") as fh:
        for taxid, parent, rank in nodes:
            fh.write(
                f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0"
                f"\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            )

    names = {
        1: "root", 131567: "cellular organisms",
        2: "Bacteria", 1224: "Pseudomonadota",
        1236: "Gammaproteobacteria", 91347: "Enterobacterales",
        543: "Enterobacteriaceae", 561: "Escherichia",
        562: "Escherichia coli",
        2759: "Eukaryota", 33154: "Opisthokonta",
        33208: "Metazoa", 6072: "Eumetazoa",
        33213: "Bilateria", 33511: "Deuterostomia",
        7711: "Chordata", 89593: "Craniata",
        7742: "Vertebrata", 7776: "Gnathostomata",
        117570: "Teleostomi", 117571: "Euteleostomi",
        8287: "Sarcopterygii", 32523: "Tetrapoda",
        32524: "Amniota", 40674: "Mammalia",
        32525: "Theria", 9347: "Eutheria",
        1437010: "Boreoeutheria", 314146: "Euarchontoglires",
        9443: "Primates", 376913: "Haplorrhini",
        314293: "Simiiformes", 9526: "Catarrhini",
        314295: "Hominoidea", 9604: "Hominidae",
        207598: "Homininae", 9605: "Homo",
        9606: "Homo sapiens",
    }
    with open(os.path.join(tax_dir, "names.dmp"), "w") as fh:
        for taxid, name in sorted(names.items()):
            fh.write(
                f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
            )

    # ── Reference sequences (2000 bp each) ──────────────────────
    ecoli_seq = _gen_seq(2000, seed_val=_ECOLI_TAXID)
    human_seq = _gen_seq(2000, seed_val=_HUMAN_TAXID)

    ref_dir = os.path.join(db, "library_src")
    os.makedirs(ref_dir, exist_ok=True)

    ecoli_fa = os.path.join(ref_dir, "ecoli.fna")
    with open(ecoli_fa, "w") as fh:
        fh.write(
            f">ecoli_ref|kraken:taxid|{_ECOLI_TAXID}"
            f" Escherichia coli mock\n"
        )
        for i in range(0, len(ecoli_seq), 70):
            fh.write(ecoli_seq[i:i + 70] + "\n")

    human_fa = os.path.join(ref_dir, "human.fna")
    with open(human_fa, "w") as fh:
        fh.write(
            f">human_ref|kraken:taxid|{_HUMAN_TAXID}"
            f" Homo sapiens mock\n"
        )
        for i in range(0, len(human_seq), 70):
            fh.write(human_seq[i:i + 70] + "\n")

    # ── Build the database ────────────────────────────────────────
    for fa in (ecoli_fa, human_fa):
        subprocess.run(
            [
                "kraken2-build", "--add-to-library", fa,
                "--db", db, "--no-masking",
            ],
            check=True, capture_output=True,
        )

    subprocess.run(
        [
            "kraken2-build", "--build", "--db", db,
            "--threads", "2", "--no-masking",
            "--kmer-len", "35", "--minimizer-len", "31",
        ],
        check=True, capture_output=True,
    )

    # Sanity: required database files must exist.
    for fname in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        assert os.path.isfile(os.path.join(db, fname)), (
            f"Database build failed: {fname} missing"
        )

    return db, ecoli_seq, human_seq


# ── Integration tests ─────────────────────────────────────────────


class TestKraken2Integration:
    """End-to-end classification tests with a real mini database."""

    def test_bacterial_reads_classified(self, mini_kraken2_db):
        """Reads from E. coli reference are classified as bacterial."""
        db, ecoli_seq, _human_seq = mini_kraken2_db
        kr = Kraken2Runner(db)

        # Extract 100 bp reads from the E. coli reference
        sequences = {}
        for i in range(3):
            sequences[f"ecoli_read_{i}"] = ecoli_seq[i * 100:(i + 1) * 100]

        result = kr.classify_sequences(sequences)

        assert result.total == 3
        assert result.classified > 0, (
            "Expected at least some bacterial reads to be classified"
        )
        assert result.bacterial_count > 0, (
            "Expected bacterial reads from E. coli sequence"
        )
        assert result.human_count == 0

    def test_human_reads_classified(self, mini_kraken2_db):
        """Reads from human reference are classified as human."""
        db, _ecoli_seq, human_seq = mini_kraken2_db
        kr = Kraken2Runner(db)

        sequences = {}
        for i in range(3):
            sequences[f"human_read_{i}"] = human_seq[i * 100:(i + 1) * 100]

        result = kr.classify_sequences(sequences)

        assert result.total == 3
        assert result.classified > 0, (
            "Expected at least some human reads to be classified"
        )
        assert result.human_count > 0, (
            "Expected human reads from Homo sapiens sequence"
        )
        assert result.bacterial_count == 0

    def test_novel_reads_unclassified(self, mini_kraken2_db):
        """Random sequences not in the database are unclassified."""
        db, _ecoli_seq, _human_seq = mini_kraken2_db
        kr = Kraken2Runner(db)

        sequences = {}
        for i in range(3):
            sequences[f"novel_{i}"] = _gen_seq(100, seed_val=99999 + i)

        result = kr.classify_sequences(sequences)

        assert result.total == 3
        assert result.unclassified > 0, (
            "Expected novel sequences to be largely unclassified"
        )

    def test_mixed_reads(self, mini_kraken2_db):
        """Mixed input correctly separates bacterial, human, and novel."""
        db, ecoli_seq, human_seq = mini_kraken2_db
        kr = Kraken2Runner(db)

        sequences = {
            "bact_1": ecoli_seq[0:100],
            "bact_2": ecoli_seq[100:200],
            "bact_3": ecoli_seq[200:300],
            "human_1": human_seq[0:100],
            "human_2": human_seq[100:200],
            "human_3": human_seq[200:300],
            "novel_1": _gen_seq(100, seed_val=88888),
            "novel_2": _gen_seq(100, seed_val=77777),
        }

        result = kr.classify_sequences(sequences)

        assert result.total == 8
        # All reads should be accounted for
        assert result.classified + result.unclassified == result.total

        # Bacterial reads should be detected
        assert result.bacterial_count > 0
        # Human reads should be detected
        assert result.human_count > 0
        # Novel reads should mostly be unclassified
        assert result.unclassified > 0

        # Bacterial reads should be named
        bact_names = {n for n in result.bacterial_read_names}
        assert all(n.startswith("bact_") for n in bact_names)

        # Bacterial fraction should be reasonable
        assert 0.0 < result.bacterial_fraction < 1.0

    def test_taxonomy_lineage_resolves_ecoli(self, mini_kraken2_db):
        """Taxonomy lineage correctly identifies E. coli as bacterial.

        The key scenario that was failing in production: kraken2
        classifies a read as taxid 562 (E. coli), but without the
        taxonomy tree, only exact taxid==2 would be counted.  With
        the tree, 562 is correctly resolved as a descendant of 2.
        """
        db, ecoli_seq, _human_seq = mini_kraken2_db

        # Verify the taxonomy file is found and parsed
        bacterial_taxids = Kraken2Runner._load_bacterial_taxids(db)
        assert bacterial_taxids is not None, (
            "Taxonomy should be loaded from the database"
        )
        assert _ECOLI_TAXID in bacterial_taxids, (
            "E. coli (562) should be identified as bacterial"
        )
        assert _HUMAN_TAXID not in bacterial_taxids, (
            "Human (9606) should NOT be identified as bacterial"
        )

    def test_confidence_threshold(self, mini_kraken2_db):
        """Higher confidence threshold reduces classifications."""
        db, ecoli_seq, human_seq = mini_kraken2_db

        sequences = {
            "bact_1": ecoli_seq[0:100],
            "human_1": human_seq[0:100],
            "novel_1": _gen_seq(100, seed_val=55555),
        }

        result_low = Kraken2Runner(db, confidence=0.0).classify_sequences(
            sequences,
        )
        result_high = Kraken2Runner(db, confidence=0.9).classify_sequences(
            sequences,
        )

        # Higher confidence should classify fewer or equal reads
        assert result_high.classified <= result_low.classified

    def test_result_summary_format(self, mini_kraken2_db):
        """Result summary string is well-formed."""
        db, ecoli_seq, _human_seq = mini_kraken2_db
        kr = Kraken2Runner(db)

        result = kr.classify_sequences({
            "r1": ecoli_seq[0:100],
        })

        summary = result.summary()
        assert "kraken2:" in summary
        assert "reads" in summary
        assert "classified" in summary
        assert "bacterial" in summary


class TestKraken2RootLevelNodesDmp:
    """Test taxonomy loading from root-level nodes.dmp (PrackenDB layout)."""

    def test_classify_with_root_level_nodes_dmp(self, mini_kraken2_db):
        """Classification works when nodes.dmp is at the DB root.

        PrackenDB places nodes.dmp directly in the database directory,
        not under a taxonomy/ subdirectory.  This test verifies that
        Kraken2Runner can still load the taxonomy and correctly
        identify bacterial reads.
        """
        db, ecoli_seq, _human_seq = mini_kraken2_db

        # Create a copy of the database with nodes.dmp at root level
        with tempfile.TemporaryDirectory() as flat_db:
            # Copy database files
            for fname in os.listdir(db):
                src = os.path.join(db, fname)
                dst = os.path.join(flat_db, fname)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)

            # Copy nodes.dmp from taxonomy/ to root, remove taxonomy/
            tax_nodes = os.path.join(db, "taxonomy", "nodes.dmp")
            shutil.copy2(tax_nodes, os.path.join(flat_db, "nodes.dmp"))
            # Don't create taxonomy/ subdirectory — flat layout only

            # Verify taxonomy loads from root-level nodes.dmp
            bacterial = Kraken2Runner._load_bacterial_taxids(flat_db)
            assert bacterial is not None
            assert _ECOLI_TAXID in bacterial

            # Run classification with the flat-layout database
            kr = Kraken2Runner(flat_db)
            result = kr.classify_sequences({
                "ecoli": ecoli_seq[0:100],
            })

            assert result.total == 1
            assert result.classified > 0
            assert result.bacterial_count > 0
