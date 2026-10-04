"""B5 prior-work evidence files: structural consistency checks (no scientific numerics)."""

import hashlib
import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("b5_evidence", ROOT / "tools" / "b5_evidence.py")
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)

B5_DOCS = ("B5_PRIOR_WORK_EXTRACTION.md", "B5_PRIOR_WORK_MATRIX.md", "B5_EVIDENCE_LEDGER.md", "B5_MATH_COMPARISON.md",
           "B5_FOLLOWUP_SOURCES.md", "B5_POTENTIAL_ISSUES.md")


@pytest.fixture(scope="module")
def ledger():
    return E.parse_ledger()


@pytest.fixture(scope="module")
def matrix():
    return E.parse_matrix()


def test_ledger_entries_complete(ledger):
    assert len(ledger) > 0
    for eid, r in ledger.items():
        assert all(r[f] is not None for f in E.FIELDS), eid
        assert r["Paper"] in ("A", "B") and eid.startswith(r["Paper"] + "-"), eid             # every row has a paper id
        assert r["Classification"] in E.ALLOWED, (eid, r["Classification"])
        locator = [r[k] for k in ("Section", "Subsection", "Equation", "Figure", "Appendix", "Page")]
        assert any(v not in E.EMPTY for v in locator), eid                                   # every row has a location
        assert r["Source version"] not in E.EMPTY and r["Source URL"].startswith("https://"), eid
        if r["Paper"] == "B":                                                                # never pretend B = published
            assert "arXiv:2507.22054v2" in r["Source version"], eid
        if r["Classification"] == "NOT LOCATED":
            assert "Not located in the reviewed version" in r["Source statement (paraphrase)"], eid
            assert "searches for:" in r["Source statement (paraphrase)"], eid


def test_matrix_cells_link_to_matching_ledger_entries(matrix, ledger):
    required = {str(i) for i in range(1, 26) if i != 14} | {"14a", "14b"}
    assert required <= {row["row"] for row in matrix}
    for row in matrix:
        for paper in ("A", "B"):
            label, ids, cell = row[paper]
            assert label in E.ALLOWED, (row["row"], paper, cell)
            assert ids, (row["row"], paper, "cell has no ledger reference")
            for i in ids:
                assert i in ledger, (row["row"], i)
                assert ledger[i]["Paper"] == paper, (row["row"], i)
                assert ledger[i]["Classification"] == label, (row["row"], i, ledger[i]["Classification"], label)
                assert row["row"] in E.matrix_rows_of(ledger[i]), (row["row"], i)


def test_every_ledger_matrix_row_is_linked_back(matrix, ledger):
    linked = {(row["row"], i) for row in matrix for p in ("A", "B") for i in row[p][1]}
    for eid, r in ledger.items():
        for rr in E.matrix_rows_of(r):
            assert (rr, eid) in linked, (eid, rr)


def test_each_candidate_has_an_evidence_trail(ledger):
    for cand in ("C1", "C2", "C3", "C4"):
        hits = [e for e, r in ledger.items() if cand in re.split(r"[;,\s]+", r["Candidate"])]
        assert len(hits) >= 2, cand          # at least one entry per paper
        assert {ledger[h]["Paper"] for h in hits} == {"A", "B"}, cand


def test_csv_matches_ledger(ledger):
    rows = E.read_csv()
    assert tuple(rows[0].keys()) == E.CSV_COLUMNS
    assert [{k: str(v) for k, v in r.items()} for r in E.ledger_to_rows(ledger)] == rows
    for r in rows:
        assert r["paper_id"] in ("thanasilp2024", "aghaeisaem2026")
        assert r["classification"] in E.ALLOWED


def test_manifest_complete_and_hashes_match():
    man = json.loads(E.MANIFEST.read_text(encoding="utf-8"))
    assert {p["paper_id"] for p in man["papers"]} == {"thanasilp2024", "aghaeisaem2026"}
    for p in man["papers"]:
        for key in ("title", "authors", "journal", "year", "doi", "arxiv_id", "version_reviewed", "publication_date",
                    "review_completed", "files"):
            assert p.get(key) not in (None, "", []), (p["paper_id"], key)
        for f in p["files"]:
            assert f["download_source"].startswith("https://") and len(f["sha256"]) == 64
            assert f["local_filename"].startswith("research_sources/")
            local = ROOT / f["local_filename"]
            if local.exists():                    # local copies are git-ignored; check only where present
                assert hashlib.sha256(local.read_bytes()).hexdigest() == f["sha256"], f["local_filename"]
    b = next(p for p in man["papers"] if p["paper_id"] == "aghaeisaem2026")
    assert b["published_version_inspected"] is False


def test_no_novelty_classification_or_language():
    for name in B5_DOCS + ("results/b5_prior_work/evidence.csv",):
        text = (ROOT / name).read_text(encoding="utf-8")
        assert "NOVEL" not in text.replace("NOVELTY", ""), name
        assert not re.search(r"\bnovel\b", text, re.I), name


def test_research_sources_are_git_ignored():
    assert "research_sources/" in (ROOT / ".gitignore").read_text(encoding="utf-8")
