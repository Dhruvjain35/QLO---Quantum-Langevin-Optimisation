"""B5 prior-work evidence files: parse the ledger and matrix, and build results/b5_prior_work/evidence.csv.

    python tools/b5_evidence.py            # rewrite evidence.csv from B5_EVIDENCE_LEDGER.md
    python tools/b5_evidence.py --check    # only verify that evidence.csv matches the ledger

The ledger is the single source of truth; the CSV is generated from it so the two cannot drift.
No scientific code is involved.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "B5_EVIDENCE_LEDGER.md"
MATRIX = ROOT / "B5_PRIOR_WORK_MATRIX.md"
CSV_PATH = ROOT / "results" / "b5_prior_work" / "evidence.csv"
MANIFEST = ROOT / "results" / "b5_prior_work" / "source_manifest.json"

ALLOWED = ("EXPLICIT", "PARTIAL / RELATED", "PARTIAL / IMPLIED", "NOT LOCATED")
PAPERS = {"A": ("thanasilp2024", "Exponential concentration in quantum kernel methods", 2024),
          "B": ("aghaeisaem2026", "Pitfalls when tackling the exponential concentration of parameterized quantum models", 2026)}
FIELDS = ("Paper", "Candidate", "Matrix rows", "Claim category", "Classification", "Section", "Subsection", "Equation",
          "Figure", "Appendix", "Page", "Source version", "Source URL", "Short quote", "Source statement (paraphrase)",
          "Mathematical expression", "Assumptions", "Scope", "Relation to Stage 7", "Does NOT establish")
CSV_COLUMNS = ("evidence_id", "paper_id", "paper_title", "year", "candidate_id", "matrix_rows", "claim_category",
               "classification", "section", "subsection", "equation", "figure", "appendix", "page", "short_quote",
               "paraphrase", "assumptions", "scope", "relation_to_stage7", "does_not_establish", "source_version",
               "source_url")
EMPTY = {"", "—", "-"}


def parse_ledger(path: Path = LEDGER) -> dict[str, dict]:
    text = path.read_text(encoding="utf-8")
    entries: dict[str, dict] = {}
    parts = re.split(r'<a id="([AB]-[A-Za-z0-9]+)"></a>', text)
    for eid, body in zip(parts[1::2], parts[2::2]):
        if eid in entries:
            raise ValueError(f"duplicate ledger id {eid}")
        head = re.search(r"^###\s+(.+)$", body, re.M)
        rec = {"id": eid, "title": head.group(1).strip() if head else ""}
        for name in FIELDS:
            m = re.search(r"^- \*\*" + re.escape(name) + r":\*\*\s*(.*)$", body, re.M)
            rec[name] = m.group(1).strip() if m else None
        entries[eid] = rec
    return entries


def matrix_rows_of(rec: dict) -> list[str]:
    raw = rec.get("Matrix rows") or ""
    return [] if raw.strip() in EMPTY else [r.strip() for r in raw.split(",") if r.strip()]


def parse_matrix(path: Path = MATRIX) -> list[dict]:
    """Rows of the main table: number, claim, and per-paper (label, [linked ids])."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not re.match(r"^\|\s*\d+[ab]?\s*\|", line):
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
        num, claim, cell_a, cell_b = cells[0], cells[1], cells[2], cells[3]
        out = {"row": num, "claim": claim}
        for key, cell in (("A", cell_a), ("B", cell_b)):
            label = next((lab for lab in sorted(ALLOWED, key=len, reverse=True) if cell.startswith(lab)), None)
            out[key] = (label, re.findall(r"\[([AB]-[A-Za-z0-9]+)\]\(B5_EVIDENCE_LEDGER\.md#\1\)", cell), cell)
        rows.append(out)
    return rows


def ledger_to_rows(entries: dict[str, dict]) -> list[dict]:
    rows = []
    for eid, r in entries.items():
        pid, title, year = PAPERS[r["Paper"]]
        rows.append({"evidence_id": eid, "paper_id": pid, "paper_title": title, "year": year,
                     "candidate_id": r["Candidate"], "matrix_rows": r["Matrix rows"], "claim_category": r["Claim category"],
                     "classification": r["Classification"], "section": r["Section"], "subsection": r["Subsection"],
                     "equation": r["Equation"], "figure": r["Figure"], "appendix": r["Appendix"], "page": r["Page"],
                     "short_quote": r["Short quote"], "paraphrase": r["Source statement (paraphrase)"],
                     "assumptions": r["Assumptions"], "scope": r["Scope"], "relation_to_stage7": r["Relation to Stage 7"],
                     "does_not_establish": r["Does NOT establish"], "source_version": r["Source version"],
                     "source_url": r["Source URL"]})
    return rows


def write_csv(rows: list[dict], path: Path = CSV_PATH) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(rows)


def read_csv(path: Path = CSV_PATH) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    rows = ledger_to_rows(parse_ledger())
    if a.check:
        on_disk = read_csv()
        same = [{k: str(v) for k, v in r.items()} for r in rows] == on_disk
        print("evidence.csv matches ledger" if same else "evidence.csv is STALE")
        return 0 if same else 1
    write_csv(rows)
    print(f"wrote {len(rows)} rows -> {CSV_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
