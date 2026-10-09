"""B5 closure audit: targeted literature search (arXiv API, OpenAlex, Semantic Scholar) + forward-citation pull.

    python tools/b5_closure_search.py

Raw API responses go to research_sources/b5_closure/search/raw/ (git-ignored) and every hit, with its query
provenance, to research_sources/b5_closure/search/hits.csv. The committed screening summary is
results/b5_prior_work/closure_search_screen.csv. No scientific code is involved.
"""
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / "research_sources" / "b5_closure" / "search"
RAW = HERE / "raw"
RAW.parent.mkdir(parents=True, exist_ok=True)
RAW.mkdir(exist_ok=True)
UA = {"User-Agent": "qlo-research-b5-audit/1.0 (literature search; no personal data)"}

# Query sets by target. T1 zero probability / exact law; T2 sign statistics; T3 SWAP vs Loschmidt;
# T4 measurement-dependent shot exponent; T5 vector / trajectory.
QUERIES = {
    # T1: exact law / zero probability of finite-shot parameter-shift gradients
    "T1a": ('abs:"parameter shift" AND abs:shot', "parameter shift gradient estimator finite shots distribution"),
    "T1b": ('abs:"parameter-shift" AND abs:binomial', "parameter shift gradient binomial shot noise distribution"),
    "T1c": ('abs:Skellam AND cat:quant-ph', "Skellam distribution quantum gradient estimator"),
    "T1d": ('abs:"finite shot" AND abs:gradient AND abs:zero', "finite shot gradient estimate exactly zero probability"),
    # T2: sign statistics of finite-shot gradients
    "T2a": ('abs:gradient AND abs:sign AND abs:shots AND cat:quant-ph', "sign of gradient estimate finite shots probability correct sign"),
    "T2b": ('abs:"sign" AND abs:"shot noise" AND abs:variational', "gradient sign accuracy shot noise variational quantum algorithm"),
    "T2c": ('abs:signSGD AND cat:quant-ph', "sign-based optimizer variational quantum shot noise"),
    # T3: SWAP vs Loschmidt (inversion / compute-uncompute / echo) fidelity estimation
    "T3a": ('abs:"swap test" AND abs:"inversion test"', "swap test inversion test fidelity estimation shots"),
    "T3b": ('abs:"swap test" AND abs:"compute-uncompute"', "swap test compute-uncompute overlap estimation comparison"),
    "T3c": ('abs:"swap test" AND abs:Loschmidt', "swap test Loschmidt echo fidelity comparison"),
    "T3d": ('abs:"swap test" AND abs:variance', "swap test variance fidelity estimation number of measurements"),
    "T3e": ('abs:"overlap estimation" AND abs:"swap test"', "overlap estimation swap test sample complexity comparison"),
    # T4: measurement-dependent shot exponent / ratio
    "T4a": ('abs:"barren plateau" AND abs:"swap test"', "barren plateau swap test fidelity cost shots"),
    "T4b": ('abs:"quantum kernel" AND abs:"swap test" AND abs:shots', "quantum kernel swap test inversion test measurement shots"),
    "T4c": ('abs:fidelity AND abs:"shot noise" AND abs:exponent', "measurement scheme changes shot complexity exponent fidelity"),
    # T5: vector / trajectory consequences
    "T5a": ('abs:"cosine similarity" AND abs:gradient AND cat:quant-ph', "cosine similarity estimated gradient shot noise variational quantum"),
    "T5b": ('abs:"random walk" AND abs:"barren plateau"', "random walk barren plateau shot noise gradient descent"),
}
# Supplemental arXiv-only queries (second pass, C5 / C2 terms).
SUPPLEMENTAL = {
    "S1": 'abs:"destructive swap test" AND abs:variance',
    "S2": 'abs:"overlap estimation" AND abs:variance',
    "S3": 'abs:"inversion test" AND abs:shots',
    "S4": 'abs:"swap test" AND abs:"barren plateau"',
    "S5": 'abs:"parameter shift" AND abs:fidelity AND abs:shots',
    "S6": 'abs:"relative gradient" AND cat:quant-ph AND abs:shots',
    "S7": 'abs:"sign" AND abs:"parameter shift" AND abs:"shot noise"',
    "S8": 'abs:"Hong-Ou-Mandel" AND abs:"overlap estimation"',
    "S9": 'abs:"shot budget" AND abs:"swap test"',
    "S10": 'abs:"measurement cost" AND abs:fidelity AND abs:gradient AND cat:quant-ph',
}
ANCHORS = {"thanasilp2024": "W4399770418", "aghaeisaem2026": "W4416437554", "teo2023": "W4366086631",
           "mari2021": "W3049070297"}


def get(url, tries=4, wait=4.0):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            code = getattr(e, "code", None)
            if code == 429 or code is None or (isinstance(code, int) and code >= 500):
                time.sleep(wait * (k + 1))
                continue
            raise
    return None


def arxiv(qid, q):
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"search_query": q, "start": 0, "max_results": 50, "sortBy": "relevance"})
    x = get(url)
    (RAW / f"arxiv_{qid}.xml").write_text(x or "", encoding="utf-8")
    out = []
    for e in re.findall(r"<entry>(.*?)</entry>", x or "", re.S):
        g = lambda t: (re.findall(r"<%s[^>]*>(.*?)</%s>" % (t, t), e, re.S) or [""])[0]
        aid = g("id").rsplit("/abs/", 1)[-1]
        out.append({"source": "arXiv", "query": qid, "id": "arXiv:" + aid, "title": re.sub(r"\s+", " ", g("title")).strip(),
                    "year": g("published")[:4], "abstract": re.sub(r"\s+", " ", g("summary")).strip()})
    time.sleep(3.2)  # arXiv API etiquette
    return out


def inv_abstract(inv):
    if not inv:
        return ""
    pos = {}
    for w, ps in inv.items():
        for p in ps:
            pos[p] = w
    return " ".join(pos[i] for i in sorted(pos))


def openalex_search(qid, q):
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(
        {"search": q, "per-page": 50, "select": "id,doi,title,publication_year,abstract_inverted_index,ids"})
    x = get(url)
    (RAW / f"openalex_{qid}.json").write_text(x or "{}", encoding="utf-8")
    d = json.loads(x or "{}")
    out = []
    for w in d.get("results", []):
        out.append({"source": "OpenAlex", "query": qid, "id": w.get("doi") or w.get("id"), "title": w.get("title") or "",
                    "year": str(w.get("publication_year") or ""), "abstract": inv_abstract(w.get("abstract_inverted_index"))})
    time.sleep(1.0)
    return out


def openalex_citers(name, wid):
    out, cursor, page = [], "*", 0
    while cursor and page < 10:
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(
            {"filter": f"cites:{wid}", "per-page": 200, "cursor": cursor,
             "select": "id,doi,title,publication_year,abstract_inverted_index"})
        x = get(url)
        (RAW / f"openalex_citers_{name}_{page}.json").write_text(x or "{}", encoding="utf-8")
        d = json.loads(x or "{}")
        for w in d.get("results", []):
            out.append({"source": "OpenAlex-citers", "query": f"cites:{name}", "id": w.get("doi") or w.get("id"),
                        "title": w.get("title") or "", "year": str(w.get("publication_year") or ""),
                        "abstract": inv_abstract(w.get("abstract_inverted_index"))})
        cursor = d.get("meta", {}).get("next_cursor")
        page += 1
        time.sleep(1.0)
    return out


def s2_search(qid, q):
    url = "https://api.semanticscholar.org/graph/v1/paper/search?" + urllib.parse.urlencode(
        {"query": q, "limit": 50, "fields": "title,year,abstract,externalIds"})
    x = get(url, tries=1, wait=0.0)
    (RAW / f"s2_{qid}.json").write_text(x or "{}", encoding="utf-8")
    try:
        d = json.loads(x or "{}")
    except json.JSONDecodeError:
        d = {}
    out = []
    for p in d.get("data", []) or []:
        ext = p.get("externalIds") or {}
        pid = ("arXiv:" + ext["ArXiv"]) if ext.get("ArXiv") else (("https://doi.org/" + ext["DOI"]) if ext.get("DOI") else p.get("paperId"))
        out.append({"source": "SemanticScholar", "query": qid, "id": pid, "title": p.get("title") or "",
                    "year": str(p.get("year") or ""), "abstract": p.get("abstract") or ""})
    time.sleep(3.5)
    return out


def main():
    hits, status = [], {}
    for qid, (aq, nq) in QUERIES.items():
        a = arxiv(qid, aq); o = openalex_search(qid, nq); s = s2_search(qid, nq)
        status[qid] = {"arxiv": len(a), "openalex": len(o), "s2": len(s), "arxiv_query": aq, "text_query": nq}
        hits += a + o + s
        print(qid, status[qid]["arxiv"], status[qid]["openalex"], status[qid]["s2"], flush=True)
    for qid, aq in SUPPLEMENTAL.items():
        a = arxiv(qid, aq)
        status[qid] = {"arxiv": len(a), "arxiv_query": aq}
        hits += a
        print(qid, len(a), flush=True)
    for name, wid in ANCHORS.items():
        c = openalex_citers(name, wid)
        status["cites:" + name] = {"openalex_citers": len(c), "openalex_id": wid}
        hits += c
        print("cites", name, len(c), flush=True)
    (HERE / "search_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    with (HERE / "hits.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["source", "query", "id", "title", "year", "abstract"])
        w.writeheader(); w.writerows(hits)
    print("total hits", len(hits))


if __name__ == "__main__":
    main()
