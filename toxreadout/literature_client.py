"""
literature_client.py
Searches PubMed for papers describing how a compound interacts with a
protein target, and pulls back real abstract text plus links.

Uses NCBI E-utilities: esearch to find matching papers (and the total
hit count), then efetch (XML) to retrieve titles, journals, years, and
abstracts. The top paper's abstract serves as a plain-English discussion
of the interaction; every top paper gets a clickable PubMed link.
"""

import re
import sys
import time
import xml.etree.ElementTree as ET

import requests

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PUBMED_BASE = "https://pubmed.ncbi.nlm.nih.gov"

# NCBI allows 3 requests/second without an API key.
MIN_INTERVAL = 0.34
MAX_RETRIES = 3
RETRY_STATUS = (500, 502, 503)

_last_request_time = 0.0


def _throttle() -> None:
    """Pause if needed so we stay under the rate limit."""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - elapsed)
    _last_request_time = time.time()


def _request(url: str, params: dict | None = None):
    """Rate-limited GET with retry/backoff on transient server errors."""
    last_response = None
    for attempt in range(MAX_RETRIES):
        _throttle()
        last_response = requests.get(url, params=params, timeout=30)
        if last_response.status_code in RETRY_STATUS:
            time.sleep(2 ** attempt)
            continue
        return last_response
    return last_response


def pubmed_url(pmid: str) -> str:
    """Build the public PubMed page URL for a PMID."""
    return f"{PUBMED_BASE}/{pmid}/"


def search_pubmed(query: str, retmax: int = 20) -> dict:
    """
    Search PubMed for a query. Returns the total hit count and the top
    PMIDs (most relevant first).
    """
    url = f"{NCBI_EUTILS}/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "sort": "relevance",
    }
    response = _request(url, params=params)
    response.raise_for_status()

    result = response.json().get("esearchresult", {})
    return {
        "total_results": int(result.get("count", 0)),
        "pmids": result.get("idlist", []),
    }


def _text(element) -> str:
    """Flatten an XML element (and any nested tags) into plain text."""
    if element is None:
        return ""
    return "".join(element.itertext()).strip()


def fetch_articles(pmids: list[str]) -> list[dict]:
    """
    Fetch full details for a list of PMIDs via efetch (XML).
    Returns one dict per article: pmid, title, journal, year, abstract, url.
    """
    if not pmids:
        return []

    url = f"{NCBI_EUTILS}/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    response = _request(url, params=params)
    response.raise_for_status()

    root = ET.fromstring(response.text)
    articles = []
    for article in root.findall(".//PubmedArticle"):
        pmid = _text(article.find(".//PMID"))
        title = _text(article.find(".//Article/ArticleTitle"))
        journal = _text(article.find(".//Article/Journal/Title"))

        year = _text(article.find(".//Article/Journal/JournalIssue/PubDate/Year"))
        if not year:
            # Some records use a free-text date like "2019 Jan".
            medline_date = _text(
                article.find(".//Article/Journal/JournalIssue/PubDate/MedlineDate")
            )
            year = medline_date[:4] if medline_date else ""

        # An abstract may be split into labelled sections; join them.
        parts = [_text(node) for node in article.findall(".//Abstract/AbstractText")]
        abstract = " ".join(p for p in parts if p)

        articles.append(
            {
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "year": year,
                "abstract": abstract,
                "url": pubmed_url(pmid),
            }
        )
    return articles


# Terms that signal a paper is actually about a compound-protein interaction
# (rather than a broad review). Used both to focus the query and to rank
# which paper's abstract becomes the discussion.
MECHANISM_TERMS = ["binding", "antagonist", "agonist", "inhibitor", "interaction"]

# Generic words to ignore when turning a protein name into match keywords.
_STOPWORDS = {"the", "and", "of", "a", "an", "for", "to"}


def _keywords(protein: str) -> list:
    """Reduce a protein name to meaningful lowercase match keywords."""
    tokens = re.findall(r"[a-z0-9]+", (protein or "").lower())
    return [t for t in tokens if len(t) >= 2 and t not in _STOPWORDS]


def _score_article(article: dict, keywords: list) -> int:
    """
    Rank how on-target an article is for the interaction: protein keywords
    in the title count most, then in the abstract, plus a bump for
    mechanism terms in the title.
    """
    title = (article.get("title") or "").lower()
    abstract = (article.get("abstract") or "").lower()
    score = 0
    for keyword in keywords:
        if keyword in title:
            score += 2
        if keyword in abstract:
            score += 1
    for term in MECHANISM_TERMS:
        if term in title:
            score += 1
    return score


def _select_discussion(articles: list, protein: str) -> dict | None:
    """Pick the most on-target article that actually has an abstract."""
    keywords = _keywords(protein)
    best = None
    best_score = -1
    for article in articles:
        if not article.get("abstract"):
            continue
        score = _score_article(article, keywords)
        if score > best_score:
            best_score = score
            best = article
    return best


def literature_lookup(compound: str, protein: str, top_n: int = 5) -> dict:
    """
    Search PubMed for papers on a compound-protein interaction and return
    a discussion (the most relevant paper's abstract) plus linked
    references.
    """
    # Focus the search on interaction/mechanism papers so the top hits are
    # about the pairing itself, not broad reviews that merely mention both.
    query = f"{compound} AND {protein} AND ({' OR '.join(MECHANISM_TERMS)})"
    search = search_pubmed(query)

    # Fetch a few extra so the on-target ranking has candidates to choose from.
    articles = fetch_articles(search["pmids"][: max(top_n, 1) + 5])

    # Use the most on-target paper (by protein/mechanism keywords) as the
    # discussion, rather than just the first one with an abstract.
    best = _select_discussion(articles, protein)
    discussion = best["abstract"] if best else ""
    discussion_source = (
        {"pmid": best["pmid"], "title": best["title"], "url": best["url"]}
        if best
        else None
    )

    references = [
        {
            "pmid": a["pmid"],
            "title": a["title"],
            "journal": a["journal"],
            "year": a["year"],
            "url": a["url"],
        }
        for a in articles[:top_n]
    ]

    return {
        "query": query,
        "total_results": search["total_results"],
        "discussion": discussion,
        "discussion_source": discussion_source,
        "top_references": references,
    }


# Runs only when you execute this file directly: a live caffeine /
# adenosine A2A receptor interaction lookup.
if __name__ == "__main__":
    # Scientific abstracts contain Unicode (e.g. thin spaces, Greek letters);
    # make sure the Windows console prints them instead of crashing.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    result = literature_lookup("caffeine", "adenosine receptor A2A")

    print(f"Query:          {result['query']}")
    print(f"Total results:  {result['total_results']:,} papers on PubMed\n")

    print("DISCUSSION OF THE INTERACTION")
    print("-" * 60)
    print(result["discussion"] or "(no abstract available)")

    source = result["discussion_source"]
    if source:
        print(f"\nSource: {source['title']}")
        print(f"        {source['url']}")

    print("\nMORE PAPERS ON THIS INTERACTION")
    print("-" * 60)
    for ref in result["top_references"]:
        print(f"[{ref['year']}] {ref['title']}")
        print(f"        {ref['journal']}")
        print(f"        {ref['url']}\n")
