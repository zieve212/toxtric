"""
test_literature_client.py
Tests for the PubMed literature client.

These never touch the network: requests.get is replaced with fakes
returning canned esearch JSON and efetch XML. We verify our own search
parsing, XML abstract extraction, and result assembly.
Run with: pytest
"""

import pytest
import requests

from toxreadout import literature_client


class FakeResponse:
    def __init__(self, status_code=200, json_data=None, text=""):
        self.status_code = status_code
        self._json = json_data
        self.text = text

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


# A minimal but realistic efetch XML payload with two articles; the first
# has a two-part (labelled) abstract, the second has none.
EFETCH_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>11111111</PMID>
      <Article>
        <Journal><Title>Journal of Caffeine Research</Title>
          <JournalIssue><PubDate><Year>2021</Year></PubDate></JournalIssue>
        </Journal>
        <ArticleTitle>Caffeine antagonism at the adenosine A2A receptor</ArticleTitle>
        <Abstract>
          <AbstractText Label="BACKGROUND">Caffeine blocks adenosine receptors.</AbstractText>
          <AbstractText Label="RESULTS">It acts as an A2A antagonist in the brain.</AbstractText>
        </Abstract>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
  <PubmedArticle>
    <MedlineCitation>
      <PMID>22222222</PMID>
      <Article>
        <Journal><Title>Neuro Letters</Title>
          <JournalIssue><PubDate><MedlineDate>2019 Jan-Feb</MedlineDate></PubDate></JournalIssue>
        </Journal>
        <ArticleTitle>A short letter with no abstract</ArticleTitle>
      </Article>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>"""


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(literature_client.time, "sleep", lambda *_: None)


def test_pubmed_url():
    assert literature_client.pubmed_url("123") == "https://pubmed.ncbi.nlm.nih.gov/123/"


def test_search_pubmed_parses_count_and_ids(monkeypatch):
    payload = {"esearchresult": {"count": "47", "idlist": ["11111111", "22222222"]}}
    monkeypatch.setattr(
        literature_client.requests, "get", lambda *a, **k: FakeResponse(200, payload)
    )
    result = literature_client.search_pubmed("caffeine AND adenosine")
    assert result["total_results"] == 47
    assert result["pmids"] == ["11111111", "22222222"]


def test_fetch_articles_parses_xml(monkeypatch):
    monkeypatch.setattr(
        literature_client.requests,
        "get",
        lambda *a, **k: FakeResponse(200, text=EFETCH_XML),
    )
    articles = literature_client.fetch_articles(["11111111", "22222222"])
    assert len(articles) == 2

    first = articles[0]
    assert first["pmid"] == "11111111"
    assert "A2A receptor" in first["title"]
    assert first["journal"] == "Journal of Caffeine Research"
    assert first["year"] == "2021"
    # Labelled abstract sections should be joined into one paragraph.
    assert "Caffeine blocks adenosine receptors." in first["abstract"]
    assert "A2A antagonist" in first["abstract"]
    assert first["url"].endswith("/11111111/")

    # Second article uses a free-text date and has no abstract.
    assert articles[1]["year"] == "2019"
    assert articles[1]["abstract"] == ""


def test_fetch_articles_empty_list_returns_empty(monkeypatch):
    # Should not even make a request for an empty PMID list.
    monkeypatch.setattr(
        literature_client.requests,
        "get",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("should not be called")),
    )
    assert literature_client.fetch_articles([]) == []


def test_literature_lookup_assembles_result(monkeypatch):
    def fake_get(url, *a, **k):
        if "esearch" in url:
            return FakeResponse(
                200, {"esearchresult": {"count": "47", "idlist": ["11111111", "22222222"]}}
            )
        return FakeResponse(200, text=EFETCH_XML)

    monkeypatch.setattr(literature_client.requests, "get", fake_get)

    result = literature_client.literature_lookup("caffeine", "adenosine receptor A2A")
    # The query is now focused with mechanism terms.
    assert result["query"].startswith("caffeine AND adenosine receptor A2A")
    assert "binding" in result["query"]
    assert result["total_results"] == 47
    # The discussion should come from the on-target article with an abstract.
    assert "Caffeine blocks adenosine receptors." in result["discussion"]
    assert result["discussion_source"]["pmid"] == "11111111"
    assert len(result["top_references"]) == 2


def test_select_discussion_prefers_on_target():
    """Given several abstracts, the most protein-relevant one should win."""
    articles = [
        {"pmid": "1", "title": "Caffeine and exercise performance",
         "abstract": "A broad review of caffeine and sports.", "url": "u1"},
        {"pmid": "2", "title": "Adenosine A2A receptor antagonist binding by caffeine",
         "abstract": "Caffeine binds the adenosine A2A receptor.", "url": "u2"},
    ]
    best = literature_client._select_discussion(articles, "adenosine A2A receptor")
    assert best["pmid"] == "2"
