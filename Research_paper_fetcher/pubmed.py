# pubmed_company_searcher/pubmed.py

from typing import List, Dict
from Bio import Entrez
import pandas as pd
import re

Entrez.email = "prajwal.ise.rymec@gmail.com"  

NON_ACADEMIC_KEYWORDS = ["Pharma", "Biotech", "Inc", "Ltd", "Corp", "Corporation"]
ACADEMIC_KEYWORDS = ["University", "Institute", "Hospital", "School", "College", "Lab"]

def search_pubmed(query: str, max_results: int = 20) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    results = Entrez.read(handle)
    return results["IdList"]

def fetch_pubmed_details(pubmed_ids: List[str]) -> List[Dict]:
    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), rettype="medline", retmode="text")
    records = handle.read()
    return parse_medline_records(records)

def parse_medline_records(records: str) -> List[Dict]:
    result = []
    papers = records.split("\n\n")
    for paper in papers:
        paper_data = {"Non-academic Author(s)": [], "Company Affiliation(s)": []}
        for line in paper.split("\n"):
            if line.startswith("PMID-"):
                paper_data["PubmedID"] = line.replace("PMID- ", "").strip()
            elif line.startswith("TI  -"):
                paper_data["Title"] = line.replace("TI  - ", "").strip()
            elif line.startswith("DP  -"):
                paper_data["Publication Date"] = line.replace("DP  - ", "").strip()
            elif line.startswith("AD  -"):
                affil = line.replace("AD  - ", "").strip()
                print("🧪 Found affiliation:", affil)
                if any(keyword in affil for keyword in NON_ACADEMIC_KEYWORDS) and not any(
                        keyword in affil for keyword in ACADEMIC_KEYWORDS):
                    paper_data["Non-academic Author(s)"].append(affil)
                    paper_data["Company Affiliation(s)"].append(affil)
                if "@" in affil:
                    paper_data["Corresponding Author Email"] = re.findall(r"\S+@\S+", affil)[0]
        if paper_data.get("Company Affiliation(s)"):
            result.append(paper_data)
    return result

def results_to_csv(results: List[Dict], filename: str):
    df = pd.DataFrame(results)
    df.to_csv(filename, index=False)
