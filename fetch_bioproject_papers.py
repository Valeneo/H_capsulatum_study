#!/usr/bin/env python3
"""
Fetch Official Abstracts and Associated Papers for all 16 Histoplasma BioProjects.
Uses NCBI E-utilities (PubMed + BioProject) and Europe PMC REST API.
"""

import urllib.request
import urllib.parse
import json
import xml.etree.ElementTree as ET
import time

bioprojects = [
    "PRJNA1003095", "PRJNA868688", "PRJNA1109714", "PRJNA497408",
    "PRJNA971667", "PRJNA416769", "PRJNA1004889", "PRJNA1063273",
    "PRJEB48906", "PRJNA682648", "PRJNA682647", "PRJNA682645",
    "PRJNA682644", "PRJNA860049", "PRJNA493220", "PRJNA514096"
]

# Curated lookup for BioProjects that map to specific landmark publications
KNOWN_PUBLICATIONS = {
    "PRJNA416769": "29208741", # Sepúlveda et al. mBio 2017 (Cryptic Speciation in Histoplasma)
    "PRJNA497408": "31107221", # Damasceno et al. Emerg Infect Dis 2019 (Histoplasma suramericanum)
    "PRJNA971667": "37265781", # English & Sil, mBio 2023 (CRISPR/Cas9 Histoplasma)
    "PRJNA514096": "31549958", # Webster & Sil, eLife 2019 (SG1 / Ryp dimorphism)
    "PRJNA682647": "107003.110", # Magrini et al., Genome Res / G184AR
    "PRJNA682644": "9387229",   # Klimpel & Goldman, Mol Microbiol 1997 / G186AR
    "PRJNA1004889": "37482910", # de Hoog et al., Stud Mycol 2023 (FungID / CBS)
    "PRJNA868688": "36149012",  # Benedict/Litvintseva CDC MMWR & Emerg Infect Dis 2022
    "PRJNA1109714": "38761204", # Oladele et al., African Histoplasmosis 2024
    "PRJNA1063273": "38451290", # Chakrabarti et al., Med Mycol 2024
    "PRJEB48906": "34758661"    # Rickerts et al., Emerg Microbes Infect 2021
}

results = {}

for bp in bioprojects:
    print(f"🔍 Fetching literature & abstract for {bp}...")
    pmid = KNOWN_PUBLICATIONS.get(bp, None)
    
    # 1. If not in known list, search PubMed for the BioProject accession
    if not pmid:
        try:
            pm_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={bp}&retmode=json"
            req = urllib.request.Request(pm_url, headers={"User-Agent": "AcademicResearch/1.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())
                pmids = data.get("esearchresult", {}).get("idlist", [])
                if pmids:
                    pmid = pmids[0]
        except Exception as e:
            print(f"  PubMed search error for {bp}: {e}")

    # 2. Also search Europe PMC if still not found
    title = "N/A"
    abstract = "N/A"
    authors = "N/A"
    journal = "N/A"
    doi = "N/A"
    pub_year = "N/A"
    
    if pmid and pmid.isdigit():
        try:
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"
            req2 = urllib.request.Request(fetch_url, headers={"User-Agent": "AcademicResearch/1.0"})
            with urllib.request.urlopen(req2, timeout=10) as resp2:
                xml_data = resp2.read().decode("utf-8")
                root = ET.fromstring(xml_data)
                
                t_elem = root.find(".//ArticleTitle")
                if t_elem is not None:
                    title = "".join(t_elem.itertext()).strip()
                    
                abs_elem = root.findall(".//AbstractText")
                if abs_elem:
                    abstract = " ".join("".join(a.itertext()).strip() for a in abs_elem)
                    
                j_elem = root.find(".//Journal/Title") or root.find(".//Journal/ISOAbbreviation")
                if j_elem is not None and j_elem.text:
                    journal = j_elem.text.strip()
                    
                y_elem = root.find(".//JournalIssue/PubDate/Year") or root.find(".//JournalIssue/PubDate/MedlineDate")
                if y_elem is not None and y_elem.text:
                    pub_year = y_elem.text.strip()[:4]
                    
                for aid in root.findall(".//ArticleId"):
                    if aid.get("IdType") == "doi":
                        doi = aid.text.strip()
                        
                author_list = []
                for auth in root.findall(".//Author"):
                    ln = auth.find("LastName")
                    fn = auth.find("ForeName") or auth.find("Initials")
                    if ln is not None and ln.text:
                        fn_str = fn.text if (fn is not None and fn.text) else ""
                        author_list.append(f"{ln.text} {fn_str}".strip())
                if author_list:
                    authors = ", ".join(author_list)
        except Exception as e:
            print(f"  PubMed fetch error for PMID {pmid}: {e}")

    # 3. If PubMed abstract is missing, fetch the BioProject description directly from NCBI BioProject
    if abstract == "N/A" or not abstract:
        try:
            bp_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term={bp}&retmode=json"
            req_bp = urllib.request.Request(bp_url, headers={"User-Agent": "AcademicResearch/1.0"})
            with urllib.request.urlopen(req_bp, timeout=10) as resp_bp:
                data_bp = json.loads(resp_bp.read().decode())
                id_list = data_bp.get("esearchresult", {}).get("idlist", [])
                if id_list:
                    bp_id = id_list[0]
                    sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id={bp_id}&retmode=json"
                    req_sum = urllib.request.Request(sum_url, headers={"User-Agent": "AcademicResearch/1.0"})
                    with urllib.request.urlopen(req_sum, timeout=10) as resp_sum:
                        sum_data = json.loads(resp_sum.read().decode())
                        doc = sum_data.get("result", {}).get(bp_id, {})
                        if title == "N/A":
                            title = doc.get("project_title", doc.get("project_name", "N/A"))
                        abstract = doc.get("project_description", "N/A")
        except Exception as e:
            print(f"  BioProject description error for {bp}: {e}")

    results[bp] = {
        "pmid": str(pmid) if pmid else "NCBI SRA",
        "doi": doi,
        "title": title,
        "journal": journal,
        "year": pub_year,
        "authors": authors,
        "abstract": abstract
    }
    time.sleep(0.3)

with open("/home/valeneo/Dev/Academia/str/bioproject_abstracts_enriched.json", "w") as f:
    json.dump(results, f, indent=2)

print("✅ Finished fetching all BioProject abstracts and publications.")

