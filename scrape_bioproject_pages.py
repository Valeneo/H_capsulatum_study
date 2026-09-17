#!/usr/bin/env python3
"""
Scrape exact Titles, Descriptions, Submitting Organizations, Data Types,
and Linked Publications from https://www.ncbi.nlm.nih.gov/bioproject/<id>
for all 16 Histoplasma BioProjects.
"""

import urllib.request
import urllib.parse
import re
from bs4 import BeautifulSoup
import json
import time

bioprojects = [
    "PRJNA1003095", "PRJNA1004889", "PRJNA1109714", "PRJNA1063273", "PRJEB48906",
    "PRJNA868688",  "PRJNA971667",  "PRJNA682648",  "PRJNA682647",  "PRJNA682645",
    "PRJNA682644",  "PRJNA860049",  "PRJNA416769",  "PRJNA497408",  "PRJNA493220",
    "PRJNA514096"
]

results = {}

for bp in bioprojects:
    url = f"https://www.ncbi.nlm.nih.gov/bioproject/{bp}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"})
    print(f"\n=======================================================")
    print(f"Scraping NCBI BioProject Page: {url}")
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            html = resp.read().decode("utf-8")
            soup = BeautifulSoup(html, "html.parser")
            
            # Find Title (in h1 or specific project header)
            h1 = soup.find("h1")
            h1_text = h1.get_text().strip() if h1 else "N/A"
            
            # Find Description
            desc_div = soup.find("div", {"id": "DescrDiv"}) or soup.find("div", {"class": "project-description"})
            desc_text = desc_div.get_text().strip() if desc_div else ""
            
            # Main Info table
            table_info = {}
            for row in soup.find_all("tr"):
                th = row.find("th")
                td = row.find("td")
                if th and td:
                    k = th.get_text().strip().replace(":", "")
                    v = td.get_text().strip()
                    table_info[k] = v
                    
            # Find all linked publication titles & PMIDs/PMCs
            pubs = []
            for a in soup.find_all("a", href=re.compile(r"/pubmed/\d+|/pmc/articles/PMC\d+")):
                p_text = a.get_text().strip()
                p_href = a.get("href")
                if p_text and p_href not in [x[1] for x in pubs]:
                    pubs.append({"text": p_text, "link": p_href})
                    
            print(f"📌 [{bp}] NCBI Page Title: {h1_text}")
            if "Organism" in table_info:
                print(f"   Organism: {table_info['Organism']}")
            if "Project data type" in table_info:
                print(f"   Data Type: {table_info['Project data type']}")
            if "Organization" in table_info or "Submitter" in table_info:
                sub = table_info.get("Organization", table_info.get("Submitter", "N/A"))
                print(f"   Submitter: {sub}")
            if desc_text:
                print(f"   Description: {desc_text[:120]}...")
            if pubs:
                print(f"   Publications on NCBI Page: {len(pubs)} link(s)")
                for p in pubs:
                    print(f"     - {p['text']} ({p['link']})")
                    
            results[bp] = {
                "bioproject_id": bp,
                "ncbi_url": url,
                "page_title": h1_text,
                "table_metadata": table_info,
                "description": desc_text,
                "linked_publications": pubs
            }
    except Exception as e:
        print(f"❌ Error scraping {bp}: {e}")
    time.sleep(0.5)

with open("/home/valeneo/Dev/Academia/str/ncbi_bioproject_scraped_data.json", "w") as f:
    json.dump(results, f, indent=2)

print("\n✅ Scraping finished successfully. Saved to ncbi_bioproject_scraped_data.json")
