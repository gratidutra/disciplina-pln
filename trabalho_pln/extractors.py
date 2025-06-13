import io
import json
import os

import pandas as pd
import pubmed_parser as pp
import requests
from Bio import Entrez
from dotenv import load_dotenv
from metapub import PubMedFetcher
from metapub.findit import FindIt
from tqdm import tqdm


def pmids_extractor(keyword, number_articles=9999, start_date="1950-01-01", end_date=pd.to_datetime('today')):
    start_date = "1950-01-01"
    end_date = pd.to_datetime('today')
    dates = pd.date_range(start=start_date, end=end_date, freq='MS')

    chunk_size = 9999
    cursor = 0
    pmids_list = []  # Esta será a lista consolidada de todos os PMIDs

    for date in tqdm(dates, desc="Extraindo do PubMed"): 
        #print(date)
        start = date.strftime("%Y/%m/%d")
        end = (date + pd.DateOffset(months=1) - pd.DateOffset(days=1)).strftime("%Y/%m/%d")
        if date.month == end_date.month and date.year == end_date.year:
            end = end_date.strftime("%Y/%m/%d")
        #AND (("2020/01/05"[Date - Publication] : "3000"[Date - Publication]))
        search_term = f'{keyword} AND (("{start}"[Date - Publication] : "{end}"[Date - Publication]))'
        try:
            pmids = Entrez.esearch(db='pubmed', term=search_term, retmax=chunk_size)
            pmid_list = Entrez.read(pmids)['IdList']
            if pmid_list:  # Só adiciona se a lista não estiver vazia
                pmids_list.extend(pmid_list)  # Usamos extend() em vez de append()
        except Exception as e:
            print(f"Erro na consulta para {start} a {end}: {e}")
    return pmids_list
 

def findIt_extractor(pmid, _):
    try:
        finder = FindIt(pmid, verify=False)
        return finder.url if finder.url else "Dont Found"
    except Exception as e:
        return f"PMID Error"
        
def pubmed_extractor(pmids_list, output_path="data/pubmed_output.jsonl"):
    fetch = PubMedFetcher()
    not_found_pmids = []

    with open(output_path, "a", encoding="utf-8") as f_out:
        for pmid in tqdm(pmids_list, desc="Extraindo do PubMed"):
            try:
                article = fetch.article_by_pmid(pmid)
                xml = article.xml
                parsed = pp.parse_medline_xml(
                    xml,
                    year_info_only=False,
                    nlm_category=False,
                    author_list=False,
                    reference_list=False,
                )
                for entry in parsed:
                    entry["url"] = findIt_extractor(pmid, _)
                    f_out.write(json.dumps(entry) + "\n")
            except Exception as e:
                tqdm.write(f"PMID {pmid}: {e}")
                not_found_pmids.append(pmid)
                continue

    # Salva os PMIDs com erro
    if not_found_pmids:
        with open("data/pmids_nao_encontrados.json", "w", encoding="utf-8") as f:
            json.dump(not_found_pmids, f, indent=2)

    print(f"Extração finalizada. Total com erro: {len(not_found_pmids)}")

#def fullText_extractor(url):
#    try:
#        response = requests.get(url, timeout=10)
#        response.raise_for_status()
#
#        with pdfplumber.open(io.BytesIO(response.content)) as pdf:
#            text = ""
#            for page in pdf.pages:
#                text += page.extract_text() or ""
#
#            return {"url": url, "content": text}
#    except requests.RequestException as req_err:
#        return {"url": url, "error": f"Erro ao baixar PDF: {req_err}"}
#    except Exception as e:
#        return {"url": url, "error": f"Erro ao processar PDF: {e}"}
