#! /usr/bin/env python

import argparse
import os
import sys
import requests
import zipfile
from time import sleep
import re

ENTREZ_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_DATASETS_REST = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id"

EMAIL = os.environ.get("NCBI_EMAIL")
API_KEY = os.environ.get("NCBI_API_KEY")

def search_gene_symbol(symbol, organism="Homo sapiens"):
    """Use Entrez ESearch to resolve a gene symbol to a Gene ID."""
    params = {
        "db": "gene",
        "term": f"{symbol}[Gene Name] AND {organism}[Organism]",
        "retmode": "json",
        "email": EMAIL,
        "api_key": API_KEY
    }
    r = requests.get(ENTREZ_ESEARCH, params=params)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None

def get_gene_symbol(gene_id):
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params=params)
    r.raise_for_status()
    result = r.json().get("result", {})
    gene_info = result.get(gene_id, {})
    return gene_info.get("name")

def download_gene_sequences(gene_id, output_zip, protein=False):
    """Download gene data from NCBI Datasets REST API."""

    headers = {}
    basic_params = {}
    basic_params["api_key"] = API_KEY

    if not protein:
        basic_params["include_annotation_type"] = "FASTA_GENE"
    else:
        basic_params["include_annotation_type"] = "FASTA_PROTEIN"

    url = f"{NCBI_DATASETS_REST}/" + gene_id + "/download"

    r = requests.get(url,
                     params=basic_params,
                     headers=headers)
    r.raise_for_status()

    with open(output_zip, "wb") as f:
        f.write(r.content)

    print(f"Downloaded: {output_zip}")

def extract_file_from_zip(zip_path, out_dir):
    """Extract files with .fa, .faa, .fna, or .fasta extensions using regex."""
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        ext_pattern = re.compile(r"\.f(ast|n|a|)?a$", re.IGNORECASE)

        for member in zip_ref.namelist():
            if ext_pattern.search(member):
                filename = os.path.basename(member)
                if not filename:
                    continue

                with zip_ref.open(member) as source, open(os.path.join(out_dir, filename), "wb") as fp:
                    fp.write(source.read())

                print(f"Extracted {filename} to {out_dir}/")

def check_zip_extension(filename):
    if not filename.lower().endswith('.zip'):
        raise argparse.ArgumentTypeError("Output filename must end with .zip")
    return filename

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genes", nargs="?", required=True, help="List of gene symbols (e.g., BRCA1 TP53)")
    p.add_argument("--organism", nargs="?")
    p.add_argument("--protein", action="store_true")
    p.add_argument("-o", "--output", type=check_zip_extension, default="ncbi_genes.zip", help="Output ZIP filename (Must end in '.zip')")
    p.add_argument("--extract", action="store_true", help="Extract sequence files from downloaded ZIP")
    args = p.parse_args()

    if not EMAIL or not API_KEY:
        sys.exit("You must set NCBI_EMAIL and NCBI_API_KEY environment variables.")

    print("Resolving gene symbols to Gene IDs...")
    gene_id = search_gene_symbol(args.genes, args.organism)
    if gene_id:
        symbol = args.genes
        print(f"✔ {symbol} → Gene ID {gene_id}")

        download_gene_sequences(gene_id, args.output, args.protein)
    else:
        gene_id = args.genes

        download_gene_sequences(gene_id, args.output, args.protein)

        symbol = get_gene_symbol(args.genes)
        print(f"✔ {gene_id} → Gene symbol {symbol}")

    sleep(0.34)  # be polite, only 3 requests per second

    if args.extract:
        extract_dir = symbol + '_' + gene_id + "_data"
        os.makedirs(extract_dir, exist_ok=True)
        extract_file_from_zip(args.output, extract_dir)

if __name__ == "__main__":
    main()
