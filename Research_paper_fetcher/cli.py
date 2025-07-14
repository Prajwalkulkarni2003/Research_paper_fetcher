# pubmed_company_searcher/cli.py

import argparse
from Research_paper_fetcher.pubmed import search_pubmed, fetch_pubmed_details, results_to_csv

def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers with at least one non-academic author.")  
    parser.add_argument("query", type=str, help="PubMed query string (in quotes)")  ## search term (required)
    parser.add_argument("-f", "--file", type=str, help="CSV output file", default=None)   ## file name to save the results 
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")  ## --debug flag to see extra info

    args = parser.parse_args()

    try:
        if args.debug:
            print(f"Running query: {args.query}")
        ids = search_pubmed(args.query)  ##This searches PubMed with the query the user gave.
        if args.debug:
            print(f"Found {len(ids)} papers")

        data = fetch_pubmed_details(ids)  ##This gets author information and filters the company-affiliated ones.
        if args.file:
            results_to_csv(data, args.file)   ##This saves the final data into a CSV file if the user gave a file name.
            print(f"Saved {len(data)} records to {args.file}")
        if not data:
            print("⚠️  No matching company-affiliated papers found.")
    
        else:
            for record in data:
                print(record)

    except Exception as e:
        print(f"[ERROR] {e}")
