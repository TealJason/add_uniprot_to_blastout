#!bin/python3

import argparse
from pathlib import Path  
import subprocess
import requests 
import pandas as pd

def blastp(query, proteome):
    # Ensure output directory exists
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    # If query is a directory, process all .fasta or .fa files in it
    if query.is_dir():
        for fasta in query.glob("*.fasta") + query.glob("*.fa"):
            output_file = output_dir / f"{fasta.stem}_blastp_full_result.txt"
            subprocess.run(f"blastp -query {fasta} -db {proteome} -outfmt 6 -out {output_file}", shell=True, check=True)
            output_tsv = extract_unique_refseq_hits(output_file, proteome, output_dir)
            return output_tsv

    # If query is a single file, process it directly
    elif query.suffix in {'.fasta', '.fa'}:
        output_file = output_dir / f"{query.stem}_blastp_result.txt"
        subprocess.run(f"blastp -query {query} -db {proteome} -outfmt 6 -out {output_file}", shell=True, check=True)
        output_tsv = extract_unique_refseq_hits(output_file, proteome, output_dir)
        return output_tsv

def extract_unique_refseq_hits(output_file, proteome, output_dir):
    """Extract the best (lowest E-value) RefSeq match for each subject from the BLAST results."""
    prefix = proteome.stem
    output_tsv = output_dir / f"{prefix}_Refseq_best_hits.tsv"
    pd.options.display.float_format = '{:.8e}'.format

    try:
        df = pd.read_csv(output_file, sep='\t', header=None, 
                         names=['Query_ID', 'Subject_ID', 'pident','length','mismatch','gap_open',
                                'qstart','qend','sstart','send', 'E-value', 'BitScore'])

        # Convert 'E-value' to float to ensure correct sorting
        df['E-value'] = pd.to_numeric(df['E-value'], errors='coerce')

        # Get the index of the row with the smallest E-value for each Subject_ID
        index_of_minimum_values = df.groupby('Subject_ID')['E-value'].idxmin()

        # Extract those rows
        best_hits = df.loc[index_of_minimum_values].sort_values(by='E-value', ascending=True)

        # Save the filtered DataFrame
        best_hits.to_csv(output_tsv, sep='\t', index=False, columns=['Query_ID', 'Subject_ID','pident','length','mismatch','gap_open',
        'qstart','qend','sstart','send', 'E-value', 'BitScore'],float_format="%.2e")

        print(f"Best unique RefSeq matches saved in {output_tsv}")

    except FileNotFoundError:
        print(f"Error: BLAST output file {output_file} not found. Make sure BLASTP ran successfully.")

    return output_tsv
 
def map_refseq_to_uniprot(refseq_file):
    """Cross-references RefSeq IDs from a file with UniProt IDs and saves the mapping."""
    output_path = "output/Refseq_to_Uniprot_mapping.tsv"

    with open(refseq_file, 'r') as refseq_file, open(output_path, 'w') as refseq_mapping_file:
        refseq_mapping_file.write("RefSeq_ID\tUniProt_ID\tUniProt_Name\n")  # Writing the header
    
        for line in refseq_file:
            # Skip the header line if present
            if line.startswith("Query_ID"):
                continue
            columns = line.strip().split('\t')
            refseq_id = columns[1]  # Assuming the second column is the RefSeq ID

            uniprot_id,uniprot_name = get_uniprot_id(refseq_id)
            if uniprot_id:
                refseq_mapping_file.write(f"{refseq_id}\t{uniprot_id}\t{uniprot_name}\n")
            else:
                print(f"No UniProt match for RefSeq ID {refseq_id}")
    return output_path

def get_uniprot_id(refseq_id):
    """Fetches the UniProt ID for a given RefSeq ID using UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/stream?query=xref:{refseq_id}&format=json&columns=primaryAccession"
    print(f"Requesting URL: {url}")  # Debugging line to show the exact URL being requested

    try:
        response = requests.get(url)
        
        if response.status_code == 200:
        
            data = response.json()  # Try parsing as JSON
            
            # Check if the "results" key exists in the response
            results = data.get('results', [])
            if not results:
                print(f"No results returned for RefSeq ID {refseq_id}")
                return "N/A", "N/A"

            # Extract the primaryAccession from the first result
            uniprot_id = results[0].get('primaryAccession')
            uniprot_name = results[0].get('uniProtkbId')

            if uniprot_id:
                return uniprot_id, uniprot_name

            print(f"No matching UniProt entry for RefSeq ID {refseq_id}")
            return "N/A", "N/A"
        else:
            print(f"Error with request for RefSeq ID {refseq_id}: {response.status_code}")
            return None

    except Exception as e:
        print(f"Error fetching UniProt ID for {refseq_id}: {e}")
        return None

def add_uniprot_mapping_to_blast_results(refseq_id_file, refseq_mapping_file):
    """Adds UniProt ID and Name to the BLAST results."""
    output_file = 'output/modified_blast_results.tsv'
    uniprot_mapping = {}

    # Read the RefSeq-to-UniProt mapping file and store it in a dictionary
    with open(refseq_mapping_file, 'r') as mapping_file:
        for line in mapping_file:
            if line.startswith("RefSeq_ID"):
                continue  # Skip header line
            columns = line.strip().split('\t')
            refseq_id = columns[0]
            uniprot_id = columns[1]
            uniprot_name = columns[2]
            uniprot_mapping[refseq_id] = (uniprot_id, uniprot_name)

    # Open the BLAST file and prepare to modify it
    with open(refseq_id_file, 'r') as blast_file, open(output_file, 'w') as out_file:
        for line in blast_file:
            if line.startswith("Query_ID"):
                # Modify the header to insert new column names in the correct positions
                columns = line.strip().split('\t')
                new_header = columns[:2] + ["UniProt_ID", "UniProt_Name"] + columns[2:]
                out_file.write("\t".join(new_header) + "\n")
                continue

            columns = line.strip().split('\t')
            refseq_id = columns[1]  # Assuming the second column is the RefSeq ID

            # Insert UniProt data at the 3rd and 4th position if available
            if refseq_id in uniprot_mapping:
                uniprot_id, uniprot_name = uniprot_mapping[refseq_id]
                columns = columns[:2] + [uniprot_id, uniprot_name] + columns[2:]

            out_file.write("\t".join(columns) + "\n")  # Write modified line

    print(f"Modified BLAST results saved in {output_file}")
  
def args_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='This script runs BLASTP and extracts the best unique RefSeq match for each query sequence.')
    parser.add_argument('--query', required=True, type=Path, help='The query file or directory')
    parser.add_argument('--proteome', required=True, type=Path, help='The proteome file')

    return parser.parse_args()

def main():
    args = args_parser()
    refseq_id_file = blastp(args.query, args.proteome)
    refseq_mapping_file_path = map_refseq_to_uniprot(refseq_id_file)
    add_uniprot_mapping_to_blast_results(refseq_id_file,refseq_mapping_file_path)

if __name__ == "__main__":
    main()
