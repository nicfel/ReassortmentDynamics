import os
import sys
from Bio import SeqIO

# check fastas have corresponding taxa for reassortment analysis and make new fasta files with corresponding taxa


def get_fasta_headers(file_path):
    headers = set()
    for record in SeqIO.parse(file_path, "fasta"):
        headers.add(record.id)
    return headers

def filter_fasta_by_headers(input_file, output_file, common_headers):
    seen_headers = set()
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id in common_headers and record.id not in seen_headers:
                SeqIO.write(record, output_handle, "fasta")
                seen_headers.add(record.id)

def main(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.fasta')]
    
    if not fasta_files:
        print("No FASTA files found in the input directory.")
        return
    
    # Get headers from the first FASTA file
    first_fasta = os.path.join(input_dir, fasta_files[0])
    common_headers = get_fasta_headers(first_fasta)
    
    # Intersect headers with other FASTA files
    for fasta_file in fasta_files[1:]:
        file_path = os.path.join(input_dir, fasta_file)
        headers = get_fasta_headers(file_path)
        common_headers.intersection_update(headers)
    
    # Filter and write new FASTA files
    for fasta_file in fasta_files:
        input_file = os.path.join(input_dir, fasta_file)
        output_file = os.path.join(output_dir, fasta_file)
        filter_fasta_by_headers(input_file, output_file, common_headers)
    
    print(f"Filtered FASTA files have been saved to {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <output_directory>")
    else:
        input_directory = sys.argv[1]
        output_directory = sys.argv[2]
        main(input_directory, output_directory)