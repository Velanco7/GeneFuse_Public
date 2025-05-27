import argparse
import subprocess
import os
import glob
import csv
import time
import datetime
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio import Entrez, SeqIO
from time import sleep


# Set your NCBI API credentials
#!---------------------------------------------------------
Entrez.email = "vf2g20@soton.ac.uk"
Entrez.api_key = "5c2771fab336c38f29ea0de197cf7f5edc09"
#!---------------------------------------------------------
BATCH_SIZE = 1000
DELAY = 0.34

def read_protein_ids(filepath, filter_file=None):
    """Extract Protein IDs from the input file, optionally filtering against another file."""
    protein_ids = set()
    filter_ids = set()
    
    # Read filter IDs and append .1 if needed
    if filter_file:
        print("Reading filter file...")
        with open(filter_file, 'r') as f:
            filter_ids = {line.strip() + ".1" for line in f if line.strip()}
        print(f"Total protein IDs in filter file (with '.1' appended): {len(filter_ids)}")
    
    # Read main protein IDs
    with open(filepath, 'r') as file:
        next(file)  # Skip header
        for line in file:
            if line.strip() == "":
                continue
            parts = line.strip().split('\t')
            if parts and parts[0] not in filter_ids:
                protein_ids.add(parts[0])
    return list(protein_ids)

def load_common_hmms(csv_path):
    common_hmms = set()
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            common_hmms.add(row["Accession"].strip())
    return common_hmms

def extract_accessions_from_summary(summary_file, common_hmms):
    accessions = set()
    with open(summary_file) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            hmm_accession = parts[1].strip()
            protein_accession = parts[2].strip()
            if hmm_accession in common_hmms:
                accessions.add(protein_accession.split('.')[0])  # Strip version
    return accessions

def fetch_sequences(protein_ids, output_fasta):
    """Fetch sequences in batches via Entrez."""
    with open(output_fasta, 'w') as fasta_out:
        for i in range(0, len(protein_ids), BATCH_SIZE):
            batch = protein_ids[i:i + BATCH_SIZE]
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    handle = Entrez.efetch(
                        db='protein', id=','.join(batch),
                        rettype='fasta', retmode='text'
                    )
                    records = list(SeqIO.parse(handle, 'fasta'))
                    SeqIO.write(records, fasta_out, 'fasta')
                    print(f"✅ Batch {i//BATCH_SIZE + 1}: {len(records)} sequences fetched")
                    handle.close()
                    break  # Success, exit retry loop
                except Exception as e:
                    if attempt == max_retries - 1:
                        print(f"❌ Failed batch {i//BATCH_SIZE + 1} after {max_retries} attempts: {e}")
                    else:
                        print(f"⚠️ Retrying batch {i//BATCH_SIZE + 1} (attempt {attempt + 1})...")
                        sleep(DELAY * 2)  # Longer delay for retries
            sleep(DELAY)

def run_makeblastdb(fasta_file, db_name):
    """Build BLAST database using makeblastdb."""
    cmd = [
        'makeblastdb',
        '-in', fasta_file,
        '-dbtype', 'prot',
        '-out', db_name,
        '-parse_seqids'
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ BLAST database '{db_name}' created successfully!")
    except subprocess.CalledProcessError as e:
        print(f"❌ makeblastdb failed: {e}")

def run_blastp(query_fasta, db_name, output_csv, threads=4):
    """Run all-vs-all BLASTP with the original format 6 output plus headers."""
    # Original format 6 columns
    headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    cmd = [
        'blastp',
        '-query', query_fasta,
        '-db', db_name,
        '-out', output_csv,
        '-outfmt', '6',  # Original tabular format
        '-evalue', '1e-5',
        '-num_threads', str(threads),
        '-max_hsps', '1'
    ]
    try:
        subprocess.run(cmd, check=True)
        
        # Add headers to the CSV file
        with open(output_csv, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('\t'.join(headers) + '\n' + content)
        print(f"✅ BLASTP results saved to {output_csv}")
    except subprocess.CalledProcessError as e:
        print(f"❌ blastp failed: {e}")


def run_blastp_simplified(query_fasta, db_name, output_csv, threads=4):
    """Run BLASTP with simplified 3-column output format (no headers)"""
    cmd = [
        'blastp',
        '-query', query_fasta,
        '-db', db_name,
        '-out', output_csv,
        '-outfmt', '6 qseqid sseqid pident',  # Only these 3 columns
        '-evalue', '1e-5',
        '-num_threads', str(threads),
        '-max_hsps', '1'
    ]
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ BLASTP results saved to {output_csv}")
    except subprocess.CalledProcessError as e:
        print(f"❌ blastp failed: {e}")
        raise


def filter_high_identity_pairs(input_csv, output_csv, threshold=95):
    """
    Filter BLAST results for high identity matches.
    Handles both simplified (3-col) and full BLAST output formats.
    """
    best_matches = defaultdict(float)
    
    with open(input_csv, 'r') as infile:
        # Read as simple TSV (qseqid, sseqid, pident)
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                continue
                
            qseqid = row[0]
            sseqid = row[1]
            try:
                pident = float(row[2])
            except ValueError:
                continue
                
            # Skip self-matches
            if qseqid == sseqid:
                continue
                
            # Create sorted tuple key
            pair = tuple(sorted((qseqid, sseqid)))
            
            # Keep highest identity
            if pident > best_matches[pair]:
                best_matches[pair] = pident
    
    # Write filtered results with headers
    with open(output_csv, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['Protein1', 'Protein2', 'Percent_Identity'])
        
        for (prot1, prot2), pident in best_matches.items():
            if pident >= threshold:
                writer.writerow([prot1, prot2, round(pident, 3)])
    
    print(f"✅ Filtered high-identity pairs saved to {output_csv}")

def find_protein_clusters(input_file, output_file, min_identity=95):
    """
    Groups interconnected proteins into clusters based on similarity relationships.
    
    Args:
        input_file: Path to TSV file with protein pairs and identities
        output_file: Path to save clustered results
        min_identity: Minimum identity threshold for considering connections
    """
    # Build adjacency list of protein connections
    graph = defaultdict(set)
    
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if float(row['Percent_Identity']) >= min_identity:
                prot1 = row['Protein1']
                prot2 = row['Protein2']
                graph[prot1].add(prot2)
                graph[prot2].add(prot1)
    
    # Find connected components
    visited = set()
    clusters = []
    
    for protein in graph:
        if protein not in visited:
            # Start a new cluster
            cluster = set()
            stack = [protein]
            
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    cluster.add(current)
                    # Add all neighbors to stack
                    stack.extend(graph[current] - visited)
            
            if cluster:
                clusters.append(cluster)
    
    # Write results
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ClusterID', 'Proteins'])
        
        for i, cluster in enumerate(sorted(clusters, key=lambda x: -len(x)), 1):
            # Sort proteins alphabetically within cluster
            sorted_proteins = sorted(cluster)
            quoted_proteins = [f'{prot}' for prot in sorted_proteins]
            writer.writerow([f'Cluster_{i}', ', '.join(quoted_proteins)])
    
    print(f"✅ Found {len(clusters)} protein clusters saved to {output_file}")


def filter_fasta_by_clusters(cluster_csv, input_fasta, output_fasta):
    """
    Remove proteins that appear after the first protein in each cluster.
    
    Args:
        cluster_csv: Path to cluster CSV file (format: ClusterID, Proteins)
        input_fasta: Path to original FASTA file
        output_fasta: Path to save filtered FASTA
    """
    # Step 1: Collect proteins to remove
    proteins_to_remove = set()
    
    with open(cluster_csv, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            # Extract proteins from quoted list like: "WP_1", "WP_2", "WP_3"
            proteins = [p.strip('" ') for p in row[1].split(',')]
            # Add all except first protein in cluster
            proteins_to_remove.update(proteins[1:])
    
    print(f"Found {len(proteins_to_remove)} proteins to remove from clusters")
    
    # Step 2: Filter FASTA file
    kept_records = []
    removed_count = 0
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in proteins_to_remove:
            removed_count += 1
        else:
            kept_records.append(record)
    
    # Step 3: Save filtered FASTA
    with open(output_fasta, 'w') as f:
        SeqIO.write(kept_records, f, "fasta")
    
    print(f"Removed {removed_count} proteins")
    print(f"Kept {len(kept_records)} proteins in {output_fasta}")
    print(f"Total proteins to remove that weren't found: {len(proteins_to_remove) - removed_count}")


def run_hmmscan(input_fasta, output_file, hmm_databases):
    
    E_VALUE_THRESH = "0.001"
    DOMAIN_E_VALUE_THRESH = "0.001"

    """
    Runs hmmscan using the input FASTA file and the given HMM databases.
    Combines results from multiple databases into a single output file.
    """
    print("Running HMMER scans...")
    print(input_fasta)

    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

    try:
        # Temporary files to hold individual database results
        temp_files = []
        for hmm_database in hmm_databases:
            if not os.path.exists(hmm_database):
                print(f"Error: HMM database not found: {hmm_database}")
                continue
            
            print(f"Scanning with database: {hmm_database}")

            # Build the hmmscan command
            temp_file = f"{hmm_database.split('/')[-1]}_result.tblout"  # Use the name of the database for temp file
            temp_files.append(temp_file)

            command = [
                "hmmscan",
                "--domtblout", temp_file,  # Output table to a temporary file
                "-E", E_VALUE_THRESH,  # Sequence E-value threshold
                "--domE", DOMAIN_E_VALUE_THRESH,  # Domain E-value threshold
                hmm_database,
                input_fasta,
            ]

            # Run the command and capture stdout and stderr
            result = subprocess.run(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            
            if result.returncode != 0:
                print(f"Error in hmmscan with database {hmm_database}: {result.stderr}")
                continue

        # Combine the temporary files into the final output file
        with open(output_file, "w") as outfile:
            for temp_file in temp_files:
                # Add header for each database
                with open(temp_file, "r") as temp:
                    outfile.write(f"# Results for {temp_file}\n")
                    outfile.write(temp.read())
                    outfile.write("\n")  # Ensure separation between results

                # Remove the temporary file after processing
                os.remove(temp_file)

    except subprocess.CalledProcessError as e:
        print("Error running hmmscan:", e)
        raise
    except FileNotFoundError:
        print("Error: hmmscan not found. Ensure HMMER is installed and in PATH.")
        raise

    print(f"HMMER scan results saved to: {output_file}")

def parse_hmmscan_output(output_file):
    """
    Parses the hmmscan output file and extracts domain hits.
    """
    domain_hits = []

    with open(output_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue

            domain_hit = {
                "target_name": parts[0],
                "accession": parts[1],
                "query_name": parts[3],
                "E-value": parts[6],
                "score": parts[7],
                "start": parts[17],  # Start of domain
                "end": parts[18],    # End of domain
                "description": " ".join(parts[22:]),
            }
            domain_hits.append(domain_hit)

    return domain_hits

def write_domain_hits(domain_hits, output_summary_file):
    """
    Writes the domain hits to a human-readable summary file.
    """
    # Open the output file for writing
    with open(output_summary_file, "w") as f:
        # Write the header
        f.write("Target Name\tAccession\tQuery Name\tE-value\tScore\tStart\tEnd\tDescription\n")

        # Process each hit and validate its fields
        for hit in domain_hits:
            try:
                # Extract and clean individual fields
                target_name = hit.get("target_name", "").strip()
                accession = hit.get("accession", "").strip()
                query_name = hit.get("query_name", "").strip()
                e_value = hit.get("E-value", "").strip()
                score = hit.get("score", "").strip()
                start = hit.get("start", "").strip()
                end = hit.get("end", "").strip()
                description = hit.get("description", "").strip()

                # Write a properly formatted line to the file
                f.write(
                    f"{target_name}\t{accession}\t{query_name}\t{e_value}\t{score}\t"
                    f"{start}\t{end}\t{description}\n"
                )
            except Exception as e:
                # Handle unexpected errors and continue
                print(f"Error processing hit: {hit}. Skipping. Error: {e}")

    print(f"Domain summary saved to: {output_summary_file}")

def summarize_results(summary_file, output_cumulative_file, fasta_file):
    """
    Summarizes the results in the summary.txt file and writes a cumulative summary.
    """
    target_counts = Counter()

    # Count total sequences from the FASTA file
    total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    print(f"Total sequences in input FASTA file: {total_sequences}")

    # Read and parse the summary file
    with open(summary_file, "r") as f:
        for line_number, line in enumerate(f, 1):  # Enumerate for debugging line numbers
            if line.startswith("Target Name"):  # Skip header line
                continue
            parts = line.strip().split("\t")
            
            # Ensure the line contains all expected fields
            if len(parts) < 8:
                print(f"Warning: Skipping malformed line {line_number}: {line.strip()}")
                continue
            
            target_name = parts[0]
            accession = parts[1]
            description = parts[7]
            
            # Aggregate counts using a tuple of (target_name, accession, description)
            target_counts[(target_name, accession, description)] += 1

    # Write the cumulative summary
    with open(output_cumulative_file, "w") as f:
        f.write("Target Name\tAccession\tDescription\tOccurrences\n")
        for (target_name, accession, description), count in target_counts.items():
            f.write(f"{target_name}\t{accession}\t{description}\t{count}\n")
        
        f.write(f"\nTotal Sequences Processed: {total_sequences}\n")
    
    print(f"Cumulative summary created: {output_cumulative_file}")


def compare_hmm_accessions_with_description(psi_summary_path, neighbourhood_summary_path, output_dir):
    """
    Compare two summary files and print common HMM profile accessions with their descriptions.
    Assumes each HMM accession has one consistent description.
    """

    def extract_hmm_accession_description(filepath):
        accession_to_description = {}
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("Target Name") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    hmm_accession_raw = parts[1]  # PFxxxxx.xx
                    hmm_accession = hmm_accession_raw.split('.')[0]  # normalize version suffix
                    description = parts[7].strip()
                    if hmm_accession not in accession_to_description:
                        accession_to_description[hmm_accession] = description
        return accession_to_description

    psi_data = extract_hmm_accession_description(psi_summary_path)
    neighbourhood_data = extract_hmm_accession_description(neighbourhood_summary_path)

    common_hmm_accessions = set(psi_data.keys()) & set(neighbourhood_data.keys())
    output_file = os.path.join(output_dir, 'common_hmm_accessions.txt')

    if common_hmm_accessions:
        print("\nCommon HMM accessions and their descriptions:\n")
        with open(output_file, 'w') as f_out:
            f_out.write("Common HMM accessions and their descriptions:\n\n")
            f_out.write("Accession\tDescription\n")
            f_out.write("----------------------------\n")

            for acc in sorted(common_hmm_accessions):
                description = psi_data.get(acc) or neighbourhood_data.get(acc)
                line = f"{acc}\t{description}"
                print(line)
                f_out.write(line + "\n")
        
        print(f"\nResults saved to: {output_file}")
    else:
        print("\nNo common HMM accessions found.\n")


def delete_blastdb_files(db_name):
    """Delete files created by makeblastdb."""
    extensions = ['.pin', '.phr', '.psq', '.pdb', '.pjs', '.pog', '.pos', '.pto', '.pot', '.ptf']
    for ext in extensions:
        file_path = f"{db_name}{ext}"
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"🗑️ Deleted {file_path}")
        else:
            print(f"⚠️ File not found: {file_path}")

def load_accessions_with_counts(filepath):
    """
    Loads a file and returns a dict:
    {accession: (description, count)}
    """
    accession_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("Target Name") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                accession = parts[1].strip()
                description = parts[2].strip()
                try:
                    count = int(parts[3].strip())
                except ValueError:
                    count = 0
                accession_data[accession] = (description, count)
    return accession_data

def find_common_hmm_accessions_with_counts(file1, file2, output_dir=None):
    """
    Finds common HMM accessions and outputs their descriptions and occurrence counts.
    Sorted by the sum of occurrences in both files (descending).
    """
    data1 = load_accessions_with_counts(file1)
    data2 = load_accessions_with_counts(file2)

    common = set(data1.keys()) & set(data2.keys())
    if not common:
        print("No common HMM accessions found.")
        return

    combined_rows = []
    for acc in common:
        desc1, count1 = data1.get(acc, ("", 0))
        desc2, count2 = data2.get(acc, ("", 0))
        description = desc1 if desc1 else desc2
        total = count1 + count2
        combined_rows.append((acc, description, count2, count1, total))  # Swapped count1 and count2 here

    # Sort by total occurrences descending
    combined_rows.sort(key=lambda x: x[4], reverse=True)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, "common_hmm_accessions.csv")
        with open(out_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Accession", "Description", "Neighbour proteins", "Homolog proteins"])  # Fixed header order
            for row in combined_rows:
                writer.writerow(row[:4])  # Drop the total column in output
        print(f"\n✅ Saved sorted results to: {out_path}")
        return out_path

def common_hmm_blast_main(common_csv, psi_summary, neigh_summary, output_dir, identity_threshold, args):

    print("[1] Loading common HMMs...")
    common_hmms = load_common_hmms(common_csv)

    print("[2] Extracting protein accessions...")
    psi_prots = extract_accessions_from_summary(psi_summary, common_hmms)
    neigh_prots = extract_accessions_from_summary(neigh_summary, common_hmms)

    print(f"  - PSI proteins: {len(psi_prots)}")
    print(f"  - Neighbour proteins: {len(neigh_prots)}")

    psi_fasta = os.path.join(output_dir, "psi_proteins.fasta")
    neigh_fasta = os.path.join(output_dir, "neighbour_proteins.fasta")

    print("[3] Downloading FASTA sequences...")
    # Convert sets to lists here
    fetch_sequences(list(psi_prots), psi_fasta)
    fetch_sequences(list(neigh_prots), neigh_fasta)


    print("[4] Creating BLAST database...")
    db_name = os.path.join(output_dir, "neighbour_db")
    run_makeblastdb(neigh_fasta, db_name)

    print("[5] Running BLASTP...")
    blast_output = os.path.join(output_dir, "blast_results.txt")
    run_blastp_simplified(psi_fasta, db_name, blast_output, threads=args.threads)

    print("[6] Filtering high-identity BLAST hits...")
    high_id_output = os.path.join(output_dir, 'filtered_pairs.tsv')
    filter_high_identity_pairs(blast_output, high_id_output, identity_threshold)

    print("[7] Clustering proteins...")
    cluster_file = os.path.join(output_dir, 'protein_clusters.csv')
    find_protein_clusters(high_id_output, cluster_file, identity_threshold)

    print("[8] Cleaning up BLAST DB files...")
    delete_blastdb_files(db_name)

    print(f"\n[✔] Pipeline complete. Results saved to {output_dir}")

    return blast_output, psi_fasta, neigh_fasta

def extract_all_pairs(txt_file):
    """Extract and normalize protein pairs"""
    pairs = []
    with open(txt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 3:
                pid1 = normalize_blast_id(parts[0])
                pid2 = normalize_blast_id(parts[1])
                try:
                    pairs.append((pid1, pid2, float(parts[2])))
                except ValueError:
                    continue
    return pairs

def extract_ids_from_fasta(fasta_file):
    """Extract all protein IDs from a FASTA file."""
    return {record.id.strip() for record in SeqIO.parse(fasta_file, 'fasta')}

def parse_domain_summary(summary_file):
    """Parse domain summary file into a dictionary {base_protein_id: list_of_domain_names}."""
    domains = defaultdict(list)
    with open(summary_file, 'r') as f:
        next(f)  # Skip header line
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                protein_id = parts[2].strip()
                domain_name = parts[0].strip()
                base_id = protein_id.split('.')[0]
                domains[base_id].append(domain_name)
    return domains

def get_domains(protein_id, domain_map):
    """Get domains for a protein ID, handling version numbers."""
    base_id = protein_id.split('.')[0]
    return domain_map.get(base_id, [])

def format_domains(domain_list):
    """Format list of domains into a semicolon-separated string."""
    return "; ".join(domain_list) if domain_list else "NA"

def build_genome_id_map(csv_directory):
    """Build a dictionary mapping protein IDs to genome IDs from CSV files."""
    protein_to_genome = {}
    for csv_file in glob.glob(os.path.join(csv_directory, "*.csv")):
        with open(csv_file, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 6:
                    genome_id = row[0].strip()
                    protein_id = row[5].strip()
                    protein_to_genome[protein_id] = genome_id
    return protein_to_genome

def normalize_blast_id(protein_id):
    """Handle both simple and NCBI-style IDs"""
    if protein_id.startswith('ref|'):
        return protein_id.split('|')[1]  # Extract WP_... from ref|WP_...|
    return protein_id  # Already in simple format

def categorize_matches(pairs, fasta1_ids, fasta2_ids, neighbour_domains, psi_domains, protein_to_genome):
    
    """Modified matching with ID normalization"""
    # Ensure FASTA IDs are also normalized (if needed)
    fasta1_ids = {normalize_blast_id(pid) for pid in fasta1_ids}
    fasta2_ids = {normalize_blast_id(pid) for pid in fasta2_ids}

    """Categorize matches into different identity ranges."""
    categories = {
        '100_to_80': {'min': 80, 'max': 100, 'matches': []},
        '80_to_60': {'min': 60, 'max': 80, 'matches': []},
        '60_to_40': {'min': 40, 'max': 60, 'matches': []},
        'below_40': {'min': 0, 'max': 40, 'matches': []}
    }
    
    for pid1, pid2, identity in pairs:
        # Check both possible directions of matching
        if (pid1 in fasta1_ids and pid2 in fasta2_ids) or (pid1 in fasta2_ids and pid2 in fasta1_ids):
            # Determine which is the neighbor and which is the PSI protein
            if pid1 in fasta1_ids:
                neighbour = pid1
                psi = pid2
            else:
                neighbour = pid2
                psi = pid1
            
            # Get genome ID from mapping
            genome_id = protein_to_genome.get(neighbour, "Unknown")
            
            match_data = {
                'psi_protein': psi,
                'neighbour_protein': neighbour,
                'percent_identity': identity,
                'genome_id': genome_id,
                'neighbour_domains': format_domains(get_domains(neighbour, neighbour_domains)),
                'psi_domains': format_domains(get_domains(psi, psi_domains))
            }
            
            # Categorize the match
            for cat_name, cat in categories.items():
                if cat['min'] < identity <= cat['max']:
                    cat['matches'].append(match_data)
                    break
    
    return categories

def write_category_csv(base_path, category_name, matches):
    """Write matches for a specific category to CSV."""
    if not matches:
        print(f"⚠️ No matches found in category {category_name}")
        return
    
    output_csv = os.path.join(base_path, f"matches_{category_name}.csv")
    fieldnames = ['psi_protein', 'neighbour_protein', 'percent_identity', 
                 'genome_id', 'neighbour_domains', 'psi_domains']
    
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matches)
    
    print(f"📝 {len(matches)} sequences written to {output_csv}")


def main():
    parser = argparse.ArgumentParser(description="Protein analysis pipeline with optimized high-identity pair filtering.")
    parser.add_argument('-b', '--Base_filepath', required=True, help="Base file path for input files")
    parser.add_argument('-d', '--database_hmm', required=True, help="Path to HMM database", nargs="+")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of CPU threads")
    parser.add_argument('--identity_threshold', type=float, default=95.0,
                      help="Minimum percent identity for filtering (default: 95)")
    args = parser.parse_args()

    base_filepath = args.Base_filepath
    
    #file path to combined_cumulative_niehgbours
    combined_cumulative_niehgbours = os.path.join(base_filepath, "Neighborhoods_results", "txt_files", "combined_cumulative.txt")

    #file path to domain summary file from psi blast search
    Psi_domain_summary = os.path.join(base_filepath, "Domains_results", f"{base_filepath}_summary.txt")

    #file path to filter file
    Psi_protein_ids = os.path.join(base_filepath, "Psi_blast_results", f"{base_filepath}_Psi_blast_ids.txt")

    #file path to cumulative domain summary file from psi blast search
    Psi_cumulative_file = os.path.join(base_filepath, "Domains_results", f"{base_filepath}_cumulative_summary.txt")

    #file path to the output directory
    output_dir = os.path.join(base_filepath, "Domain_analysis")
    os.makedirs(output_dir, exist_ok=True)

    #file path to common_hmm dir
    common_hmm_dir = os.path.join(base_filepath, "Domain_analysis", "Common_HMMs")
    os.makedirs(common_hmm_dir, exist_ok=True)

    #file path to matches directory
    matches_dir = os.path.join(base_filepath, "Domain_analysis", "Common_HMMs", "Matches")
    os.makedirs(matches_dir, exist_ok=True)

    #file path to csv directory
    csv_directory = os.path.join(base_filepath, "Neighborhoods_results", "csv_files")





    # Step 1: Extract and filter Protein IDs
    protein_ids = read_protein_ids(combined_cumulative_niehgbours, Psi_protein_ids)
    print(f"📊 Final Protein IDs after filtering: {len(protein_ids)}")

    # Step 2: Fetch sequences
    fasta_file = os.path.join(output_dir, 'sequences.fasta')
    fetch_sequences(protein_ids, fasta_file)

    # Step 3: Build BLAST DB
    db_name = os.path.join(output_dir, 'blast_db')
    run_makeblastdb(fasta_file, db_name)

    # Step 4: Run BLASTP
    blast_output = os.path.join(output_dir, 'blast_results.tsv')
    run_blastp(fasta_file, db_name, blast_output, threads=args.threads)

    # Step 5: Filter high identity pairs
    high_id_output = os.path.join(output_dir, 'filtered_pairs.tsv')
    filter_high_identity_pairs(blast_output, high_id_output, args.identity_threshold)

    # Step 6: Cluster proteins
    clusters = os.path.join(output_dir, 'protein_clusters.csv')
    find_protein_clusters(high_id_output, clusters, args.identity_threshold)

    # Step 7: Filter FASTA by clusters
    filtered_fasta_file = os.path.join(output_dir, 'filtered_sequences.fasta')
    filter_fasta_by_clusters(clusters, fasta_file, filtered_fasta_file)

    # Step 8: Run hmmscan
    output_tblout_file = os.path.join(output_dir, 'hmmscan_results.tblout')
    run_hmmscan(filtered_fasta_file, output_tblout_file, args.database_hmm)
    domain_hits = parse_hmmscan_output(output_tblout_file)

    output_summary_file = os.path.join(output_dir, "domain_summary.txt")
    write_domain_hits(domain_hits, output_summary_file)

    output_cumulative_file = os.path.join(output_dir, "cumulative_summary.txt")
    summarize_results(output_summary_file, output_cumulative_file, filtered_fasta_file)

    # Step 9: Compare HMM accessions between neighbours and the psi
    compare_hmm_accessions_with_description(Psi_domain_summary, output_summary_file, output_dir)

    # Step 10: Clean up BLAST DB files
    delete_blastdb_files(db_name)

    # Step 11: Find common HMM accessions with counts
    out_path = find_common_hmm_accessions_with_counts(output_cumulative_file, Psi_cumulative_file, common_hmm_dir)

    # Step 12: common hmm blast
    blast_output, psi_fasta, neigh_fasta = common_hmm_blast_main(out_path, Psi_domain_summary, output_summary_file, common_hmm_dir, args.identity_threshold, args)

    # Get all data
    all_pairs = extract_all_pairs(blast_output)
    print(f"\nFound {len(all_pairs)} protein pairs total")
    
    fasta1_ids = extract_ids_from_fasta(neigh_fasta)
    print(f"Found {len(neigh_fasta)} IDs in neighbour proteins")
    
    fasta2_ids = extract_ids_from_fasta(psi_fasta)
    print(f"Found {len(psi_fasta)} IDs in PSI proteins")

    # Load domain information
    neighbour_domains = parse_domain_summary(output_summary_file)
    psi_domains = parse_domain_summary(Psi_domain_summary)
    print(f"Loaded domains for {len(neighbour_domains)} neighbour proteins")
    print(f"Loaded domains for {len(psi_domains)} PSI proteins")

    # Build genome ID mapping
    protein_to_genome = build_genome_id_map(csv_directory)
    print(f"Loaded genome IDs for {len(protein_to_genome)} proteins")

    # Categorize matches by identity ranges
    categories = categorize_matches(all_pairs, fasta1_ids, fasta2_ids, 
                                  neighbour_domains, psi_domains, protein_to_genome)
    
    # Write separate CSV files for each category
    for cat_name, cat_data in categories.items():
        #if cat_name == 'below_40':  # Skip if you don't want this category
            #continue
        print(f"\nProcessing category {cat_name} ({cat_data['min']}-{cat_data['max']}% identity)")
        write_category_csv(matches_dir, cat_name, cat_data['matches'])

    print("\n✅ All categories processed") 

    
    print("✨ Analysis complete!")
if __name__ == '__main__':
    main()
