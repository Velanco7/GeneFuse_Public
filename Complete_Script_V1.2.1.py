#!/usr/bin/env python3
"""
GeneFuse: Advanced Genomic Neighborhood and Domain Analysis Pipeline

Personal Project

Author: Velanco Fernandes
Email: Velancofernandes7@googlemail.com
"""
########################################################################################################################

import os
import argparse
import subprocess
import time
import re
import shutil
import pandas as pd
import csv
import logging
import requests
import gzip
import multiprocessing
import configparser
import urllib.request
import tarfile
import glob
from datetime import datetime
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from collections import Counter, defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from openpyxl import Workbook
from BCBio import GFF
from lxml import etree
from time import sleep
from threading import Event
from pathlib import Path

# Constants
DEFAULT_CONFIG = {
    'general': {
        'max_retries': '3',
        'timeout': '30',
        'log_level': 'INFO'
    },
    'ncbi': {
        'BATCH_SIZE': '2000',
        'batch_size': '2000',
        'delay': '0.34'
    },
    'analysis': {
        'e_value_threshold': '0.001',
        'domain_e_value_threshold': '0.001'
    }
}

# Add these new constants for NCBI operations
BATCH_SIZE = 2000  # Number of sequences to fetch in each batch
DELAY = 0.34       # Delay between NCBI requests (seconds)

class PipelineLogger:
    """Enhanced logging system with dual file/console output and colored formatting"""
    
    COLORS = {
        'DEBUG': '\033[36m',    # Cyan
        'INFO': '\033[32m',     # Green
        'WARNING': '\033[33m',   # Yellow
        'ERROR': '\033[31m',     # Red
        'CRITICAL': '\033[31m'   # Red
    }
    RESET = '\033[0m'

    def __init__(self, output_dir):
        self.log_file = os.path.join(output_dir, "genomic_neighborhood_execution.log")
        self.logger = self._setup_logging()

    def _setup_logging(self):
        logger = logging.getLogger('GeneFuse')
        logger.setLevel(logging.DEBUG)

        # Clear existing handlers
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        # File handler
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = self.ColoredFormatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        return logger

    # Standard logging methods
    def debug(self, message):
        self.logger.debug(message)
        print(f"{self.COLORS['DEBUG']}[{datetime.now().strftime('%H:%M:%S')}] DEBUG: {message}{self.RESET}")

    def info(self, message):
        self.logger.info(message)
        print(f"{self.COLORS['INFO']}[{datetime.now().strftime('%H:%M:%S')}] INFO: {message}{self.RESET}")

    def warning(self, message):
        self.logger.warning(message)
        print(f"{self.COLORS['WARNING']}[{datetime.now().strftime('%H:%M:%S')}] WARNING: {message}{self.RESET}")

    def error(self, message):
        self.logger.error(message)
        print(f"{self.COLORS['ERROR']}[{datetime.now().strftime('%H:%M:%S')}] ERROR: {message}{self.RESET}")

    def critical(self, message):
        self.logger.critical(message)
        print(f"{self.COLORS['CRITICAL']}[{datetime.now().strftime('%H:%M:%S')}] CRITICAL: {message}{self.RESET}")

    # Backwards-compatible log method
    def log(self, message, level='info'):
        log_method = getattr(self, level.lower(), self.info)
        log_method(message)

    class ColoredFormatter(logging.Formatter):
        def format(self, record):
            message = super().format(record)
            return f"{PipelineLogger.COLORS.get(record.levelname, '')}{message}{PipelineLogger.RESET}"

class ConfigManager:
    """Centralized configuration management with file support"""
    
    def __init__(self, config_path=None):
        self.config = configparser.ConfigParser()
        self.config.read_dict(DEFAULT_CONFIG)
        
        if config_path:
            self.load_config(config_path)

    def load_config(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Config file not found: {path}")
        self.config.read(path)

    def get(self, section, key, fallback=None):
        try:
            return self.config.get(section, key)
        except (configparser.NoSectionError, configparser.NoOptionError):
            return fallback or DEFAULT_CONFIG.get(section, {}).get(key)

class NetworkManager:
    """Manages network connectivity state and wait operations"""
    def __init__(self, logger):
        self.logger = logger
        self._connection_available = True  # Changed from multiprocessing.Value to regular boolean
        self.test_url = "http://www.google.com"
        self.check_interval = 30

    @property
    def connection_available(self):
        return self._connection_available

    @connection_available.setter
    def connection_available(self, value):
        self._connection_available = value

    def check_connection(self):
        """Check if network connection is available"""
        try:
            urllib.request.urlopen(self.test_url, timeout=5)
            if not self.connection_available:
                self.logger.log("Network connection restored", "info")
            self.connection_available = True
            return True
        except Exception:
            if self.connection_available:
                self.logger.log("Network connection lost", "warning")
            self.connection_available = False
            return False
            
    def wait_for_connection(self):
        """Block until network connection is available"""
        self.logger.log("Waiting for network connection...", "info")
        while not self.connection_available:
            time.sleep(1)
        self.logger.log("Network connection available, resuming", "info")

def network_required(func):
    """Decorator to pause functions when network is unavailable"""
    def wrapper(*args, **kwargs):
        # The first argument should be the class instance (self)
        network_mgr = args[0].network_mgr if hasattr(args[0], 'network_mgr') else None
        
        if network_mgr is None:
            # If no network manager, just run the function
            return func(*args, **kwargs)
            
        while True:
            if network_mgr.check_connection():
                try:
                    return func(*args, **kwargs)
                except urllib.error.URLError as e:
                    network_mgr.logger.log(f"Network error during operation: {str(e)}", "warning")
                    network_mgr.wait_for_connection()
                    continue
            else:
                network_mgr.wait_for_connection()
    return wrapper

def parse_arguments():
    parser = argparse.ArgumentParser(description="GeneFuse: gene fusion detection tool")
    parser.add_argument("input_file", help="Input FASTA file containing query sequence")
    
    parser.add_argument("-u", "--user_email", required=True, help="User Email address (required)")
    parser.add_argument("-e1", "--E_value_psi", help="User defined E-Value cutoff for Psiblast, default = 0.00001", default="0.00001")
    parser.add_argument("-iterations", "--iterations", help="User defined number of iterations for Psiblast, default = 3 iterations", default="3")
    parser.add_argument("-db1", "--database_psi", required=True, help="Database to be used for psiblast search")
    parser.add_argument("-db2", "--database_hmm", required=True, help="HMM databases (e.g., Pfam-A.hmm)", nargs="+")
    parser.add_argument("-g", "--num_genes", help="Number of upstream and downstream genes to extract (default: 5)", default=5, type=int)
    parser.add_argument("-num_threads", help="Number of threads (default: 4)", default=4, type=int)
    parser.add_argument("-api_key", help="Entrez api key (optional)")
    parser.add_argument('--identity_threshold', type=float, default=95.0,
                      help="Minimum percent identity for filtering (default: 95)")
    # Config file
    parser.add_argument("-c", "--config", default="config.ini", help="Path to configuration file")
    
    return parser.parse_args()

def output_directory_initialisation(input_name, config):
    base_name = os.path.splitext(os.path.basename(input_name))[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = f"GeneFuse_{base_name}_{timestamp}"

    # Create main directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Setup subdirectories
    subdirs = {
        'psi_results': os.path.join(output_dir, "Psi_blast_results"),
        'neighborhoods': os.path.join(output_dir, "Neighborhoods_results"),
        'domains': os.path.join(output_dir, "Domains_results"),
        'temp': os.path.join(output_dir, "Temp_files"),
        'domain_analysis': os.path.join(output_dir, "Domain_analysis")
    }
    
    for path in subdirs.values():
        os.makedirs(path, exist_ok=True)

    # Copy input file
    input_copy_path = os.path.join(output_dir, "input_file.txt")
    shutil.copy2(input_name, input_copy_path)

    return output_dir, subdirs

def run_psiblast(sequence_file, E_value, iterations, database, output_file, num_threads, PSSM_file, logger):
    """
    Run PSI-BLAST for the specified number of iterations, saving only the results of the final iteration.
    """
    try:
        logger.log(f"Running PSI-BLAST with E-Value {E_value} and for {iterations} iterations on {database} database", "info")

        max_target_seqs = 100000  # Maximum number of target sequences to report

        command = [
            'psiblast',
            '-query', sequence_file,
            '-db', database,
            '-evalue', str(E_value),
            '-num_iterations', str(iterations),
            '-outfmt', '5',  # XML output format
            '-out', output_file,
            '-max_target_seqs', str(max_target_seqs),
            '-num_threads', str(num_threads),
            '-out_ascii_pssm', PSSM_file,
            '-save_pssm_after_last_round'
        ]

        if os.name == 'posix' and os.uname().sysname == 'Linux':
            command = ['ulimit', '-v', '4000000'] + command  # 4GB memory limit
            
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        logger.log("PSI-BLAST completed successfully", "info")
        return output_file

    except subprocess.CalledProcessError as e:
        error_msg = f"PSI-BLAST failed with error: {e.stderr}"
        if "SIGBUS" in str(e):
            error_msg += "\nThis is typically a memory error. Try reducing threads (-num_threads) or using a smaller database."
        logger.log(error_msg, "error")
        raise
    except Exception as e:
        logger.log(f"Unexpected error in PSI-BLAST: {str(e)}", "error")
        raise

def process_blast_xml(xml_file, logger):
    """Process BLAST XML output and extract hits (one sequence per unique hit_def)"""
    try:
        logger.log(f"Processing BLAST XML file: {xml_file}", "info")
        
        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            
            blast_results = []
            protein_ids = set()
            fasta_sequences = []
            seen_hit_defs = set()  # Track unique hit definitions
            
            for record in blast_records:
                if not record.alignments:
                    logger.log(f"No alignments found for query: {record.query}", "debug")
                    continue
                
                for alignment in record.alignments:
                    hit_def = f"{alignment.accession} {alignment.hit_def}"
                    
                    if hit_def in seen_hit_defs:
                        continue  # Skip duplicates
                    
                    for hsp in alignment.hsps:
                        try:
                            protein_id = alignment.accession.split('|')[0] if alignment.accession else "Unknown"
                            protein_ids.add(protein_id)
                            
                            # Prepare CSV row
                            blast_results.append([
                                record.query,
                                hit_def,
                                str(round((hsp.identities / float(hsp.align_length)) * 100, 3)),
                                str(hsp.align_length),
                                str(hsp.align_length - hsp.identities),
                                str(hsp.gaps if hasattr(hsp, 'gaps') else 0),
                                str(hsp.query_start),
                                str(hsp.query_end),
                                str(hsp.sbjct_start),
                                str(hsp.sbjct_end),
                                str(hsp.expect),
                                str(hsp.bits),
                                hsp.sbjct
                            ])
                            
                            # Prepare FASTA entry
                            fasta_sequences.append(f">{hit_def}\n{hsp.sbjct}")
                            
                            seen_hit_defs.add(hit_def)
                            break
                            
                        except Exception as e:
                            logger.log(f"Error processing alignment: {str(e)}", "warning")
                            continue
            
            return blast_results, protein_ids, fasta_sequences
            
    except Exception as e:
        logger.log(f"Error processing BLAST XML: {str(e)}", "error")
        raise

def clean_fasta(input_fasta, cleaned_fasta, logger):
    """
    Cleans the input FASTA file by removing invalid characters (e.g., '-') from sequences.
    """
    try:
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

        valid_sequences = []
        logger.log(f"Cleaning FASTA file: {input_fasta}", "info")

        with open(input_fasta, "r") as infile, open(cleaned_fasta, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                cleaned_sequence = str(record.seq).replace("-", "")
                record.seq = Seq(cleaned_sequence)

                if cleaned_sequence:
                    valid_sequences.append(record)
                else:
                    logger.log(f"Invalid sequence removed: {record.id}", "warning")

            SeqIO.write(valid_sequences, outfile, "fasta")

        if not valid_sequences:
            raise ValueError("All sequences were removed during cleaning. Check your input file.")
        
        logger.log(f"Cleaned FASTA file created: {cleaned_fasta}", "info")
        return cleaned_fasta

    except Exception as e:
        logger.log(f"Error cleaning FASTA file: {str(e)}", "error")
        raise

def run_hmmscan(input_fasta, output_file, hmm_databases, logger):
    """
    Runs hmmscan using the input FASTA file and the given HMM databases.
    Combines results from multiple databases into a single output file.
    """
    try:
        logger.log("Running HMMER scans...", "info")
        logger.log(f"Input FASTA: {input_fasta}", "debug")
    
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

        # Temporary files to hold individual database results
        temp_files = []
        for hmm_database in hmm_databases:
            if not os.path.exists(hmm_database):
                logger.log(f"Error: HMM database not found: {hmm_database}", "error")
                continue
            
            logger.log(f"Scanning with database: {hmm_database}", "info")

            # Build the hmmscan command
            temp_file = f"{hmm_database.split('/')[-1]}_result.tblout"
            temp_files.append(temp_file)

            command = [
                "hmmscan",
                "--domtblout", temp_file,
                "-E", "0.001",
                "--domE", "0.001",
                hmm_database,
                input_fasta,
            ]

            result = subprocess.run(
                command, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True
            )
            
            if result.returncode != 0:
                logger.log(f"Error in hmmscan with database {hmm_database}: {result.stderr}", "error")
                continue

        # Combine the temporary files into the final output file
        with open(output_file, "w") as outfile:
            for temp_file in temp_files:
                try:
                    with open(temp_file, "r") as temp:
                        outfile.write(f"# Results for {temp_file}\n")
                        outfile.write(temp.read())
                        outfile.write("\n")
                    os.remove(temp_file)
                except Exception as e:
                    logger.log(f"Error processing temp file {temp_file}: {str(e)}", "error")

        logger.log(f"HMMER scan results saved to: {output_file}", "info")
    
    except subprocess.CalledProcessError as e:
        logger.log(f"Error running hmmscan: {str(e)}", "error")
        raise
    except FileNotFoundError:
        logger.log("Error: hmmscan not found. Ensure HMMER is installed and in PATH.", "error")
        raise

def parse_hmmscan_output(output_file, logger):
    """
    Parses the hmmscan output file and extracts domain hits.
    """
    try:
        domain_hits = []
        logger.log(f"Parsing hmmscan output: {output_file}", "info")

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
                    "start": parts[17],
                    "end": parts[18],
                    "description": " ".join(parts[22:]),
                }
                domain_hits.append(domain_hit)

        logger.log(f"Found {len(domain_hits)} domain hits", "info")
        return domain_hits

    except Exception as e:
        logger.log(f"Error parsing hmmscan output: {str(e)}", "error")
        raise

def write_domain_hits(domain_hits, output_summary_file, logger):
    """
    Writes the domain hits to a human-readable summary file.
    """
    try:
        logger.log(f"Writing domain hits to: {output_summary_file}", "info")

        with open(output_summary_file, "w") as f:
            f.write("Target Name\tAccession\tQuery Name\tE-value\tScore\tStart\tEnd\tDescription\n")

            for hit in domain_hits:
                try:
                    target_name = hit.get("target_name", "").strip()
                    accession = hit.get("accession", "").strip()
                    query_name = hit.get("query_name", "").strip()
                    e_value = hit.get("E-value", "").strip()
                    score = hit.get("score", "").strip()
                    start = hit.get("start", "").strip()
                    end = hit.get("end", "").strip()
                    description = hit.get("description", "").strip()

                    f.write(
                        f"{target_name}\t{accession}\t{query_name}\t{e_value}\t{score}\t"
                        f"{start}\t{end}\t{description}\n"
                    )
                except Exception as e:
                    logger.log(f"Error processing hit: {hit}. Error: {e}", "warning")

        logger.log(f"Domain summary saved to: {output_summary_file}", "info")

    except Exception as e:
        logger.log(f"Error writing domain hits: {str(e)}", "error")
        raise 

def write_gb_file(domain_hits, output_gb_file, query_sequences, logger):
    """
    Writes the domain hits to a .gb file with detailed information.
    """
    try:
        logger.log(f"Writing GenBank file: {output_gb_file}", "info")
    
        with open(output_gb_file, "w") as f:
            for seq_record in query_sequences:
                f.write(f"LOCUS\t{seq_record.id}\t{len(seq_record.seq)} bp\tPROTEIN\tUNDEFINED\n")
                f.write("FEATURES\tLocation/Qualifiers\n")
                f.write(f"     SOURCE\t1..{len(seq_record.seq)}\n")
                f.write(f"         /organism=\"Unknown\"\n")
                f.write(f"         /sequence=\"{seq_record.seq}\"\n")
                
                for i, hit in enumerate(domain_hits, 1):
                    if hit["query_name"] == seq_record.id:
                        f.write(f"     DOMAIN\t{hit['start']}..{hit['end']}\n")
                        f.write(f"         /ID=\"{i}\"\n")
                        f.write(f"         /Target=\"{hit['target_name']}\"\n")
                        f.write(f"         /Accession=\"{hit['accession']}\"\n")
                        f.write(f"         /E-value=\"{hit['E-value']}\"\n")
                        f.write(f"         /Score=\"{hit['score']}\"\n")
                        f.write(f"         /Description=\"{hit['description']}\"\n")
                f.write("//\n")

        logger.log(f"GenBank file created: {output_gb_file}", "info")
        
    except Exception as e:
        logger.log(f"Error writing GenBank file: {str(e)}", "error")
        raise

def summarize_results(summary_file, output_cumulative_file, fasta_file, logger):
    """
    Summarizes the results in the summary.txt file and writes a cumulative summary.
    """
    try:
        logger.log(f"Summarizing results from: {summary_file}", "info")
    
        target_counts = Counter()
        total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
        logger.log(f"Total sequences in input FASTA file: {total_sequences}", "info")

        with open(summary_file, "r") as f:
            for line_number, line in enumerate(f, 1):
                if line.startswith("Target Name"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 8:
                    logger.log(f"Skipping malformed line {line_number}", "warning")
                    continue
                
                target_name = parts[0]
                accession = parts[1]
                description = parts[7]
                target_counts[(target_name, accession, description)] += 1

        with open(output_cumulative_file, "w") as f:
            f.write("Target Name\tAccession\tDescription\tOccurrences\n")
            for (target_name, accession, description), count in target_counts.items():
                f.write(f"{target_name}\t{accession}\t{description}\t{count}\n")
            
            f.write(f"\nTotal Sequences Processed: {total_sequences}\n")
        
        logger.log(f"Cumulative summary created: {output_cumulative_file}", "info")
        
    except Exception as e:
        logger.log(f"Error summarizing results: {str(e)}", "error")
        raise

def read_accessions(file_path, logger):
    try:
        logger.log(f"Reading accessions from: {file_path}", "info")
        with open(file_path, "r") as file:
            accessions = [line.strip() for line in file]
        logger.log(f"Found {len(accessions)} accessions", "info")
        return accessions
    except Exception as e:
        logger.log(f"Error reading accessions: {str(e)}", "error")
        raise

def download_with_retries(url, file_name, max_retries=3, timeout=30, logger=None):
    """Download a file with retries and timeout handling."""
    session = requests.Session()
    for attempt in range(max_retries):
        try:
            if logger:
                logger.log(f"Attempt {attempt + 1} to download {url}", "info")
                
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            
            with open(file_name, "wb") as f:
                f.write(response.content)
                
            if logger:
                logger.log(f"Downloaded {file_name}", "info")
            return file_name
            
        except requests.exceptions.RequestException as e:
            if logger:
                logger.log(f"Attempt {attempt + 1} failed: {e}", "warning")
            if attempt < max_retries - 1:
                time.sleep(5)
            else:
                if logger:
                    logger.log(f"Failed to download {url} after {max_retries} attempts", "error")
                raise
    return None

@network_required
def download_assembly_summary_files(logger=None):
    try:
        genbank_url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
        refseq_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

        if logger:
            logger.log("Downloading assembly summary files", "info")

        genbank_file = download_with_retries(genbank_url, "assembly_summary_genbank.txt", logger=logger)
        refseq_file = download_with_retries(refseq_url, "assembly_summary_refseq.txt", logger=logger)
        
        if logger:
            logger.log("Successfully downloaded assembly summary files", "info")
            
        return genbank_file, refseq_file
        
    except Exception as e:
        if logger:
            logger.log(f"Error downloading assembly summaries: {str(e)}", "error")
        raise

def load_assembly_summary(file_path, logger=None):
    try:
        if logger:
            logger.log(f"Loading assembly summary: {file_path}", "info")
            
        column_names = ["#assembly_accession", "bioproject", "biosample", "ftp_path"]    
        df = pd.read_csv(file_path, sep="\t", skiprows=1, usecols=[0, 1, 2, 19], names=column_names, dtype=str)
        
        if logger:
            logger.log(f"Loaded {len(df)} assembly records", "info")
            
        return df
        
    except Exception as e:
        if logger:
            logger.log(f"Error loading assembly summary: {str(e)}", "error")
        raise

@network_required
def get_assembly_from_protein(protein_id, logger=None):
    """
    Retrieve assembly accessions associated with a given protein ID using NCBI Entrez API.
    """
    try:
        if logger:
            logger.log(f"Fetching nucleotide links for protein ID: {protein_id}", "info")
        
        handle = Entrez.elink(dbfrom="protein", id=protein_id, db="nuccore", retmode="xml")
        raw_xml = handle.read()
        handle.close()

        root = etree.fromstring(raw_xml)
        nucl_ids = []
        for link_set_db in root.xpath(".//LinkSetDb[DbTo='nuccore']"):
            for link in link_set_db.xpath(".//Link"):
                nucl_id = link.findtext("Id")
                if nucl_id:
                    nucl_ids.append(nucl_id)

        if not nucl_ids:
            if logger:
                logger.log(f"No nucleotide IDs found for protein ID {protein_id}", "warning")
            return set()

        if logger:
            logger.log(f"Found nucleotide IDs for protein ID {protein_id}: {nucl_ids}", "debug")

        assembly_uids = set()
        for nucl_id in nucl_ids:
            try:
                if logger:
                    logger.log(f"Fetching assembly links for nucleotide ID: {nucl_id}", "debug")

                handle = Entrez.elink(dbfrom="nuccore", id=nucl_id, db="assembly", retmode="xml")
                raw_xml = handle.read()
                handle.close()

                root = etree.fromstring(raw_xml)
                for link_set_db in root.xpath(".//LinkSetDb[DbTo='assembly']"):
                    for link in link_set_db.xpath(".//Link"):
                        assembly_uid = link.findtext("Id")
                        if assembly_uid:
                            assembly_uids.add(assembly_uid)

            except Exception as e:
                if logger:
                    logger.log(f"Error processing nucleotide ID {nucl_id}: {e}", "warning")
                continue

        if not assembly_uids:
            if logger:
                logger.log(f"No assembly UIDs found for nucleotide IDs: {nucl_ids}", "warning")
            return set()

        if logger:
            logger.log(f"Found assembly UIDs for protein ID {protein_id}: {assembly_uids}", "debug")

        assembly_accessions = set()
        try:
            if logger:
                logger.log(f"Fetching assembly summaries for assembly UIDs: {assembly_uids}", "debug")

            handle = Entrez.esummary(db="assembly", id=",".join(assembly_uids), retmode="xml")
            raw_xml = handle.read()
            handle.close()

            root = etree.fromstring(raw_xml)
            for doc in root.xpath(".//DocumentSummary"):
                assembly_accession = doc.findtext("AssemblyAccession")
                if assembly_accession:
                    assembly_accessions.add(assembly_accession)

        except Exception as e:
            if logger:
                logger.log(f"Error fetching assembly summaries: {e}", "error")
            return set()

        if not assembly_accessions:
            if logger:
                logger.log(f"No assembly accessions found for assembly UIDs: {assembly_uids}", "warning")
            return set()

        if logger:
            logger.log(f"Found assembly accessions for protein ID {protein_id}: {assembly_accessions}", "info")
            
        return assembly_accessions

    except Exception as e:
        if logger:
            logger.log(f"Unexpected error processing protein ID {protein_id}: {e}", "error")
        return set()

def find_ftp_link(assembly_id, genbank_df, refseq_df, logger=None):
    try:
        if logger:
            logger.log(f"Finding FTP link for assembly: {assembly_id}", "debug")
            
        genbank_match = genbank_df[genbank_df["#assembly_accession"] == assembly_id]
        if not genbank_match.empty:
            ftp_path = genbank_match.iloc[0]["ftp_path"]
            if logger:
                logger.log(f"Found GenBank FTP link for {assembly_id}: {ftp_path}", "debug")
            return ftp_path
            
        refseq_match = refseq_df[refseq_df["#assembly_accession"] == assembly_id]
        if not refseq_match.empty:
            ftp_path = refseq_match.iloc[0]["ftp_path"]
            if logger:
                logger.log(f"Found RefSeq FTP link for {assembly_id}: {ftp_path}", "debug")
            return ftp_path
            
        if logger:
            logger.log(f"No FTP link found for assembly: {assembly_id}", "warning")
        return None
        
    except Exception as e:
        if logger:
            logger.log(f"Error finding FTP link: {str(e)}", "error")
        return None

@network_required
def download_gff3_file(ftp_link, output_dir, assembly_id, logger=None):
    try:
        if logger:
            logger.log(f"Downloading GFF3 for {assembly_id} from {ftp_link}", "info")
            
        file_name = os.path.basename(ftp_link)
        gff3_url = f"{ftp_link}/{file_name}_genomic.gff.gz"
        output_file_path = os.path.join(output_dir, f"{assembly_id}.gff3.gz")

        response = requests.get(gff3_url, timeout=30)
        if response.status_code == 200:
            with open(output_file_path, "wb") as output_file:
                output_file.write(response.content)
                
            if logger:
                logger.log(f"Downloaded .gff3 file for {assembly_id}", "info")
            return output_file_path
        else:
            if logger:
                logger.log(f"Failed to download .gff3 file for {assembly_id}", "error")
            return None
    except Exception as e:
        if logger:
            logger.log(f"Error downloading .gff3 file: {str(e)}", "error")
        return None

def gunzip_file(gz_file_path, output_dir, logger=None):
    try:
        if logger:
            logger.log(f"Decompressing {gz_file_path}", "info")
            
        file_name = os.path.basename(gz_file_path).replace(".gz", "")
        output_file_path = os.path.join(output_dir, file_name)

        if not os.path.exists(gz_file_path):
            raise FileNotFoundError(f"Input file {gz_file_path} does not exist")

        with gzip.open(gz_file_path, "rb") as f_in, open(output_file_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        if os.path.exists(output_file_path):
            if logger:
                logger.log(f"Successfully decompressed to {output_file_path}", "info")
            os.remove(gz_file_path)
            if logger:
                logger.log(f"Deleted compressed file: {gz_file_path}", "debug")
            return output_file_path
        else:
            raise IOError(f"Decompression failed: {output_file_path} does not exist")
            
    except Exception as e:
        if logger:
            logger.log(f"Error decompressing file: {str(e)}", "error")
        raise

def find_neighboring_proteins(gff_file, target_protein_id, num_genes, logger=None):
    """
    Find proteins upstream and downstream of a target protein in a GFF file.
    """
    try:
        if logger:
            logger.log(f"Finding neighbors for {target_protein_id} in {gff_file}", "info")

        if not target_protein_id.endswith(".1"):
            target_protein_id += ".1"

        target_gene = None

        with open(gff_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")
                if len(columns) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = columns
                
                if feature_type != "CDS":
                    continue

                attr_dict = {}
                for attr in attributes.split(";"):
                    if "=" in attr:
                        key, value = attr.split("=", 1)
                        attr_dict[key] = value

                if "protein_id" in attr_dict and attr_dict["protein_id"] == target_protein_id:
                    target_gene = (seqid, int(start), int(end), strand, attr_dict)
                    break

        if not target_gene:
            if logger:
                logger.log(f"Target protein ID '{target_protein_id}' not found", "warning")
            return []

        neighboring_genes = []
        target_seqid, target_start, target_end, target_strand, target_attrs = target_gene

        with open(gff_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")
                if len(columns) < 9:
                    continue

                seqid, source, feature_type, start, end, score, strand, phase, attributes = columns

                if seqid != target_seqid or feature_type != "CDS":
                    continue

                attr_dict = {}
                for attr in attributes.split(";"):
                    if "=" in attr:
                        key, value = attr.split("=", 1)
                        attr_dict[key] = value

                start = int(start)
                end = int(end)

                if start == target_start and end == target_end:
                    position = "target"
                elif end < target_start:
                    position = "upstream"
                else:
                    position = "downstream"

                protein_id = attr_dict.get("protein_id", "Unknown")
                product = attr_dict.get("product", "Unknown")

                neighboring_genes.append({
                    "Sequence ID": seqid,
                    "Position": position,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "Protein ID": protein_id,
                    "Product": product
                })

        neighboring_genes.sort(key=lambda x: x["Start"])

        target_index = next((i for i, gene in enumerate(neighboring_genes) if gene["Position"] == "target"), None)
        if target_index is None:
            if logger:
                logger.log("Target gene not found in sorted gene list", "warning")
            return []

        result = neighboring_genes[max(0, target_index - num_genes):target_index + num_genes + 1]

        return result

    except Exception as e:
        if logger:
            logger.log(f"Error finding neighboring proteins: {str(e)}", "error")
        raise

def process_csv_files(csv_dir, output_dir, logger=None):
    """
    Process all .csv files in the given directory to extract protein IDs and descriptions,
    create cumulative frequency tables for each .csv file, and generate a combined cumulative file.
    """
    try:
        if logger:
            logger.log(f"Processing CSV files in {csv_dir}", "info")
        
        txt_files_dir = os.path.join(output_dir, "txt_files")
        os.makedirs(txt_files_dir, exist_ok=True)

        combined_cumulative_freq = {}
        description_groups = {}
        hypothetical_protein_groups = {}

        for csv_file in os.listdir(csv_dir):
            if not csv_file.endswith(".csv"):
                continue
                
            csv_file_path = os.path.join(csv_dir, csv_file)
            cumulative_freq = {}

            if logger:
                logger.log(f"Processing file: {csv_file}", "debug")

            with open(csv_file_path, "r") as file:
                reader = csv.reader(file)
                header = next(reader)

                for row_number, row in enumerate(reader, start=2):
                    try:
                        if len(row) < 7:
                            if logger:
                                logger.log(f"Skipping row {row_number}: insufficient columns", "warning")
                            continue

                        protein_id = row[5]
                        description = row[6]

                        if protein_id != "Unknown":
                            cumulative_freq[protein_id] = cumulative_freq.get(protein_id, 0) + 1
                            combined_cumulative_freq[protein_id] = combined_cumulative_freq.get(protein_id, 0) + 1

                        if description != "Unknown":
                            if description.lower() == "hypothetical protein":
                                if protein_id not in hypothetical_protein_groups:
                                    hypothetical_protein_groups[protein_id] = []
                                hypothetical_protein_groups[protein_id].append((protein_id, description))
                            else:
                                if description not in description_groups:
                                    description_groups[description] = {}
                                description_groups[description][protein_id] = description_groups[description].get(protein_id, 0) + 1

                    except Exception as e:
                        if logger:
                            logger.log(f"Error processing row {row_number}: {e}", "error")
                        continue

            cumulative_file_path = os.path.join(txt_files_dir, f"{os.path.splitext(csv_file)[0]}_cumulative.txt")
            with open(cumulative_file_path, "w") as file:
                file.write("Protein ID\tFrequency\tDescription\n")
                for protein_id, frequency in sorted(cumulative_freq.items(), key=lambda x: x[1], reverse=True):
                    description = next((row[6] for row in csv.reader(open(csv_file_path)) if len(row) > 6 and row[5] == protein_id), "Unknown")
                    file.write(f"{protein_id}\t{frequency}\t{description}\n")
                    
            if logger:
                logger.log(f"Saved cumulative file: {cumulative_file_path}", "debug")

        combined_cumulative_file_path = os.path.join(txt_files_dir, "combined_cumulative.txt")
        with open(combined_cumulative_file_path, "w") as file:
            file.write("Protein ID\tFrequency\tDescription\n")
            for protein_id, frequency in sorted(combined_cumulative_freq.items(), key=lambda x: x[1], reverse=True):
                description = next((row[6] for csv_file in os.listdir(csv_dir) if csv_file.endswith(".csv") 
                                   for row in csv.reader(open(os.path.join(csv_dir, csv_file))) if len(row) > 6 and row[5] == protein_id), "Unknown")
                file.write(f"{protein_id}\t{frequency}\t{description}\n")
                
        if logger:
            logger.log(f"Saved combined cumulative file: {combined_cumulative_file_path}", "info")

        grouped_file_path = os.path.join(txt_files_dir, "grouped_by_description.txt")
        with open(grouped_file_path, "w") as file:
            file.write("Most Frequent Protein ID\tFrequency\tDescription\n")

            combined_groups = {}

            for description, protein_ids in description_groups.items():
                total_freq = sum(protein_ids.values())
                most_frequent_protein_id = max(protein_ids, key=protein_ids.get)
                combined_groups[description] = {
                    "most_frequent_protein_id": most_frequent_protein_id,
                    "total_freq": total_freq,
                    "description": description
                }

            for protein_id, entries in hypothetical_protein_groups.items():
                total_freq = len(entries)
                combined_groups[f"hypothetical protein ({protein_id})"] = {
                    "most_frequent_protein_id": protein_id,
                    "total_freq": total_freq,
                    "description": "hypothetical protein"
                }

            for group_key, group_data in sorted(combined_groups.items(), key=lambda x: x[1]["total_freq"], reverse=True):
                file.write(f"{group_data['most_frequent_protein_id']}\t{group_data['total_freq']}\t{group_data['description']}\n")
                
        if logger:
            logger.log(f"Saved grouped results: {grouped_file_path}", "info")
            
    except Exception as e:
        if logger:
            logger.log(f"Error processing CSV files: {str(e)}", "error")
        raise

def domain_analysis(args, output_dir, psi_output_file, domains_sub_dir, logger):
    try:
        logger.log("Starting domain analysis", "info")
        
        cleaned_fasta_file = os.path.join(domains_sub_dir, f"{os.path.basename(output_dir)}_cleaned.fasta")
        output_tblout_file = os.path.join(domains_sub_dir, f"{os.path.basename(output_dir)}.tblout")
        output_summary_file = os.path.join(domains_sub_dir, f"{os.path.basename(output_dir)}_summary.txt")
        output_gb_file = os.path.join(domains_sub_dir, f"{os.path.basename(output_dir)}.gb")
        output_cumulative_file = os.path.join(domains_sub_dir, f"{os.path.basename(output_dir)}_cumulative_summary.txt")

        # Ensure the directory exists
        os.makedirs(domains_sub_dir, exist_ok=True)

        cleaned_fasta = clean_fasta(psi_output_file, cleaned_fasta_file, logger)
        query_sequences = list(SeqIO.parse(cleaned_fasta, "fasta"))

        run_hmmscan(cleaned_fasta_file, output_tblout_file, args.database_hmm, logger)
        
        domain_hits = parse_hmmscan_output(output_tblout_file, logger)
        write_domain_hits(domain_hits, output_summary_file, logger)
        write_gb_file(domain_hits, output_gb_file, query_sequences, logger)
        summarize_results(output_summary_file, output_cumulative_file, cleaned_fasta, logger)
        
        logger.log("Domain analysis completed successfully", "info")
        
    except Exception as e:
        logger.log(f"Domain analysis failed: {str(e)}", "error")
        raise

def neighborhood_analysis(args, output_dir, psi_id_output_file, hoods_sub_dir, num_genes, logger, genbank_df=None, refseq_df=None):
    try:
        logger.log("Starting neighborhood analysis", "info")

        domain_analysis_dir = os.path.join(output_dir, "Domain_analysis")
        os.makedirs(domain_analysis_dir, exist_ok=True)

        gff_dir = os.path.join(hoods_sub_dir, "gff_files")
        csv_dir = os.path.join(hoods_sub_dir, "csv_files")
        txt_dir = os.path.join(hoods_sub_dir, "txt_files")
        
        os.makedirs(gff_dir, exist_ok=True)
        os.makedirs(csv_dir, exist_ok=True)
        os.makedirs(txt_dir, exist_ok=True)
        
        if genbank_df is None or refseq_df is None:
            logger.log("Downloading assembly summary files", "info")
            genbank_file, refseq_file = download_assembly_summary_files(logger)
            genbank_df = load_assembly_summary(genbank_file, logger)
            refseq_df = load_assembly_summary(refseq_file, logger)
        else:
            logger.log("Using pre-loaded assembly data", "info")

        processed_assemblies_per_protein = {}
        accessions = read_accessions(psi_id_output_file, logger)

        for accession in accessions:
            if accession not in processed_assemblies_per_protein:
                processed_assemblies_per_protein[accession] = set()

            assembly_ids = get_assembly_from_protein(accession, logger)
            for assembly_id in assembly_ids:
                if assembly_id in processed_assemblies_per_protein[accession]:
                    continue

                processed_assemblies_per_protein[accession].add(assembly_id)
                ftp_link = find_ftp_link(assembly_id, genbank_df, refseq_df, logger)
                
                if not ftp_link:
                      continue
                
                logger.log(f"Processing assembly {assembly_id} for protein {accession}", "info")
                        
                gff_gz_file = download_gff3_file(ftp_link, hoods_sub_dir, assembly_id, logger)
                if not gff_gz_file:
                    continue

                try:
                    gff_file_path = gunzip_file(gff_gz_file, gff_dir, logger)
                    results = find_neighboring_proteins(gff_file_path, accession, num_genes, logger)
                    
                    if results:
                        output_csv = os.path.join(csv_dir, f"{accession}_neighbors.csv")
                        file_exists = os.path.exists(output_csv)

                        with open(output_csv, "a", newline="") as csvfile:
                            writer = csv.writer(csvfile)
                            
                            if not file_exists:
                                writer.writerow(["Sequence ID", "Position", "Start", "End", "Strand", "Protein ID", "Product"])
                                writer.writerow([])

                            if file_exists:
                                writer.writerow([])
                                writer.writerow(["----------------------------------------------------------------"])
                                writer.writerow([])

                            upstream_written = False
                            target_written = False
                            downstream_written = False

                            for result in results:
                                if result["Position"] == "target" and not target_written:
                                    if upstream_written:
                                        writer.writerow([])
                                    target_written = True

                                if result["Position"] == "downstream" and not downstream_written:
                                    if target_written:
                                        writer.writerow([])
                                    downstream_written = True

                                writer.writerow([
                                    result["Sequence ID"],
                                    result["Position"],
                                    result["Start"],
                                    result["End"],
                                    result["Strand"],
                                    result["Protein ID"],
                                    result["Product"]
                                ])

                                if result["Position"] == "upstream":
                                    upstream_written = True

                        logger.log(f"Saved neighborhood results to {output_csv}", "info")

                    os.remove(gff_file_path)
                    logger.log(f"Removed temporary GFF file: {gff_file_path}", "debug")
                    
                except Exception as e:
                    logger.log(f"Error processing GFF for {assembly_id}: {str(e)}", "error")
                    continue

        process_csv_files(csv_dir, hoods_sub_dir, logger)

        try:
            shutil.rmtree(gff_dir)
            logger.log(f"Removed temporary directory: {gff_dir}", "debug")
        except Exception as e:
            logger.log(f"Error removing directory {gff_dir}: {str(e)}", "warning")
            
        logger.log("Neighborhood analysis completed successfully", "info")
        
    except Exception as e:
        logger.log(f"Neighborhood analysis failed: {str(e)}", "error")
        raise

def read_protein_ids(filepath, filter_file=None, logger=None):
    """Extract Protein IDs from the input file, optionally filtering against another file."""
    try:
        protein_ids = set()
        filter_ids = set()
        
        if filter_file:
            logger.log("Reading filter file...", "info")
            with open(filter_file, 'r') as f:
                filter_ids = {line.strip() + ".1" for line in f if line.strip()}
            logger.log(f"Total protein IDs in filter file (with '.1' appended): {len(filter_ids)}", "info")
        
        with open(filepath, 'r') as file:
            next(file)
            for line in file:
                if line.strip() == "":
                    continue
                parts = line.strip().split('\t')
                if parts and parts[0] not in filter_ids:
                    protein_ids.add(parts[0])
        return list(protein_ids)

    except Exception as e:
        logger.log(f"Neighbourhood proteins filtering failed: {str(e)}", "error")
        raise

@network_required
def fetch_sequences(protein_ids, output_fasta, logger):
    """Fetch sequences in batches via Entrez."""
    try:
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
                        logger.log(f"✅ Batch {i//BATCH_SIZE + 1}: {len(records)} sequences fetched", "info")
                        handle.close()
                        break
                    except Exception as e:
                        if attempt == max_retries - 1:
                            print(f"❌ Failed batch {i//BATCH_SIZE + 1} after {max_retries} attempts: {e}")
                        else:
                            print(f"⚠️ Retrying batch {i//BATCH_SIZE + 1} (attempt {attempt + 1})...")
                            sleep(DELAY * 2)
                sleep(DELAY)

    except Exception as e:
        logger.log(f"Error while fetching neighbourhood protein sequences: {str(e)}", "error")
        raise

def run_makeblastdb(fasta_file, db_name, logger):
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
        logger.log(f"✅ BLAST database '{db_name}' created successfully!", "info")
    except subprocess.CalledProcessError as e:
        logger.log(f"❌ makeblastdb failed: {e}", "error")

def run_blastp(query_fasta, db_name, output_csv, threads, logger):
    """Run all-vs-all BLASTP with the original format 6 output plus headers."""
    headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    cmd = [
        'blastp',
        '-query', query_fasta,
        '-db', db_name,
        '-out', output_csv,
        '-outfmt', '6',
        '-evalue', '1e-5',
        '-num_threads', str(threads),
        '-max_hsps', '1'
    ]
    try:
        subprocess.run(cmd, check=True)
        
        with open(output_csv, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('\t'.join(headers) + '\n' + content)
        logger.log(f"✅ BLASTP results saved to {output_csv}", "info")
    except subprocess.CalledProcessError as e:
        logger.log(f"❌ blastp failed: {e}", "error")

def filter_high_identity_pairs(input_csv, output_csv, threshold, logger):
    """
    Filter BLAST results for high identity matches.
    """
    try:
        best_matches = defaultdict(float)
        
        with open(input_csv, 'r') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            for row in reader:
                qseqid = row['qseqid']
                sseqid = row['sseqid']
                pident = float(row['pident'])
                
                if qseqid == sseqid:
                    continue
                
                pair = tuple(sorted((qseqid, sseqid)))
                
                if pident > best_matches[pair]:
                    best_matches[pair] = pident
        
        with open(output_csv, 'w') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow(['Protein1', 'Protein2', 'Percent_Identity'])
            
            for (prot1, prot2), pident in best_matches.items():
                if pident >= threshold:
                    writer.writerow([prot1, prot2, pident])
        
        logger.log(f"✅ Filtered high-identity pairs saved to {output_csv}", "info")

    except Exception as e:
        logger.log(f"Error filtering high identity pairs: {str(e)}", "error")
        raise

def find_protein_clusters(input_file, output_file, min_identity, logger):
    """
    Groups interconnected proteins into clusters based on similarity relationships.
    """
    try:
        graph = defaultdict(set)
        
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if float(row['Percent_Identity']) >= min_identity:
                    prot1 = row['Protein1']
                    prot2 = row['Protein2']
                    graph[prot1].add(prot2)
                    graph[prot2].add(prot1)
        
        visited = set()
        clusters = []
        
        for protein in graph:
            if protein not in visited:
                cluster = set()
                stack = [protein]
                
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        cluster.add(current)
                        stack.extend(graph[current] - visited)
                
                if cluster:
                    clusters.append(cluster)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ClusterID', 'Proteins'])
            
            for i, cluster in enumerate(sorted(clusters, key=lambda x: -len(x)), 1):
                sorted_proteins = sorted(cluster)
                quoted_proteins = [f'{prot}' for prot in sorted_proteins]
                writer.writerow([f'Cluster_{i}', ', '.join(quoted_proteins)])
        
        logger.log(f"✅ Found {len(clusters)} protein clusters saved to {output_file}", "info")

    except Exception as e:
        logger.log(f"Protein Clustering failed: {str(e)}", "error")
        raise

def filter_fasta_by_clusters(cluster_csv, input_fasta, output_fasta, logger):
    """
    Remove proteins that appear after the first protein in each cluster.
    """
    proteins_to_remove = set()
    
    with open(cluster_csv, 'r') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            proteins = [p.strip('" ') for p in row[1].split(',')]
            proteins_to_remove.update(proteins[1:])
    
    logger.log(f"Found {len(proteins_to_remove)} proteins to remove from clusters", "debug")
    
    kept_records = []
    removed_count = 0
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in proteins_to_remove:
            removed_count += 1
        else:
            kept_records.append(record)
    
    with open(output_fasta, 'w') as f:
        SeqIO.write(kept_records, f, "fasta")
    
    logger.log(f"Removed {removed_count} proteins", "debug")
    logger.log(f"Kept {len(kept_records)} proteins in {output_fasta}","info")
    logger.log(f"Total proteins to remove that weren't found: {len(proteins_to_remove) - removed_count}", "warning")

def compare_hmm_accessions_with_description(psi_summary_path, neighbourhood_summary_path, output_dir, logger):
    """Compare two summary files and print common HMM profile accessions with their descriptions."""
    def extract_hmm_accession_description(filepath):
        accession_to_description = {}
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("Target Name") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    hmm_accession_raw = parts[1]
                    hmm_accession = hmm_accession_raw.split('.')[0]
                    description = parts[7].strip()
                    if hmm_accession not in accession_to_description:
                        accession_to_description[hmm_accession] = description
        return accession_to_description

    psi_data = extract_hmm_accession_description(psi_summary_path)
    neighbourhood_data = extract_hmm_accession_description(neighbourhood_summary_path)

    common_hmm_accessions = set(psi_data.keys()) & set(neighbourhood_data.keys())
    output_file = os.path.join(output_dir, 'common_hmm_accessions.txt')

    if not common_hmm_accessions:
        logger.log("No common HMM accessions found", "info")
        return False
    
    logger.log(f"Found {len(common_hmm_accessions)} common HMM accessions", "info")
    with open(output_file, 'w') as f_out:
        f_out.write("Common HMM accessions and their descriptions:\n\n")
        f_out.write("Accession\tDescription\n")
        f_out.write("----------------------------\n")

        for acc in sorted(common_hmm_accessions):
            description = psi_data.get(acc) or neighbourhood_data.get(acc)
            line = f"{acc}\t{description}"
            logger.log(line, "debug")
            f_out.write(line + "\n")
    
    logger.log(f"Results saved to: {output_file}", "info")
    return True

def delete_blastdb_files(db_name, logger):
    """Delete files created by makeblastdb."""
    extensions = ['.pin', '.phr', '.psq', '.pdb', '.pjs', '.pog', '.pos', '.pto', '.pot', '.ptf']
    for ext in extensions:
        file_path = f"{db_name}{ext}"
        if os.path.exists(file_path):
            os.remove(file_path)
            logger.log(f"🗑️ Deleted {file_path}", "info")
        else:
            logger.log(f"⚠️ File not found: {file_path}", "cyan")

class GeneFusePipeline:
    def __init__(self, args): #! <<<
        self.output_dir, self.subdirs = output_directory_initialisation(args.input_file, None)
        self.logger = PipelineLogger(self.output_dir)
        
        self.args = args
        self.logger.log("Initializing GeneFuse pipeline", "info")
        
        Entrez.email = args.user_email
        if args.api_key:
            Entrez.api_key = args.api_key
            self.logger.log("Using provided API key", "debug")

    def run(self):
        try:
            self.network_mgr = NetworkManager(self.logger)
            self._run_initial_steps()
            
            self.analyser = DataAnalyser(self.args, self.output_dir, self.subdirs, self.logger)
            self.analyser.run_analysis()
            self.analyser.cleanup()
        except Exception as e:
            self.logger.log(f"Pipeline failed: {str(e)}", "error")
            raise

    def _run_initial_steps(self):
        """Run PSI-BLAST and assembly downloads in parallel"""
        with multiprocessing.Pool(processes=2) as pool:
            psi_result = pool.apply_async(self._run_psiblast_and_process)
            assembly_result = pool.apply_async(self._download_assembly_data)
            
            psi_result.get()
            self.genbank_df, self.refseq_df = assembly_result.get()
            
        self._run_parallel_analyses()

    def _run_psiblast_and_process(self):
        """Run PSI-BLAST and process results"""
        self._run_psiblast()
        self._process_blast_results()

    def _download_assembly_data(self):
        """Download and load assembly summary files"""
        self.logger.log("Downloading assembly summary files...", "info")
        genbank_file, refseq_file = download_assembly_summary_files(self.logger)
        genbank_df = load_assembly_summary(genbank_file, self.logger)
        refseq_df = load_assembly_summary(refseq_file, self.logger)
        return genbank_df, refseq_df

    def _run_parallel_analyses(self):
        """Run domain and neighborhood analyses in parallel"""
        psi_fasta = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast.fasta")
        psi_ids = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast_ids.txt")

        processes = [
            multiprocessing.Process(
                target=domain_analysis,
                args=(self.args, self.output_dir, psi_fasta, self.subdirs['domains'], self.logger)
            ),
            multiprocessing.Process(
                target=neighborhood_analysis,
                args=(self.args, self.output_dir, psi_ids, self.subdirs['neighborhoods'], 
                      self.args.num_genes, self.logger, self.genbank_df, self.refseq_df)
            )
        ]

        for p in processes:
            p.start()
        for p in processes:
            p.join()

        self.logger.log("Analysis completed successfully", "info")

    def _run_psiblast(self):
        """Run PSI-BLAST (local operation)"""
        psi_output = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast.xml")
        pssm_output = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_PSSM.txt")
        
        run_psiblast(
            self.args.input_file,
            self.args.E_value_psi,
            self.args.iterations,
            self.args.database_psi,
            psi_output,
            self.args.num_threads,
            pssm_output,
            self.logger
        )

    def _process_blast_results(self):
        """Process BLAST results (local operation)"""
        psi_output = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast.xml")
        psi_csv = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast.csv")
        psi_fasta = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast.fasta")
        psi_ids = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast_ids.txt")

        blast_results, protein_ids, fasta_sequences = process_blast_xml(psi_output, self.logger)
        
        if blast_results:
            columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq"]
            pd.DataFrame(blast_results, columns=columns).to_csv(psi_csv, index=False)
            self.logger.log(f"Saved {len(blast_results)} BLAST hits to {psi_csv}", "info")
        
        if fasta_sequences:
            with open(psi_fasta, 'w') as f:
                f.write("\n".join(fasta_sequences) + "\n")
            self.logger.log(f"Saved {len(fasta_sequences)} sequences to {psi_fasta}", "info")
        
        if protein_ids:
            with open(psi_ids, 'w') as f:
                f.write("\n".join(protein_ids) + "\n")
            self.logger.log(f"Saved {len(protein_ids)} protein IDs to {psi_ids}", "info")

class DataAnalyser:
    """Handles post-processing and comparative analysis of pipeline results"""
    
    def __init__(self, args, output_dir, subdirs, logger):
        self.args = args
        self.output_dir = output_dir
        self.subdirs = subdirs
        self.logger = logger
        
        # Constants
        self.BATCH_SIZE = 2000
        self.DELAY = 0.34  # NCBI recommended delay
        
        # Setup directories
        self.domain_analysis_dir = os.path.join(output_dir, "Domain_analysis")
        self.common_hmm_dir = os.path.join(self.domain_analysis_dir, "Common_HMMs")
        self.matches_dir = os.path.join(self.common_hmm_dir, "matches")
        
        os.makedirs(self.domain_analysis_dir, exist_ok=True)
        os.makedirs(self.common_hmm_dir, exist_ok=True)
        os.makedirs(self.matches_dir, exist_ok=True)

    def run_analysis(self):
        """Execute the complete analysis workflow"""
        try:
            self.logger.log("Starting comprehensive data analysis", "info")
            
            # Step 1: Run Rosetta Stone analysis (neighborhood domain analysis)
            self._run_rosetta_stone_analysis()
            
            # Step 2: Compare HMM accessions between PSI and neighborhood results
            self._compare_hmm_accessions()
            
            # Step 3: Run BLAST analysis on common HMMs
            self._run_common_hmm_blast()
            
            # Step 4: Analyze matches between PSI and neighborhood proteins
            self._analyze_matches()
            
            self.logger.log("Data analysis completed successfully", "info")
            return True
            
        except Exception as e:
            self.logger.log(f"Data analysis failed: {str(e)}", "error")
            raise

    def _run_rosetta_stone_analysis(self):
        """Run the equivalent of rosetta_stone.py"""
        try:
            self.logger.log("Running Rosetta Stone analysis (neighborhood domain analysis)", "info")
            
            # Input files
            neighborhood_csv = os.path.join(self.subdirs['neighborhoods'], "txt_files", "combined_cumulative.txt")
            psi_id_file = os.path.join(self.subdirs['psi_results'], f"{os.path.basename(self.output_dir)}_Psi_blast_ids.txt")
            
            # Output files
            fasta_file = os.path.join(self.domain_analysis_dir, 'neighborhood_sequences.fasta')
            output_tblout = os.path.join(self.domain_analysis_dir, 'hmmscan_results.tblout')
            output_summary = os.path.join(self.domain_analysis_dir, "domain_summary.txt")
            output_cumulative = os.path.join(self.domain_analysis_dir, "cumulative_summary.txt")
            
            # Step 1: Extract protein IDs from neighborhood results
            protein_accessions = self._extract_neighborhood_proteins(neighborhood_csv)
            
            if not protein_accessions:
                raise ValueError("No protein accessions found in neighborhood results")
            
            # Step 2: Fetch sequences for these accessions
            self._fetch_sequences(protein_accessions, fasta_file)
            
            # Step 3: Run domain analysis on neighborhood sequences
            run_hmmscan(fasta_file, output_tblout, self.args.database_hmm, self.logger)
            domain_hits = parse_hmmscan_output(output_tblout, self.logger)
            
            # Generate outputs
            write_domain_hits(domain_hits, output_summary, self.logger)
            summarize_results(output_summary, output_cumulative, fasta_file, self.logger)
            
            self.logger.log("Rosetta Stone analysis completed", "info")
            
        except Exception as e:
            self.logger.log(f"Rosetta Stone analysis failed: {str(e)}", "error")
            raise

    def _compare_hmm_accessions(self):
        """Compare HMM accessions between PSI and neighborhood results and create a proper CSV"""
        try:
            self.logger.log("Comparing HMM accessions between PSI and neighborhood results", "info")
            
            # Input files
            psi_summary = os.path.join(self.subdirs['domains'], f"{os.path.basename(self.output_dir)}_summary.txt")
            neigh_summary = os.path.join(self.domain_analysis_dir, "domain_summary.txt")
            
            # Output file
            output_file = os.path.join(self.common_hmm_dir, 'common_hmm_accessions.csv')
            
            # Parse both summary files properly
            psi_data = self._parse_domain_summary(psi_summary)
            neigh_data = self._parse_domain_summary(neigh_summary)
            
            # Find common HMM accessions
            common_hmms = set(psi_data.keys()) & set(neigh_data.keys())
            
            if not common_hmms:
                self.logger.log("No common HMM accessions found between PSI and neighborhood results", "warning")
                return False
            
            # Prepare the CSV data
            csv_data = []
            for hmm_acc in common_hmms:
                # Get counts and descriptions
                psi_count = psi_data[hmm_acc]['count']
                neigh_count = neigh_data[hmm_acc]['count']
                
                # Use PSI description if available, otherwise neighborhood description
                description = psi_data[hmm_acc].get('description', 
                                                neigh_data[hmm_acc].get('description', 'Unknown'))
                
                csv_data.append({
                    'Accession': hmm_acc,
                    'Description': description,
                    'Homolog_proteins': psi_count,
                    'Neighbour_proteins': neigh_count
                })
            
            # Sort by total occurrences (descending)
            csv_data.sort(key=lambda x: -(x['Homolog_proteins'] + x['Neighbour_proteins']))
            
            # Write to CSV
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['Accession', 'Description', 'Homolog_proteins', 'Neighbour_proteins'])
                writer.writeheader()
                writer.writerows(csv_data)
            
            self.logger.log(f"Saved {len(csv_data)} common HMM accessions to {output_file}", "info")
            return True
            
        except Exception as e:
            self.logger.log(f"Error comparing HMM accessions: {str(e)}", "error")
            raise

    def _parse_domain_summary(self, summary_file):
        """Parse domain summary file into a structured format"""
        hmm_data = {}
        
        try:
            with open(summary_file, 'r') as f:
                # Skip header line
                next(f)
                
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        target_name = parts[0].strip()
                        accession = parts[1].strip().split('.')[0]  # Remove version number
                        query_name = parts[2].strip()
                        description = parts[7].strip()
                        
                        if accession not in hmm_data:
                            hmm_data[accession] = {
                                'description': description,
                                'count': 0,
                                'proteins': set()
                            }
                        
                        # Add the protein ID (remove version if present)
                        protein_id = query_name.split('.')[0]
                        if protein_id not in hmm_data[accession]['proteins']:
                            hmm_data[accession]['proteins'].add(protein_id)
                            hmm_data[accession]['count'] += 1
                        
        except Exception as e:
            self.logger.log(f"Error parsing domain summary file {summary_file}: {str(e)}", "error")
            raise
        
        return hmm_data

    def _parse_hmm_summary(self, summary_file):
        """Parse HMM summary file into a dictionary {accession: {'description': str, 'count': int}}"""
        data = {}
        with open(summary_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    accession = parts[1].split('.')[0]  # Remove version number
                    description = parts[-1]  # Description is last column
                    count = 1  # Each line represents one occurrence
                    
                    if accession in data:
                        data[accession]['count'] += 1
                    else:
                        data[accession] = {
                            'description': description,
                            'count': count
                        }
        return data

    def _run_common_hmm_blast(self):
        try:
            self.logger.log("Running BLAST analysis on common HMM proteins", "info")
            
            # Input files
            common_csv = os.path.join(self.common_hmm_dir, 'common_hmm_accessions.csv')
            psi_summary = os.path.join(self.subdirs['domains'], f"{os.path.basename(self.output_dir)}_summary.txt")
            neigh_summary = os.path.join(self.domain_analysis_dir, "domain_summary.txt")
            
            # Output files
            psi_fasta = os.path.join(self.common_hmm_dir, "psi_proteins.fasta")
            neigh_fasta = os.path.join(self.common_hmm_dir, "neighbour_proteins.fasta")
            blast_output = os.path.join(self.common_hmm_dir, "blast_results.txt")
            high_id_output = os.path.join(self.common_hmm_dir, 'filtered_pairs.tsv')
            cluster_file = os.path.join(self.common_hmm_dir, 'protein_clusters.csv')
            
            # Step 1: Load common HMMs
            common_hmms = self._load_common_hmms(common_csv)
            if not common_hmms:
                self.logger.log("No common HMM accessions found - skipping BLAST analysis", "warning")
                return True
            
            # Step 2: Extract protein accessions with proper validation
            psi_prots = self._extract_proteins_from_domains(psi_summary, common_hmms)
            neigh_prots = self._extract_proteins_from_domains(neigh_summary, common_hmms)
            
            self.logger.log(f"Found {len(psi_prots)} PSI proteins and {len(neigh_prots)} neighbor proteins", "info")
            
            if not psi_prots or not neigh_prots:
                self.logger.log("Insufficient proteins for BLAST analysis - skipping", "warning")
                return True
            
            # Step 3: Download sequences with validation
            if not self._fetch_sequences_with_validation(psi_prots, psi_fasta, "PSI"):
                return True
            if not self._fetch_sequences_with_validation(neigh_prots, neigh_fasta, "neighbor"):
                return True
            
            # Step 4: Create BLAST database and run BLASTP
            db_name = os.path.join(self.common_hmm_dir, "neighbour_db")
            self._make_blast_db(neigh_fasta, db_name)
            self._run_blastp(psi_fasta, db_name, blast_output)
            
            # Step 5: Filter and cluster results
            self._filter_high_identity_pairs(blast_output, high_id_output, self.args.identity_threshold)
            self._find_protein_clusters(high_id_output, cluster_file, self.args.identity_threshold)
            
            # Clean up
            self._delete_blastdb_files(db_name)
            
            self.logger.log("Common HMM BLAST analysis completed", "info")
            return True
            
        except Exception as e:
            self.logger.log(f"Common HMM BLAST analysis failed: {str(e)}", "error")
            raise

    def _extract_proteins_from_domains(self, summary_file, common_hmms):
        """Extract protein IDs from domain summary file that match common HMMs"""
        proteins = set()
        try:
            with open(summary_file, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        hmm_acc = parts[1].split('.')[0]  # Get base accession without version
                        protein_id = parts[2].strip()
                        if hmm_acc in common_hmms and protein_id and protein_id != "Unknown":
                            proteins.add(protein_id.split('.')[0])  # Store base protein ID
        except Exception as e:
            self.logger.log(f"Error extracting proteins from {summary_file}: {str(e)}", "error")
            raise
        return list(proteins)

    def _fetch_sequences_with_validation(self, protein_ids, output_fasta, protein_type):
        """Fetch sequences with validation of the input IDs"""
        if not protein_ids:
            self.logger.log(f"No {protein_type} protein IDs to fetch", "warning")
            return False
        
        try:
            # First check if we have valid protein IDs
            valid_ids = [pid for pid in protein_ids if pid and len(pid) > 3]  # Basic validation
            if not valid_ids:
                self.logger.log(f"No valid {protein_type} protein IDs to fetch", "warning")
                return False
            
            # Then fetch the sequences
            self._fetch_sequences(valid_ids, output_fasta)
            
            # Verify the output file was created and has content
            if not os.path.exists(output_fasta) or os.path.getsize(output_fasta) == 0:
                self.logger.log(f"Failed to fetch {protein_type} protein sequences - empty output file", "error")
                return False
                
            return True
        except Exception as e:
            self.logger.log(f"Error fetching {protein_type} protein sequences: {str(e)}", "error")
            return False

    def _analyze_matches(self):
        """Run the equivalent of matches.py"""
        try:
            self.logger.log("Analyzing matches between PSI and neighborhood proteins", "info")
            
            # Input files
            csv_directory = os.path.join(self.subdirs['neighborhoods'], "csv_files")
            blast_results = os.path.join(self.common_hmm_dir, "blast_results.txt")
            psi_fasta = os.path.join(self.common_hmm_dir, "psi_proteins.fasta")
            neigh_fasta = os.path.join(self.common_hmm_dir, "neighbour_proteins.fasta")
            psi_domains_path = os.path.join(self.subdirs['domains'], f"{os.path.basename(self.output_dir)}_summary.txt")
            neigh_domains_path = os.path.join(self.domain_analysis_dir, "domain_summary.txt")
            
            # Get all data
            all_pairs = self._extract_all_pairs(blast_results)
            self.logger.log(f"Found {len(all_pairs)} protein pairs total", "info")
            
            fasta1_ids = self._extract_ids_from_fasta(neigh_fasta)
            fasta2_ids = self._extract_ids_from_fasta(psi_fasta)
            
            # Load domain information
            neighbour_domains = self._parse_domain_summary(neigh_domains_path)
            psi_domains = self._parse_domain_summary(psi_domains_path)
            
            # Build genome ID mapping
            protein_to_genome = self._build_genome_id_map(csv_directory)
            
            # Categorize matches by identity ranges
            categories = self._categorize_matches(
                all_pairs, fasta1_ids, fasta2_ids, 
                neighbour_domains, psi_domains, protein_to_genome
            )
            
            # Write separate CSV files for each category
            for cat_name, cat_data in categories.items():
                self._write_category_csv(cat_name, cat_data['matches'])
            
            self.logger.log("Match analysis completed", "info")
            
        except Exception as e:
            self.logger.log(f"Match analysis failed: {str(e)}", "error")
            raise

    # Helper methods for _run_rosetta_stone_analysis
    def _extract_neighborhood_proteins(self, neighborhood_file):
        """Extract unique protein accessions from neighborhood results"""
        protein_accessions = set()
        
        with open(neighborhood_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # Skip header
            for row in reader:
                if len(row) > 0:
                    protein_id = row[0].strip()
                    if protein_id and protein_id != "Unknown":
                        # Remove version number if present (e.g., WP_123456789.1)
                        base_id = protein_id.split('.')[0]
                        protein_accessions.add(base_id)
        
        self.logger.log(f"Extracted {len(protein_accessions)} protein accessions from neighborhood results", "info")
        return list(protein_accessions)

    # Helper methods for _compare_hmm_accessions
    def _load_accessions_with_counts(self, filepath):
        """Loads a file and returns a dict: {accession: (description, count)}"""
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

    # Helper methods for _run_common_hmm_blast
    def _load_common_hmms(self, csv_path):
        common_hmms = set()
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                common_hmms.add(row["Accession"].strip())
        return common_hmms

    def _extract_accessions_from_summary(self, summary_file, common_hmms):
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

    def _make_blast_db(self, fasta_file, db_name):
        subprocess.run([
            "makeblastdb", "-in", fasta_file, 
            "-dbtype", "prot", "-out", db_name
        ], check=True)

    def _run_blastp(self, query_fasta, db_name, output_file):
        subprocess.run([
            "blastp", "-query", query_fasta, "-db", db_name,
            "-out", output_file, "-outfmt", "6 qseqid sseqid pident",
            "-evalue", "1e-5", "-num_threads", str(self.args.num_threads)
        ], check=True)

    def _filter_high_identity_pairs(self, input_csv, output_csv, threshold):
        best_matches = defaultdict(float)
        with open(input_csv, 'r') as infile:
            reader = csv.DictReader(infile, delimiter='\t', fieldnames=["qseqid", "sseqid", "pident"])
            for row in reader:
                qseqid = row['qseqid']
                sseqid = row['sseqid']
                pident = float(row['pident'])
                if qseqid == sseqid:
                    continue
                pair = tuple(sorted((qseqid, sseqid)))
                if pident > best_matches[pair]:
                    best_matches[pair] = pident
        
        with open(output_csv, 'w') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow(['Protein1', 'Protein2', 'Percent_Identity'])
            for (prot1, prot2), pident in best_matches.items():
                if pident >= threshold:
                    writer.writerow([prot1, prot2, pident])
        
        self.logger.log(f"Filtered high-identity pairs saved to {output_csv}", "info")

    def _find_protein_clusters(self, input_file, output_file, min_identity):
        graph = defaultdict(set)
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if float(row['Percent_Identity']) >= min_identity:
                    prot1 = row['Protein1']
                    prot2 = row['Protein2']
                    graph[prot1].add(prot2)
                    graph[prot2].add(prot1)
        
        visited = set()
        clusters = []
        for protein in graph:
            if protein not in visited:
                cluster = set()
                stack = [protein]
                while stack:
                    current = stack.pop()
                    if current not in visited:
                        visited.add(current)
                        cluster.add(current)
                        stack.extend(graph[current] - visited)
                if cluster:
                    clusters.append(cluster)
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ClusterID', 'Proteins'])
            for i, cluster in enumerate(sorted(clusters, key=lambda x: -len(x)), 1):
                sorted_proteins = sorted(cluster)
                writer.writerow([f'Cluster_{i}', ', '.join(sorted_proteins)])
        
        self.logger.log(f"Found {len(clusters)} protein clusters saved to {output_file}", "info")

    def _delete_blastdb_files(self, db_name):
        extensions = ['.pin', '.phr', '.psq', '.pdb', '.pjs', '.pog', '.pos', '.pto', '.pot', '.ptf']
        for ext in extensions:
            file_path = f"{db_name}{ext}"
            if os.path.exists(file_path):
                os.remove(file_path)
                self.logger.log(f"Deleted {file_path}", "debug")

    # Helper methods for _analyze_matches
    def _extract_all_pairs(self, txt_file):
        """Extract all protein pairs with their identity scores."""
        pairs = []
        with open(txt_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    pid1, pid2, identity = parts
                    try:
                        pairs.append((pid1.strip(), pid2.strip(), float(identity)))
                    except ValueError:
                        continue  # Skip invalid float conversion
        return pairs

    def _extract_ids_from_fasta(self, fasta_file):
        """Extract all protein IDs from a FASTA file."""
        return {record.id.strip() for record in SeqIO.parse(fasta_file, 'fasta')}

    def _parse_domain_summary(self, summary_file):
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

    def _build_genome_id_map(self, csv_directory):
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

    def _categorize_matches(self, pairs, fasta1_ids, fasta2_ids, neighbour_domains, psi_domains, protein_to_genome):
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
                    'neighbour_domains': self._format_domains(self._get_domains(neighbour, neighbour_domains)),
                    'psi_domains': self._format_domains(self._get_domains(psi, psi_domains))
                }
                
                # Categorize the match
                for cat_name, cat in categories.items():
                    if cat['min'] < identity <= cat['max']:
                        cat['matches'].append(match_data)
                        break
        
        return categories

    def _get_domains(self, protein_id, domain_map):
        """Get domains for a protein ID, handling version numbers."""
        base_id = protein_id.split('.')[0]
        return domain_map.get(base_id, [])

    def _format_domains(self, domain_list):
        """Format list of domains into a semicolon-separated string."""
        return "; ".join(domain_list) if domain_list else "NA"

    def _write_category_csv(self, category_name, matches):
        """Write matches for a specific category to CSV."""
        if not matches:
            self.logger.log(f"No matches found in category {category_name}", "warning")
            return
        
        output_csv = os.path.join(self.matches_dir, f"matches_{category_name}.csv")
        fieldnames = ['psi_protein', 'neighbour_protein', 'percent_identity', 
                     'genome_id', 'neighbour_domains', 'psi_domains']
        
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(matches)
        
        self.logger.log(f"{len(matches)} sequences written to {output_csv}", "info")

    # Shared helper methods
    def _fetch_sequences(self, protein_ids, output_fasta):
        """Fetch protein sequences from NCBI in batches"""
        with open(output_fasta, 'w') as fasta_out:
            for i in range(0, len(protein_ids), self.BATCH_SIZE):
                batch = protein_ids[i:i + self.BATCH_SIZE]
                attempt = 1
                max_attempts = 3
                
                while attempt <= max_attempts:
                    try:
                        handle = Entrez.efetch(
                            db='protein',
                            id=','.join(batch),
                            rettype='fasta',
                            retmode='text'
                        )
                        records = list(SeqIO.parse(handle, 'fasta'))
                        SeqIO.write(records, fasta_out, 'fasta')
                        self.logger.log(f"Fetched batch {i//self.BATCH_SIZE + 1}: {len(records)} sequences", "info")
                        handle.close()
                        break
                    except Exception as e:
                        if attempt == max_attempts:
                            self.logger.log(f"Failed batch {i//self.BATCH_SIZE + 1} after {max_attempts} attempts: {e}", "error")
                            raise
                        self.logger.log(f"Retrying batch {i//self.BATCH_SIZE + 1} (attempt {attempt})", "warning")
                        time.sleep(5 * attempt)  # Exponential backoff
                        attempt += 1
                
                time.sleep(self.DELAY)  # Respect NCBI rate limits
        
        return output_fasta

    def cleanup(self):
        """Clean up temporary files"""
        try:
            # Clean up any temporary BLAST databases
            db_files = glob.glob(os.path.join(self.common_hmm_dir, 'blast_db*'))
            for db_file in db_files:
                try:
                    os.remove(db_file)
                except Exception as e:
                    self.logger.log(f"Error removing {db_file}: {str(e)}", "warning")
            
            self.logger.log("Cleanup completed", "info")
        except Exception as e:
            self.logger.log(f"Error during cleanup: {str(e)}", "warning")

def main():
    """Main execution function"""
    try:
        args = parse_arguments()
        pipeline = GeneFusePipeline(args)
        pipeline.run()
        print("✨ Analysis complete!")
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        raise

if __name__ == "__main__":
    print("""
   ██████╗  ███████╗ ███╗   ██╗ ███████╗ ███████╗ ██╗   ██╗ ███████╗ ███████╗
  ██╔════╝  ██╔════╝ ████╗  ██║ ██╔════╝ ██╔════╝ ██║   ██║ ██╔════╝ ██╔════╝
  ██║  ███╗ █████╗   ██╔██╗ ██║ █████╗   █████╗   ██║   ██║ ███████╗ █████╗
  ██║   ██║ ██╔══╝   ██║╚██╗██║ ██╔══╝   ██╔══╝   ██║   ██║ ╚════██║ ██╔══╝
  ╚██████╔╝ ███████╗ ██║ ╚████║ ███████╗ ██║      ╚██████╔╝ ███████║ ███████╗
   ╚═════╝  ╚══════╝ ╚═╝  ╚═══╝ ╚══════╝ ╚═╝       ╚═════╝  ╚══════╝ ╚══════╝
    """)
    print("GeneFuse - Genomic Neighborhood Analysis Pipeline\n")
    main()
  
