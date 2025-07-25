# 🧬 GeneFuse

**GeneFuse** is a bioinformatics tool (currently optimised for Linux) that identifies and analyses potential gene fusion events. The pipeline integrates protein domain analysis, genomic neighbourhood exploration, and sequence similarity to help infer the function of uncharacterized proteins or fusion events.

> **Version**: 1.3.6  
> **Author**: Velanco Fernandes   
> **License** This project is licensed under the GNU General Public License v3 (GPLv3). See [LICENSE.txt](LICENSE.txt) for details.  

---

## 🚀 Features

- Genomic neighbourhood and protein clustering
- Integration of domain information via HMMER and PSI-BLAST
- Rosetta-style "fusion logic" to flag potential fused proteins (complete_rosetta.py)
- Customisable, modular architecture for further expansion

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `Complete_Script_V1.2.1.py` | Main script integrating multiple modules for fusion detection |
| `Complete_rosetta.py` | Rosetta-style analysis logic used to identify potential protein fusions |
| `README.md` | Project overview and setup guide |
| `requirements.txt` | Python dependencies for running GeneFuse (to be added) |

---

## Example Datasets

We include two real-world analyses to help you get started:

**GeneFuse_PRTase_20250522_134512** 
  - Focus: Phosphoribosyltransferase (PRTase) gene neighbourhood analysis 
  - Method: Uses legacy complete_rosetta.py (deprecated standalone script) 
  - Use Case: Demonstrates early pipeline capabilities (retained for archival purposes) 

**Place Holder** 
   - Focus: ************ fusion detection 
   - Method: Uses modern DataAnalyser class (fully integrated in complete_script_V1.3.6.py) 
   - Key Improvements: 
     - Unified error handling 
     - Optimised clustering algorithms 
     - Automated report generation 

---

## 🛠️ Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Velanco7/GeneFuse_Public.git
   cd GeneFuse
   ```

2. (Optional) Create a virtual environment:
   ```bash
   python -m venv env
   source env/bin/activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
# External Tools

The following external tools are not Python packages and must be installed separately on your system:

BLAST+:
Download and install from NCBI BLAST executables

HMMER (hmmscan):
Download and install from HMMER official site

Databases:
RefSeq protein database and Pfam HMM databases are data files that need to be downloaded independently. Make sure to provide the correct paths to these files when running GeneFuse.


---

## 🧪 Usage


📦 GeneFuse – Command Line Usage
Run Complete_Script_V1.2.1.py with the following arguments:

```bash
python Complete_Script_V1.2.1.py \
  input.fasta \
  -u your_email@example.com \
  -db1 path/to/psiblast_db \
  -db2 path/to/Pfam-A.hmm path/to/another_db.hmm \
  -g 5 \
  -e1 0.00001 \
  -iterations 3 \
  -num_threads 4 \
  -api_key your_entrez_api_key \
  -c config.ini \
  --identity_threshold 95.0
```

🧾 Required Arguments
| Argument                 | Description                                    |
| ------------------------ | ---------------------------------------------- |
| `input_file`             | FASTA file of your query protein               |
| `-u`, `--user_email`     | Email address (used with Entrez queries)       |
| `-db1`, `--database_psi` | Path to PSI-BLAST database                     |
| `-db2`, `--database_hmm` | One or more HMMER databases (e.g., Pfam-A.hmm) |


⚙️ Optional Arguments
| Argument               | Default    | Description                                    |
| ---------------------- | ---------- | ---------------------------------------------- |
| `-g`, `--num_genes`    | 5          | Number of upstream/downstream genes to extract |
| `-e1`, `--E_value_psi` | 0.00001    | E-value cutoff for PSI-BLAST                   |
| `--iterations`         | 3          | Number of PSI-BLAST iterations                 |
| `--num_threads`        | 4          | Threads for processing                         |
| `--api_key`            | None       | Optional Entrez API key                        |
| `-c`, `--config`       | config.ini | Path to configuration file                     |
| `--identity_threshold` | 95.0       | Minimum percent identity for filtering         |


> Add flags and input specifications as needed once CLI is implemented.

---

## ✅ To-Do

- [ ] Improve the Data Analyser class with further clustering, and matching of the final.csv and reintroduced.csv
- [ ] Add CLI support using `argparse` or `typer`
- [ ] Write unit tests
- [ ] Add Docker support


---


##🤝 Contributions
This is a personal research tool, but contributions are welcome — particularly from those with experience in:

- Protein domain analysis
- Comparative genomics
- Functional annotation pipelines
- Python and bioinformatics tooling

Feel free to fork the repository, open issues, or suggest improvements!


