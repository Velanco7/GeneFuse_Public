# 🧬 GeneFuse

**GeneFuse** is a bioinformatics tool under active development that identifies and analyses potential gene fusion events. The pipeline integrates protein domain analysis, genomic neighbourhood exploration, and sequence similarity to help infer the function of uncharacterized proteins or fusion events.

> **Status**: In development 🚧  
> Author: Velanco Fernandes

---

## 🚀 Features

- Genomic neighbourhood and protein clustering
- Integration of domain information via HMMER and PSI-BLAST
- Rosetta-style "fusion logic" to flag potential fused proteins
- Customizable, modular architecture for further expansion

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `Complete_Script_V1.2.1.py` | Main script integrating multiple modules for fusion detection |
| `Complete_rosetta.py` | Rosetta-style analysis logic used to identify potential protein fusions |
| `README.md` | Project overview and setup guide |
| `requirements.txt` | Python dependencies for running GeneFuse (to be added) |

---

## 🛠️ Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Velanco7/GeneFuse_shared.git
   cd GeneFuse_shared
   ```

2. (Optional) Create a virtual environment:
   ```bash
   python -m venv env
   source env/bin/activate  # or env\Scripts\activate on Windows
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---

## 🧪 Usage

Run the main script:
```bash
python Complete_Script_V1.2.1.py --input example.fasta
```

> Add flags and input specifications as needed once CLI is implemented.

---

## ✅ To-Do

- [ ] Integrate `Complete_rosetta.py` fully into a class/module structure
- [ ] Add CLI support using `argparse` or `typer`
- [ ] Add test dataset and output examples
- [ ] Write unit tests
- [ ] Add Docker support

---

## 📄 License

To be added (suggested: MIT or GPL).

---

## 🤝 Contributions

This project is personal, but collaborators with relevant bioinformatics experience are welcome to fork, suggest features, or raise issues.
