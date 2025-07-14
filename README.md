from pathlib import Path

readme_text = """
#  Research Paper Fetcher (Command-Line Tool)

This is a beginner-friendly Python CLI tool that helps you **search research papers from PubMed** and filters them to show only those with **authors affiliated with pharmaceutical or biotech companies** (like Pfizer, Roche, AstraZeneca, etc.).

---

##  What Is This Tool?

This tool connects to [PubMed](https://pubmed.ncbi.nlm.nih.gov/), a free search engine for life science and biomedical research papers.

Instead of showing *all* papers, this tool:
- Finds authors who **work in companies**, not just universities or hospitals
- Filters and shows only those research papers
- Saves the filtered results into a **CSV file (Excel-friendly)**

---

##  Why Should You Use This?

Let’s say you want to:
- Research **industry involvement** in cancer, diabetes, or COVID
- Track what **pharma companies** are working on
- Avoid manual searching and filtering on the PubMed website

This tool makes it super easy. Just run a command, and it gives you a ready-to-use CSV with filtered results.

---

##  Features

-  Search any PubMed topic (e.g., "cancer AND 2023[dp]")
-  Detect non-academic (company) affiliations
-  Save results as `.csv` (Excel compatible)
-  CLI-based and beginner-friendly
-  Built with Biopython, Pandas, Poetry

---

##  Project Structure

Research_paper_fetcher/
├── Research_paper_fetcher/
│   ├── __init__.py
│   ├── pubmed.py         # Core logic to query and filter papers
│   └── cli.py            # Command-line interface
├── tests/
│   └── test_pubmed.py    # (Optional) Unit tests
├── pyproject.toml        # Poetry config
└── README.md             # You're reading this

---

##  Setup (Step-by-Step)

This project uses **Poetry** for dependency management.

### Step 1: Install Poetry

    pip install poetry

### Step 2: Clone the Repository

    git clone https://https://github.com/Prajwalkulkarni2003/Research_paper_fetcher.git
    cd Research-paper-fetcher

Or download the ZIP file and extract it.

### Step 3: Install Dependencies

    poetry install

This installs:
- `biopython` → to access the PubMed API
- `pandas` → to generate the CSV file

### Step 4: Set Your Email

Open `pubmed.py` and update this line with your real email:

    Entrez.email = "your_email@example.com"

PubMed requires this for API usage tracking.

---

##  How to Use the Tool

### Basic Usage

    poetry run get-papers-list "YOUR SEARCH TERM" -f output.csv

### Example:

    poetry run get-papers-list "cancer AND 2023[dp]" -f cancer_results.csv

This:
- Searches cancer papers from 2023
- Filters only papers with company-affiliated authors
- Saves results to `cancer_results.csv`

### Debug Mode

Print results to terminal (instead of saving):

    poetry run get-papers-list "diabetes AND insulin" -d

---

##  Example CSV Output

| PubMed ID | Title                  | Publication Date | Company Affiliation   | Email              |
|-----------|------------------------|------------------|------------------------|--------------------|
| 12345678  | New Drug for Cancer    | 2023             | Pfizer Inc             | john@pfizer.com    |
| 12345679  | Biotech Vaccine Study  | 2023             | Genentech Ltd          | dr.smith@genetech.com |

---

##  What Company Names Are Detected?

The tool filters based on keywords found in the affiliation line:

**Included if it contains:**  
"Inc", "Ltd", "Corp", "Pharma", "Biotech", "Pfizer", "Roche", etc.

**Excluded if it contains:**  
"University", "Hospital", "School", "Institute", "Lab", "College"

You can edit the list in `pubmed.py`:

    NON_ACADEMIC_KEYWORDS = ["Pharma", "Biotech", "Inc", "Ltd", "Corp"]
    ACADEMIC_KEYWORDS = ["University", "Institute", "Hospital", "School", "College", "Lab"]

---

##  FAQ

**Can I use this without knowing Python?**  
 Yes! Just copy-paste commands in your terminal.

**Will it work on Windows?**  
 Yes. Works on Windows, macOS, and Linux (Python + Poetry required).

**Does it include university papers?**  
 No. It only includes papers with **company-affiliated** authors.

**Can I add more filters or export to Excel?**  
 Yes! You can extend it using `pandas` to filter more, or save as `.xlsx`.

---

##  Example Commands

    poetry run get-papers-list "diabetes AND 2023[dp]" -f diabetes.csv
    poetry run get-papers-list "covid AND 2022[dp]" -d

---

##  How It Works (In Simple Terms)

1. You enter a PubMed search query  
2. The tool sends this query to PubMed using Biopython  
3. It fetches the list of matching papers  
4. It checks each author’s affiliation  
5. If any author is from a company, it saves the paper  
6. You get a filtered CSV file with results

---

##  Technologies Used

- BioPython – PubMed API access
- Pandas – Data processing
- Poetry – Dependency management

---

##  Author

**Prajwal R K**  
Built for educational, research, and technical demonstration purposes.

---

##  License

Free to use for research and educational purposes.  
Follow PubMed’s API guidelines and provide a valid email in all requests.
"""

# Save to text file
readme_path = Path("/mnt/data/README_PubMedFetcher.txt")
readme_path.write_text(readme_text.strip())

readme_path

