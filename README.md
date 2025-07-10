#  Research_paper_fetche(Command Line Tool)

This is a beginner-friendly command-line Python tool that helps you **search research papers from PubMed** and filters them to show only those with **authors affiliated with pharmaceutical or biotech companies** (like Pfizer, Roche, AstraZeneca, etc.).

---

##  What Is This Tool?

This tool connects to [PubMed](https://pubmed.ncbi.nlm.nih.gov/), a free search engine for life science and biomedical research papers.

Instead of showing **all papers**, this tool:
- Finds authors who **work in companies**, not just universities or hospitals
- Filters and shows only those research papers
- Saves the filtered results into a **CSV file (Excel-style)** for easy reading

---

##  Why Should You Use This?

Let’s say you want to:
- Research **industry involvement** in cancer, diabetes, or COVID
- Track what **pharma companies** are working on
- Avoid manual searching and filtering on the PubMed website

 This tool makes it super easy. Just type one command, and it gives you a ready-to-use CSV file with clean results.

---

##  Features

-  Search any PubMed topic (e.g., `"cancer AND 2023[dp]"`)
-  Filters authors who work in **non-academic companies**
-  Exports results as `.csv` (Excel compatible)
-  Command-line friendly, fast, and beginner-safe

---

##  Project Structure

Research_paper_fetcher/
├── Research_paper_fetche/
│ ├── init.py
│ ├── pubmed.py ← Logic to query PubMed and filter company authors
│ └── cli.py ← Command-line interface (CLI)
├── tests/
│ └── test_pubmed.py ← (Optional) Test file
├── pyproject.toml ← Poetry config file (dependencies, commands)
└── README.md ← You are here!

yaml
Copy
Edit

---

##  How to Set It Up (Step by Step)

This project uses **[Poetry](https://python-poetry.org/)** — a tool that manages Python projects easily.

###  Step 1: Install Poetry

If not installed, run:

```bash
pip install poetry
Step 2: Download or Clone the Project
If using Git:

bash
Copy
Edit
git clone https://github.com/yourusername/pubmed-company-paper-fetcher.git
cd pubmed-company-paper-fetcher
Or download the ZIP file and extract it.

Step 3: Install Dependencies
Inside the project folder:

bash
Copy
Edit
poetry install
This installs required libraries:

biopython → used to connect with PubMed API

pandas → used to save data into CSV

Step 4: Set Your Email (required for PubMed access)
Open the file pubmed.py and update this line:

python
Copy
Edit
Entrez.email = "your_email@example.com"
Change it to your real email. PubMed uses it for API usage tracking.

How to Use the Tool
Run the tool from your terminal like this:

bash
Copy
Edit
poetry run get-papers-list "YOUR SEARCH TERM" -f output.csv
Example:
bash
Copy
Edit
poetry run get-papers-list "cancer AND 2023[dp]" -f cancer_results.csv
This will:

Search for cancer-related papers published in 2023

Filter only papers where at least one author is from a company

Save the results to cancer_results.csv

Debug Mode (Optional)
If you want to see results in the terminal instead of saving to file:

bash
Copy
Edit
poetry run get-papers-list "diabetes AND insulin" -d
The -d or --debug option prints all filtered records to the screen.

Example CSV Output
After running, you’ll get a file like this:

PubMed ID	Title	Publication Date	Company Affiliation	Email
12345678	New Drug for Cancer	2023	Pfizer Inc	john@pfizer.com
12345679	Biotech Vaccine Study	2023	Genentech Ltd	dr.smith@gentech.com

 What Company Names Are Detected?
This tool checks for keywords in the author's affiliation such as:

"Inc", "Ltd", "Corp", "Pharma", "Biotech"

"Pfizer", "Roche", "AstraZeneca", etc.

You can edit this list inside pubmed.py under:

python
Copy
Edit
NON_ACADEMIC_KEYWORDS = [...]


 Can I use this without knowing Python?
 Yes. You only need to copy-paste commands in the terminal.

 Can I use this on Windows?
 Yes. Works on Windows, macOS, or Linux with Python installed.

 Will it show university papers?
 No. It filters out academic authors and only keeps company-affiliated ones.

 Can I add more filters or save as Excel?
 Yes! This tool uses pandas, so you can easily enhance it later.

### Example Commands
bash
Copy
Edit
# Search papers from 2023 related to diabetes and save to file
poetry run get-papers-list "diabetes AND 2023[dp]" -f diabetes.csv

# Show filtered papers in terminal only (no file)
poetry run get-papers-list "covid AND 2022[dp]" -d

How It Works (In Simple Language)

You type a search topic
The tool sends this to PubMed's database
It receives a list of research papers
It looks at each paper’s authors and affiliations
It filters out authors from companies, not universities
It saves this filtered list into a CSV file for you



This tool is built using:

BioPython
Pandas
PubMed API
Poetry