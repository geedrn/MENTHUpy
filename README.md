# MENTHU Python version

MENTHUpy is a Python tool for analyzing CRISPR target sites with consideration of Microhomology-mediated End joining (MMEJ) patterns. It helps researchers identify and evaluate potential CRISPR target sites based on their microhomology patterns and predicted repair outcomes.
MENTHU (https://github.com/FriedbergLab/MENTHU) is originally implemented in R and this software mainly relies on the idea from MENTHU (but not identical in some points). 

## Features

- Analyze CRISPR target sites for MMEJ patterns
- Calculate MENTHU scores for target site evaluation
- Support for various PAM sequences
- Integration with GenBank for sequence retrieval
- CDS region analysis
- Detailed output reports

## Installation

### Prerequisites

- Python 3.7 (Use conda and set the virtual environment like conda create -n menthu -c conda-forge python=3.7)
- pip (Python package installer)

### Method to install the software

1. Clone the repository:
```bash
git clone https://github.com/yourusername/menthu.git
cd menthu
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Install the package:
```bash
pip install .
```

## Usage

### Basic Command

```bash
menthu --genbank_file input.gb
```

### Command Line Options

```
Required arguments (one of):
  --genbank_accession ACCESSION  GenBank accession number
  --genbank_file FILE           Path to local GenBank file

Optional arguments:
  --pam PAM                     PAM sequence (default: NGG)
  --min_microhomology_length N  Minimum microhomology length (default: 3)
  --start START                 Start coordinate (optional)
  --end END                     End coordinate (optional)
```

### Example Usage

1. Analyze a local GenBank file:
```bash
menthu --genbank_file example.gb --pam NGG
```

2. Analyze a sequence from NCBI:
```bash
menthu --genbank_accession NC_000001.11 --start 1000 --end 2000
```

