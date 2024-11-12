"""
MENTHU (Microhomology-mediated End joining Target HSting Utility)
A tool for analyzing CRISPR target sites with MMEJ consideration.
"""

from .calculators.menthu_calculator import MenthuScoreCalculator, calculate_menthu_score
from .handlers.genbank_handler import retrieve_genbank_file, trim_genbank_record
from .handlers.cds_handler import extract_cds_ranges, add_cds_information
from .handlers.sequence_handler import find_pam_sites

__all__ = [
    'MenthuScoreCalculator',
    'calculate_menthu_score',
    'retrieve_genbank_file',
    'trim_genbank_record',
    'extract_cds_ranges',
    'add_cds_information',
    'find_pam_sites'
]

__version__ = "1.0.0"