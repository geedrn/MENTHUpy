"""
Handlers for sequence-related operations.
"""

from Bio.Seq import Seq
import re
from typing import List, Tuple

def find_pam_sites(sequence: str, start_coord: int = 1, pam: str = "NGG", spacer_length: int = 20) -> Tuple[List[Tuple[int, str, str]], List[Tuple[int, str, str]]]:
    """
    Identify PAM sites in the sequence and return positions with strand info.
    
    Args:
        sequence: The DNA sequence to analyze
        start_coord: Starting coordinate in reference sequence
        pam: PAM sequence pattern (e.g., NGG)
        spacer_length: Length of the spacer sequence (default 20)
    
    Returns:
        forward_sites: List of (position, spacer, PAM) for forward strand
        reverse_sites: List of (position, spacer, PAM) for reverse strand
    """
    # IUPAC codes for the PAM sequence
    iupac_dict = {
        "N": "[ATGC]",
        "R": "[AG]",
        "Y": "[CT]",
        "S": "[GC]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]"
    }
    
    # Convert PAM to regex pattern
    pam_pattern = "".join([iupac_dict.get(base, base) for base in pam])
    pam_regex = re.compile(f'(?=({pam_pattern}))')  # Use lookahead to match PAM without consuming it
    
    # Forward strand
    forward_sites = []
    for match in pam_regex.finditer(sequence):
        pam_start = match.start()
        spacer_start = pam_start - spacer_length
        if spacer_start >= 0:
            spacer = sequence[spacer_start:pam_start]
            pam_seq = sequence[pam_start:pam_start + len(pam)]
            position = start_coord + pam_start
            forward_sites.append((position, spacer, pam_seq))
    
    # Reverse strand
    reverse_sequence = str(Seq(sequence).reverse_complement())
    reverse_sites = []
    for match in pam_regex.finditer(reverse_sequence):
        pam_start = match.start()
        spacer_start = pam_start - spacer_length
        if spacer_start >= 0:
            spacer = reverse_sequence[spacer_start:pam_start]
            pam_seq = reverse_sequence[pam_start:pam_start + len(pam)]
            # Adjust position to original sequence coordinates
            # Adjust to a first base of the PAM sequence
            position = start_coord + len(sequence) - pam_start
            reverse_sites.append((position, spacer, pam_seq))
    
    return forward_sites, reverse_sites

def extract_spacer_and_pam(sequence: Seq, position: int,
                          pam_length: int, strand: str) -> Tuple[str, str]:
    """
    Extract spacer and PAM sequences based on strand information.
    
    Args:
        sequence: The DNA sequence
        position: Position in the sequence
        pam_length: Length of the PAM sequence
        strand: Strand orientation ('+' or '-')
        
    Returns:
        Tuple containing (spacer_sequence, pam_sequence)
    """
    if strand == '-':
        seq = sequence.reverse_complement()
        position = len(sequence) - position - 1
    else:
        seq = sequence

    spacer_sequence = seq[position - 20 : position]
    pam_sequence = seq[position : position + pam_length]

    return str(spacer_sequence), str(pam_sequence)

def validate_pam_sequence(pam_sequence: str, strand: str) -> bool:
    """
    Check if PAM sequence is valid for the given strand.
    
    Args:
        pam_sequence: PAM sequence to validate
        strand: Strand orientation ('+' or '-')
        
    Returns:
        bool: True if PAM sequence is valid
    """
    if strand == '-' and not pam_sequence.endswith('GG'):
        return False
    return True

class CasConfig:
    def __init__(self, pam: str, pam_length: int, cut_offset: int):
        self.pam = pam
        self.pam_length = pam_length
        self.cut_offset = cut_offset

# 各Cas酵素の設定
SpCas9 = CasConfig(pam="NGG", pam_length=3, cut_offset=3)
AsCas12a = CasConfig(pam="TTTN", pam_length=4, cut_offset=18)