"""
Handlers for sequence-related operations.
"""

from typing import List, Tuple, Dict
from Bio.Seq import Seq
import re

def find_pam_sites(sequence: str, start_coord: int = 1, pam: str = "NGG", min_spacer_length: int = 20) -> List[Tuple[int, str, str]]:
    """
    Identify PAM sites in the sequence and return positions with strand info.
    
    Args:
        sequence: The DNA sequence to analyze
        start_coord: Starting coordinate in reference sequence
        pam: PAM sequence pattern
        min_spacer_length: Minimum required spacer length (default: 20)
        
    Returns:
        List of tuples containing (position, spacer_sequence, pam_sequence)
    """
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
    
    # Generate regex pattern with positive lookahead for overlapping matches
    pam_pattern = "".join([iupac_dict.get(base, base) for base in pam])
    regex = re.compile(f'(?=({pam_pattern}))')
    
    def find_strand_sites(seq: str, is_reverse: bool = False) -> List[Tuple[int, str, str]]:
        sites = []
        for match in regex.finditer(seq):
            position = match.start()
            if is_reverse:
                position = len(sequence) - position - len(pam)
            
            absolute_position = start_coord + position
            
            # Get spacer sequence
            spacer_start = max(0, match.start()-min_spacer_length)
            spacer = seq[spacer_start:match.start()]
            pam_seq = seq[match.start():match.start()+len(pam)]
            
            # Include sites even with shorter spacers, but flag them
            if len(spacer) < min_spacer_length:
                spacer = "!" + spacer  # Flag short spacers
            
            sites.append((absolute_position, spacer, pam_seq))
        return sites
    
    # Forward strand sites
    pam_sites = find_strand_sites(sequence)
    
    # Reverse strand sites
    rev_comp = str(Seq(sequence).reverse_complement())
    pam_sites.extend(find_strand_sites(rev_comp, is_reverse=True))
    
    # Sort by position
    return sorted(pam_sites)

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