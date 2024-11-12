"""
Handlers for sequence-related operations.
"""

from typing import List, Tuple, Dict
from Bio.Seq import Seq
import re

def find_pam_sites(sequence: str, start_coord: int = 1, pam: str = "NGG") -> List[Tuple[int, str, str]]:
    """
    Identify PAM sites in the sequence and return positions with strand info.
    
    Args:
        sequence: The DNA sequence to analyze
        start_coord: Starting coordinate in reference sequence
        pam: PAM sequence pattern
        
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
    
    pam_pattern = "".join([iupac_dict.get(base, base) for base in pam])
    regex = re.compile(pam_pattern)
    
    # Forward strand
    pam_sites = []
    for match in regex.finditer(sequence):
        absolute_position = start_coord + match.start()
        spacer_start = max(0, match.start()-20)
        spacer = sequence[spacer_start:match.start()]
        pam_seq = sequence[match.start():match.start()+len(pam)]
        
        if len(spacer) == 20:
            pam_sites.append((absolute_position, spacer, pam_seq))
    
    # Reverse strand
    rev_comp = str(Seq(sequence).reverse_complement())
    for match in regex.finditer(rev_comp):
        position = len(sequence) - match.start() - len(pam)
        absolute_position = start_coord + position
        spacer_start = max(0, match.start()-20)
        spacer = rev_comp[spacer_start:match.start()]
        pam_seq = rev_comp[match.start():match.start()+len(pam)]
        
        if len(spacer) == 20:
            pam_sites.append((absolute_position, spacer, pam_seq))
    
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