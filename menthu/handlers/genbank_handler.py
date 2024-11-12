"""
Handlers for GenBank file operations.
"""

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import os
from typing import Optional

def retrieve_genbank_file(accession_id: Optional[str] = None,
                         file_path: Optional[str] = None,
                         email: str = "your_email@example.com") -> SeqRecord:
    """
    Retrieve a GenBank file either from a local file or NCBI database.
    
    Args:
        accession_id: NCBI accession ID for fetching from the database
        file_path: Local file path for the GenBank file
        email: Email address for NCBI queries
        
    Returns:
        The GenBank file content as a Biopython SeqRecord
        
    Raises:
        ValueError: If neither accession_id nor file_path is provided or if file not found
    """
    try:
        if file_path:
            # ファイルパスの正規化（空白の除去など）
            file_path = file_path.strip()
            
            if not os.path.exists(file_path):
                raise ValueError(f"GenBank file not found: {file_path}")
                
            print(f"Loading local GenBank file: {file_path}")
            
            try:
                with open(file_path, "r") as handle:
                    genbank_record = SeqIO.read(handle, "genbank")
                    
                    if not genbank_record or not genbank_record.seq:
                        raise ValueError("Invalid GenBank record: Empty sequence")
                    
                    return genbank_record
            except Exception as e:
                raise ValueError(f"Error reading GenBank file: {str(e)}")
                
        elif accession_id:
            print(f"Fetching from NCBI. Accession: {accession_id}")
            Entrez.email = email
            
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
            genbank_record = SeqIO.read(handle, "genbank")
            handle.close()
            
            if not genbank_record or not genbank_record.seq:
                raise ValueError("Invalid GenBank record from NCBI")
            
            return genbank_record
            
        else:
            raise ValueError("Either file_path or accession_id must be provided")
            
    except Exception as e:
        print(f"Error retrieving GenBank file: {str(e)}")
        import traceback
        traceback.print_exc()
        raise

def trim_genbank_record(genbank_record: SeqRecord,
                       start_coord: Optional[int] = None,
                       end_coord: Optional[int] = None) -> SeqRecord:
    """
    Trim a GenBank record to a specified range based on chromosome coordinates.
    
    Args:
        genbank_record: The GenBank file content as a Biopython SeqRecord
        start_coord: The starting coordinate for trimming
        end_coord: The ending coordinate for trimming
        
    Returns:
        A trimmed SeqRecord containing only the specified range
    """
    try:
        if not genbank_record or not hasattr(genbank_record, 'seq'):
            raise ValueError("Invalid input: Record has no sequence data")
            
        # Handle default coordinates
        if start_coord is None:
            start_coord = 1
        if end_coord is None:
            end_coord = len(genbank_record.seq)
            
        # Convert to 0-based indexing
        start_idx = start_coord - 1 if start_coord > 0 else 0
        end_idx = end_coord
        
        # Validate coordinates
        if start_idx < 0:
            raise ValueError("Start coordinate must be >= 1")
        if end_idx > len(genbank_record.seq):
            raise ValueError(f"End coordinate exceeds sequence length")
        if start_idx >= end_idx:
            raise ValueError("Invalid coordinate range")
            
        # Trim sequence
        trimmed_seq = genbank_record.seq[start_idx:end_idx]
        
        # Create new record
        trimmed_record = SeqRecord(
            trimmed_seq,
            id=genbank_record.id,
            name=genbank_record.name,
            description=f"Trimmed {genbank_record.description} from {start_coord} to {end_coord}"
        )
        
        # Copy annotations
        if hasattr(genbank_record, 'annotations'):
            trimmed_record.annotations = genbank_record.annotations.copy()
        
        if 'molecule_type' not in trimmed_record.annotations:
            trimmed_record.annotations['molecule_type'] = 'DNA'
        
        # Process features
        trimmed_record.features = []
        for feature in genbank_record.features:
            if feature.location.start < end_idx and feature.location.end > start_idx:
                new_start = max(0, feature.location.start - start_idx)
                new_end = min(len(trimmed_seq), feature.location.end - start_idx)
                
                if new_start < new_end:
                    from Bio.SeqFeature import SeqFeature, FeatureLocation
                    new_feature = SeqFeature(
                        FeatureLocation(new_start, new_end, feature.location.strand),
                        type=feature.type,
                        qualifiers=feature.qualifiers.copy()
                    )
                    trimmed_record.features.append(new_feature)
        
        return trimmed_record
        
    except Exception as e:
        print(f"Error in trim_genbank_record: {str(e)}")
        import traceback
        traceback.print_exc()
        raise