"""
Handlers for CDS (Coding Sequence) related operations.
"""

from typing import List, Tuple
import pandas as pd
from Bio.SeqRecord import SeqRecord

def extract_cds_ranges(genbank_record: SeqRecord) -> List[Tuple[int, int, int, str]]:
    """
    Extract CDS ranges from a GenBank record with proper coordinate handling.
    
    Args:
        genbank_record: The GenBank record from which to extract CDS features
        
    Returns:
        List of tuples containing (start, end, strand, gene_name)
    """
    cds_ranges = []
    for feature in genbank_record.features:
        if feature.type == "CDS":
            # Convert to 1-based coordinates
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = int(feature.location.strand)
            
            # Get gene name or product if available
            gene_name = "Unknown"
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]
            elif 'product' in feature.qualifiers:
                gene_name = feature.qualifiers['product'][0]
                
            cds_ranges.append((start, end, strand, gene_name))
            
    print(f"Found {len(cds_ranges)} CDS regions")
    return cds_ranges

def add_cds_information(summary_results: pd.DataFrame, cds_ranges: list) -> pd.DataFrame:
    """
    Add CDS overlap information to the summary results.
    
    Args:
        summary_results: DataFrame containing MENTHU analysis summary
        cds_ranges: List of CDS ranges from extract_cds_ranges()
        
    Returns:
        DataFrame with added CDS information
    """
    if summary_results.empty:
        return summary_results
        
    # CDS情報を追加
    summary_results['is_within_CDS'] = False
    summary_results['affected_genes'] = ''
    
    for idx, row in summary_results.iterrows():
        is_in_cds = False
        affected_genes = []
        
        # 各CDSについて重複をチェック
        for start, end, strand, gene in cds_ranges:
            if row['del_start'] <= end and row['del_end'] >= start:
                is_in_cds = True
                affected_genes.append(gene)
        
        summary_results.at[idx, 'is_within_CDS'] = is_in_cds
        summary_results.at[idx, 'affected_genes'] = '; '.join(affected_genes) if affected_genes else ''
    
    return summary_results

def print_cds_summary(summary_results: pd.DataFrame) -> None:
    """Print summary of CDS overlaps in the results."""
    if 'is_within_CDS' not in summary_results.columns:
        print("\nNo CDS information available")
        return
        
    cds_overlaps = summary_results['is_within_CDS'].sum()
    total = len(summary_results)
    
    print("\nCDS Analysis Summary:")
    print(f"Total sites analyzed: {total}")
    print(f"Sites overlapping CDS: {cds_overlaps} ({(cds_overlaps/total)*100:.1f}%)")
    
    if cds_overlaps > 0:
        print("\nAffected genes:")
        affected = summary_results[summary_results['is_within_CDS']]['affected_genes'].unique()
        for genes in affected:
            if genes:  # Skip empty strings
                print(f"- {genes}")