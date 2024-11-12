"""
Handlers for output formatting and display.
"""

from typing import Dict, Optional
import pandas as pd

def print_analysis_summary(detailed_results: pd.DataFrame, summary_results: pd.DataFrame) -> None:
    """Print summary of MENTHU analysis results."""
    print("\nResults Overview:")
    print(f"Total PAM sites analyzed: {len(summary_results)}")
    print(f"Total microhomology patterns found: {len(detailed_results)}")
    
    if not summary_results.empty:
        avg_patterns = len(detailed_results) / len(summary_results)
        print(f"Average patterns per PAM site: {avg_patterns:.1f}")
        
        # Only count actual numerical MENTHU scores for the high score count
        high_scores = len(summary_results[
            summary_results['MENTHU score'].apply(lambda x: isinstance(x, (int, float)) and x > 1.5)
        ])
        print(f"Sites with MENTHU score > 1.5: {high_scores}")
        
        print("\nTop scoring PAM sites:")
        top_sites = summary_results.head(3)
        for _, site in top_sites.iterrows():
            print(f"\nPAM site {site['PAM site']} ({site['Strand']})")
            print(f"Spacer: {site['Spacer Sequence']}")
            print(f"PAM: {site['PAM Sequence']}")
            print(f"Total patterns found: {site['Total patterns']}")
            
            # Best pattern information
            print(f"Best pattern:")
            print(f"  mH sequence: {site['Best pattern mH']}")
            print(f"  Length: {site['Best pattern length']}")
            print(f"  Deletion: {site['Best pattern deletion']}")
            print(f"  Score: {site['Best pattern score']:.2f}")
            
            # Second best pattern information (if available)
            if site['Second best mH'] != "NA":
                print(f"Second best pattern:")
                print(f"  mH sequence: {site['Second best mH']}")
                print(f"  Length: {site['Second best length']}")
                print(f"  Deletion: {site['Second best deletion']}")
                print(f"  Score: {float(site['Second best score']):.2f}")
                print(f"MENTHU score: {float(site['MENTHU score']):.2f}")
            else:
                print("Second best pattern: NA")
                print("MENTHU score: NA")
            
            print(f"Frameshift: {site['frameShift']}")
            if 'affected_genes' in site and site['affected_genes']:
                print(f"Affected genes: {site['affected_genes']}")

def format_sequence_visualization(wt_seq: str, del_seq: str, 
                               highlight_indices: Optional[Dict[str, int]] = None) -> str:
    """
    Format sequence visualization with optional highlighting.
    
    Args:
        wt_seq: Wild-type sequence
        del_seq: Sequence with deletion marked
        highlight_indices: Dictionary containing positions to highlight
        
    Returns:
        str: Formatted sequence visualization
    """
    result = []
    result.append("Wild-type sequence:")
    if highlight_indices:
        # Add position markers
        positions = " " * len(wt_seq)
        for label, pos in highlight_indices.items():
            if 0 <= pos < len(positions):
                positions = positions[:pos] + "^" + positions[pos+1:]
        result.append(positions)
    result.append(wt_seq)
    
    result.append("\nDeletion sequence:")
    result.append(del_seq)
    
    return "\n".join(result)

def export_detailed_results(results: pd.DataFrame, output_file: str,
                          include_sequences: bool = True) -> None:
    """
    Export detailed analysis results to CSV file.
    
    Args:
        results: DataFrame containing analysis results
        output_file: Path to output file
        include_sequences: Whether to include full sequence information
    """
    export_df = results.copy()
    
    if not include_sequences:
        # 配列情報を除外
        sequence_columns = ['WT Sequence', 'Del Sequence']
        export_df = export_df.drop(columns=sequence_columns, errors='ignore')
    
    # CSV形式で出力
    export_df.to_csv(output_file, index=False)
    print(f"\nResults exported to: {output_file}")

def create_analysis_report(summary_results: pd.DataFrame, 
                         output_file: str,
                         include_cds: bool = True) -> None:
    """Create comprehensive analysis report."""
    with open(output_file, 'w') as f:
        f.write("MENTHU Analysis Report\n")
        f.write("=====================\n\n")
        
        f.write("Analysis Statistics\n")
        f.write("-----------------\n")
        f.write(f"Total PAM sites analyzed: {len(summary_results)}\n")
        
        # Count high scores (only numerical values)
        high_scores = len(summary_results[
            summary_results['MENTHU score'].apply(lambda x: isinstance(x, (int, float)) and x > 1.5)
        ])
        f.write(f"Sites with MENTHU score > 1.5: {high_scores}\n\n")
        
        f.write("Top Scoring Sites\n")
        f.write("---------------\n")
        top_sites = summary_results.head(5)
        for _, site in top_sites.iterrows():
            f.write(f"\nPAM site {site['PAM site']} ({site['Strand']})\n")
            f.write(f"Spacer: {site['Spacer Sequence']}\n")
            f.write(f"PAM: {site['PAM Sequence']}\n")
            f.write(f"Total patterns: {site['Total patterns']}\n")
            
            f.write("Best pattern:\n")
            f.write(f"  mH sequence: {site['Best pattern mH']}\n")
            f.write(f"  Length: {site['Best pattern length']}\n")
            f.write(f"  Deletion: {site['Best pattern deletion']}\n")
            f.write(f"  Score: {site['Best pattern score']:.2f}\n")
            
            if site['Second best mH'] != "NA":
                f.write("Second best pattern:\n")
                f.write(f"  mH sequence: {site['Second best mH']}\n")
                f.write(f"  Length: {site['Second best length']}\n")
                f.write(f"  Deletion: {site['Second best deletion']}\n")
                f.write(f"  Score: {float(site['Second best score']):.2f}\n")
                f.write(f"MENTHU score: {float(site['MENTHU score']):.2f}\n")
            else:
                f.write("Second best pattern: NA\n")
                f.write("MENTHU score: NA\n")
            
            f.write(f"Frameshift: {site['frameShift']}\n")
            f.write(f"Deletion range: {site['del_start']}-{site['del_end']}\n")
            
            if include_cds and 'affected_genes' in site:
                if site['is_within_CDS']:
                    f.write(f"Affected genes: {site['affected_genes']}\n")
                else:
                    f.write("No CDS overlap\n")
            
            f.write("\n")
        
        if include_cds and 'is_within_CDS' in summary_results.columns:
            f.write("\nCDS Analysis\n")
            f.write("-----------\n")
            cds_overlaps = summary_results['is_within_CDS'].sum()
            total = len(summary_results)
            f.write(f"Sites overlapping CDS: {cds_overlaps} ({(cds_overlaps/total)*100:.1f}%)\n")