"""
Main calculator class for MENTHU analysis.
"""

from typing import Dict, List, Tuple, Optional
import pandas as pd
from Bio.SeqRecord import SeqRecord
from menthu.calculators.pattern_calculator import PatternScoreCalculator
from menthu.handlers.sequence_handler import find_pam_sites
from menthu.calculators.indelphi_calculator import get_indelphi_predictions

class MenthuScoreCalculator:
    def __init__(self,
                 trimmed_record: SeqRecord, 
                 cds_ranges: List[Tuple],
                 pam: str = "NGG",
                 min_microhomology_length: int = 3,
                 start_coord: int = 1,
                 weight: float = 20.0,
                 search_range: int = 50,
                 min_mh_distance: int = 0,
                 max_mh_distance: Optional[int] = None,
                 min_menthu_score: float = 1.5):
        """
        Initialize MENTHU Score Calculator.
        
        Args:
            trimmed_record: The sequence record to analyze
            cds_ranges: List of tuples containing CDS information (start, end, strand, gene_name)
            pam: PAM sequence pattern
            min_microhomology_length: Minimum length for microhomology
            start_coord: Starting coordinate
            weight: Weight factor for pattern score calculation
            search_range: Search range around DSB site
            min_mh_distance: Minimum distance between microhomology sequences
            max_mh_distance: Maximum distance between microhomology sequences
            min_menthu_score: Minimum MENTHU score to include in summary
        """
        self.trimmed_record = trimmed_record
        self.cds_ranges = cds_ranges
        self.pam = pam
        self.start_coord = start_coord
        self.results = []
        self.min_menthu_score = min_menthu_score
        
        self.pattern_calculator = PatternScoreCalculator(
            search_range=search_range,
            min_microhomology_length=min_microhomology_length,
            weight=weight,
            min_mh_distance=min_mh_distance,
            max_mh_distance=max_mh_distance
        )

    def _process_forward_strand(self) -> None:
        """Process PAM sites on forward strand."""
        sequence_str = str(self.trimmed_record.seq)
        pam_sites = find_pam_sites(sequence_str, start_coord=self.start_coord, pam=self.pam)
        
        for absolute_position, spacer, pam_seq in pam_sites:
            if 'N' not in spacer:
                self._process_pam_site(absolute_position, '+', spacer, pam_seq)

    def _process_reverse_strand(self) -> None:
        """Process PAM sites on reverse strand."""
        rev_comp = str(self.trimmed_record.seq.reverse_complement())
        pam_sites = find_pam_sites(rev_comp, start_coord=self.start_coord, pam=self.pam)
        
        for absolute_position, spacer, pam_seq in pam_sites:
            position = len(self.trimmed_record.seq) - (absolute_position - self.start_coord) - len(self.pam)
            absolute_position = self.start_coord + position
            
            if 'N' not in spacer:
                self._process_pam_site(absolute_position, '-', spacer, pam_seq)

    def _process_pam_site(self, pos: int, strand: str, spacer_seq: str, pam_seq: str) -> None:
        """Process individual PAM site and store results."""
        patterns = self.pattern_calculator.process_pam_site(
            pos, strand, spacer_seq, pam_seq,
            self.trimmed_record, self.start_coord
        )
        
        for pattern in patterns:
            self._store_result(
                pattern['pos'],
                pattern['strand'],
                pattern['dsb_site'],
                pattern['spacer_seq'],
                pattern['pam_seq'],
                pattern['microhomology'],
                pattern['pattern_info']
            )

    def _store_result(self, pos: int, strand: str, dsb_site: int,
                 spacer_seq: str, pam_seq: str,
                 mh: Dict, pattern_info: Dict) -> None:
        """Store processed results without CDS check."""
        self.results.append({
            "PAM site": pos,
            "Strand": strand,
            "DSB site": dsb_site,
            "Spacer Sequence": spacer_seq,
            "PAM Sequence": pam_seq,
            "mH sequence": mh['sequence'],
            "mH length": mh['length'],
            "GC content": mh['gc_content'],
            "delLength": mh['del_length'],
            "patternScore": pattern_info['pattern_score'],
            "frameShift": "Yes" if pattern_info['frameshift'] else "No",
            "WT Sequence": pattern_info['visualization']['wt_seq'],
            "Del Sequence": pattern_info['visualization']['del_seq_viz'],
            "Real Del Sequence": pattern_info['visualization']['del_seq'],
            "mh_start_up": mh['upstream_start'],
            "mh_end_up": mh['upstream_end'],
            "mh_start_down": mh['downstream_start'],
            "mh_end_down": mh['downstream_end'],
            "mh_distance": pattern_info['mh_distance'],
            "del_start": dsb_site - (self.pattern_calculator.search_range - mh['upstream_end']),
            "del_end": dsb_site - (self.pattern_calculator.search_range - mh['upstream_end']) + mh['del_length'] - 1
        })

    def _calculate_scores(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate MENTHU scores and return both detailed and summary results."""
        try:
            self._process_forward_strand()
            self._process_reverse_strand()
            
            if not self.results:
                return pd.DataFrame(), pd.DataFrame()

            detailed_df = pd.DataFrame(self.results)
            summary_df = self._create_summary(detailed_df)
            
            return detailed_df, summary_df
            
        except Exception as e:
            print(f"Error in calculate_scores: {str(e)}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame(), pd.DataFrame()
    
    def _check_cds_overlap(self, del_start: int, del_end: int) -> Tuple[bool, List[str]]:
        """
        Check if deletion overlaps with any CDS regions.
        
        Args:
            del_start: Deletion start position
            del_end: Deletion end position
            
        Returns:
            Tuple of (is_in_cds, affected_genes)
        """
        is_in_cds = False
        affected_genes = []
        
        for start, end, cds_strand, gene in self.cds_ranges:
            if del_start <= end and del_end >= start:
                is_in_cds = True
                affected_genes.append(gene)
        
        return is_in_cds, affected_genes

    def _create_summary(self, detailed_df: pd.DataFrame) -> pd.DataFrame:
        """Create summary DataFrame with MENTHU scores and inDelphi predictions."""
        summary_data = []
        search_range = self.pattern_calculator.search_range

        # まずDSBサイトとストランドでグループ化
        for (dsb_site, strand), group in detailed_df.groupby(["DSB site", "Strand"]):
            # 3bp以上のmicrohomologyでフィルタ
            valid_patterns = group[group['mH length'] >= 3]
            if valid_patterns.empty:
                continue

            # Real Del Sequenceでさらにグループ化
            del_seq_groups = valid_patterns.groupby("Real Del Sequence").agg({
                "patternScore": "max",
                "mH sequence": "first",
                "mH length": "first",
                "GC content": "first",
                "delLength": "first",
                "WT Sequence": "first",
                "Del Sequence": "first",
                "del_start": "first",
                "del_end": "first",
                "frameShift": "first",
                "PAM site": "first",
                "Spacer Sequence": "first",
                "PAM Sequence": "first"
            }).reset_index()

            del_seq_groups = del_seq_groups.sort_values("patternScore", ascending=False)
            
            if len(del_seq_groups) < 2:
                continue

            best_pattern = del_seq_groups.iloc[0]
            second_pattern = del_seq_groups.iloc[1]
            
            menthu_score = best_pattern["patternScore"] / second_pattern["patternScore"]
            
            if menthu_score < self.min_menthu_score:
                continue

            # CDS判定（ベストパターンのみ）
            is_in_cds, affected_genes = self._check_cds_overlap(
                best_pattern["del_start"],
                best_pattern["del_end"]
            )

            # inDelphi予測
            dsb_idx = dsb_site - self.start_coord
            upstream = str(self.trimmed_record.seq[max(0, dsb_idx - search_range):dsb_idx])
            downstream = str(self.trimmed_record.seq[dsb_idx:min(dsb_idx + search_range, len(self.trimmed_record.seq))])
            
            if len(upstream) < search_range:
                upstream = "N" * (search_range - len(upstream)) + upstream
            if len(downstream) < search_range:
                downstream = downstream + "N" * (search_range - len(downstream))

            indelphi_results = get_indelphi_predictions(upstream, downstream)
            indelphi_data = self._process_indelphi_results(indelphi_results, best_pattern["Real Del Sequence"])

            # サマリーデータの作成
            summary_data.append({
                "PAM site": best_pattern["PAM site"],
                "Strand": strand,
                "DSB site": dsb_site,
                "Spacer Sequence": best_pattern["Spacer Sequence"],
                "PAM Sequence": best_pattern["PAM Sequence"],
                "Total patterns": len(valid_patterns),
                "Unique deletion patterns": len(del_seq_groups),
                "Best mH": best_pattern["mH sequence"],
                "Best mH length": best_pattern["mH length"],
                "Best mH GC": best_pattern["GC content"],
                "Best deletion length": best_pattern["delLength"],
                "Best pattern score": best_pattern["patternScore"],
                "Second best mH": second_pattern["mH sequence"],
                "Second best mH length": second_pattern["mH length"],
                "Second best deletion length": second_pattern["delLength"],
                "Second best score": second_pattern["patternScore"],
                "MENTHU score": menthu_score,
                "Query sequence": upstream + downstream,
                "WT sequence with markers": best_pattern["WT Sequence"],
                "Del Sequence": best_pattern["Del Sequence"],
                "Second best WT sequence with markers": second_pattern["WT Sequence"],
                "Second best Del Sequence": second_pattern["Del Sequence"],
                "del_start": best_pattern["del_start"],
                "del_end": best_pattern["del_end"],
                "frameShift": best_pattern["frameShift"],
                "is_within_CDS": is_in_cds,
                "affected_genes": '; '.join(affected_genes) if affected_genes else '',
                **indelphi_data
            })

        return pd.DataFrame(summary_data)

    def _process_indelphi_results(self, indelphi_results: Optional[Dict], 
                                best_del_sequence: str) -> Dict:
        """Process inDelphi prediction results."""
        if not indelphi_results:
            return {
                "inDelphi Top1 pattern": "NA", "inDelphi Top1 frequency": "NA",
                "inDelphi Top1 category": "NA", "inDelphi Top1 length": "NA",
                "inDelphi Top2 pattern": "NA", "inDelphi Top2 frequency": "NA",
                "inDelphi Top2 category": "NA", "inDelphi Top2 length": "NA",
                "inDelphi Top3 pattern": "NA", "inDelphi Top3 frequency": "NA",
                "inDelphi Top3 category": "NA", "inDelphi Top3 length": "NA",
                "inDelphi Top1 match": "NA"
            }
        
        return {
            "inDelphi Top1 pattern": indelphi_results["Top3 genotypes"][0],
            "inDelphi Top1 frequency": indelphi_results["Top3 frequencies"][0],
            "inDelphi Top1 category": indelphi_results["Top3 categories"][0],
            "inDelphi Top1 length": indelphi_results["Top3 lengths"][0],
            "inDelphi Top2 pattern": indelphi_results["Top3 genotypes"][1],
            "inDelphi Top2 frequency": indelphi_results["Top3 frequencies"][1],
            "inDelphi Top2 category": indelphi_results["Top3 categories"][1],
            "inDelphi Top2 length": indelphi_results["Top3 lengths"][1],
            "inDelphi Top3 pattern": indelphi_results["Top3 genotypes"][2],
            "inDelphi Top3 frequency": indelphi_results["Top3 frequencies"][2],
            "inDelphi Top3 category": indelphi_results["Top3 categories"][2],
            "inDelphi Top3 length": indelphi_results["Top3 lengths"][2],
            "inDelphi Top1 match": "TRUE" if best_del_sequence == indelphi_results["Top3 genotypes"][0] else "FALSE"
        }

def calculate_menthu_score(trimmed_record: SeqRecord, 
                         cds_ranges: List[Tuple],  # CDS情報を追加
                         pam: str = "NGG", 
                         min_microhomology_length: int = 3, 
                         start_coord: int = 1, 
                         output_prefix: str = "MENTHU_output",
                         search_range: int = 50,
                         min_mh_distance: int = 0,
                         max_mh_distance: Optional[int] = None,
                         weight: float = 20.0,
                         min_menthu_score: float = 1.5) -> pd.DataFrame:
    """
    Calculate MENTHU scores and export results.
    
    Args:
        trimmed_record: Biopython SeqRecord object
        cds_ranges: List of tuples containing CDS information (start, end, strand, gene_name)
        pam: PAM sequence pattern (default: NGG)
        min_microhomology_length: Minimum length for microhomology sequences
        start_coord: Starting coordinate in the reference sequence
        output_prefix: Prefix for output files
        search_range: Search range around DSB site (default: 50)
        min_mh_distance: Minimum distance between microhomology sequences (default: 0)
        max_mh_distance: Maximum distance between microhomology sequences (optional)
        weight: Weight factor for pattern score calculation (default: 20.0)
        min_menthu_score: Minimum MENTHU score to include in summary (default: 1.5)
        
    Returns:
        pd.DataFrame: Analysis results filtered by minimum MENTHU score
    """
    try:
        if not isinstance(trimmed_record, SeqRecord):
            raise ValueError("trimmed_record must be a Bio.SeqRecord object")
        
        calculator = MenthuScoreCalculator(
            trimmed_record=trimmed_record,
            cds_ranges=cds_ranges, 
            pam=pam,
            min_microhomology_length=min_microhomology_length,
            start_coord=start_coord,
            weight=weight,
            search_range=search_range,
            min_mh_distance=min_mh_distance,
            max_mh_distance=max_mh_distance,
            min_menthu_score=min_menthu_score
        )
        
        detailed_results, summary_results = calculator._calculate_scores()
        
        if not detailed_results.empty:
            detailed_file = f"{output_prefix}_pattern_score.csv"
            summary_file = f"{output_prefix}_menthu_score.csv"
            
            detailed_results.to_csv(detailed_file, index=False)
            summary_results.to_csv(summary_file, index=False)
            
            print(f"\nDetailed results exported to: {detailed_file}")
            print(f"Summary results exported to: {summary_file}")
            print(f"Number of sites with MENTHU score >= {min_menthu_score}: {len(summary_results)}")
            print(f"Number of sites within CDS: {summary_results['is_within_CDS'].sum()}")
            
            return summary_results
        else:
            print("No microhomology patterns found")
            return pd.DataFrame()
            
    except Exception as e:
        print(f"Error in calculate_menthu_score: {str(e)}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()