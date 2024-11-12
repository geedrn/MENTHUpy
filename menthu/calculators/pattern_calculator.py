"""
Pattern score calculation and related utilities.
"""

from typing import Dict, List, Tuple, Optional, Set
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math

class PatternScoreCalculator:
    def __init__(self,
                 search_range: int = 50,
                 min_microhomology_length: int = 3,
                 weight: float = 20.0,
                 min_mh_distance: int = 0,
                 max_mh_distance: Optional[int] = None):
        self.search_range = search_range
        self.min_microhomology_length = min_microhomology_length
        self.weight = weight
        self.min_mh_distance = min_mh_distance
        self.max_mh_distance = max_mh_distance

    def find_common_substrings(self, upstream: str, downstream: str, 
                            min_length: int = 3) -> List[Dict]:
        """
        Find valid microhomology sequences with consistent boundary detection.
        
        Args:
            upstream: Upstream sequence (50bp)
            downstream: Downstream sequence (50bp)
            min_length: Minimum length for microhomology
            
        Returns:
            List of dictionaries containing microhomology information
        """
        if len(upstream) < self.search_range or len(downstream) < self.search_range:
            return []
        
        microhomologies = []
        
        for length in range(min_length, min(20, len(upstream), len(downstream)) + 1):
            for i in range(len(upstream) - length + 1):
                substring = upstream[i:i+length]
                
                for j in range(len(downstream) - length):
                    if downstream[j:j+length] == substring:
                        # 1-based coordinates relative to the DSB site
                        up_start = i + 1
                        up_end = up_start + length - 1
                        down_start = j + 1
                        down_end = down_start + length - 1
                        
                        # 削除長の計算を修正
                        # DSBサイトを中心に前後50bpを考慮
                        del_length = self.search_range - up_end + down_end

                        if down_end < len(downstream):
                            microhomologies.append({
                                'sequence': substring,
                                'length': length,
                                'upstream_start': up_start,
                                'upstream_end': up_end,
                                'downstream_start': down_start,
                                'downstream_end': down_end,
                                'del_length': del_length,
                                'gc_content': substring.count('G') + substring.count('C')
                            })
        
        return sorted(microhomologies, 
                    key=lambda x: (x['del_length'], -x['length'], -x['gc_content']))

    def process_pam_site(self, pos: int, strand: str, spacer_seq: str, pam_seq: str,
                        seq_record: SeqRecord, start_coord: int) -> List[Dict]:
        """
        Process PAM site and identify microhomology patterns.
        """
        if strand == '-':
            dsb_site = pos + 6
        else:
            dsb_site = pos - 3
        
        dsb_idx = dsb_site - start_coord
        seq_len = len(seq_record.seq)
        
        # DSBサイトを中心に前後50bpを取得
        upstream_start = max(0, dsb_idx - self.search_range)
        upstream_end = dsb_idx
        downstream_start = dsb_idx
        downstream_end = min(seq_len, dsb_idx + self.search_range)
        
        upstream = seq_record.seq[upstream_start:upstream_end]
        downstream = seq_record.seq[downstream_start:downstream_end]
        
        # シーケンス長の検証
        if not self.validate_sequence_lengths(
            upstream=upstream, 
            downstream=downstream, 
            pos=pos, 
            dsb_site=dsb_site, 
            strand=strand
        ):
            return []
        
        patterns = []
        microhomologies = self.find_common_substrings(
            str(upstream),
            str(downstream),
            min_length=self.min_microhomology_length
        )
        
        for mh in microhomologies:
            pattern_info = self.process_microhomology(mh, upstream, downstream)
            if pattern_info is None:
                continue
                
            patterns.append({
                'pos': pos,
                'strand': strand,
                'dsb_site': dsb_site,
                'spacer_seq': spacer_seq,
                'pam_seq': pam_seq,
                'microhomology': mh,
                'pattern_info': pattern_info,
                'search_range': self.search_range,
            })
            
        return patterns

    def calculate_pattern_score(self, del_length: int, gc_count: int, mh_length: int) -> float:
        """Calculate pattern score based on deletion length, GC content, and microhomology length."""
        if del_length < 0:
            return 0.0
        return 100 * (1 / math.exp(del_length / self.weight)) * (gc_count + mh_length)

    def get_sequence_lengths(self, dsb_idx: int, seq_length: int) -> Tuple[int, int, int, int]:
        """Calculate sequence ranges around DSB site."""
        upstream_start = max(0, dsb_idx - self.search_range)
        upstream_end = dsb_idx
        downstream_start = dsb_idx
        downstream_end = min(seq_length, dsb_idx + self.search_range)
        
        return upstream_start, upstream_end, downstream_start, downstream_end

    def validate_sequence_lengths(self, upstream: Seq, downstream: Seq, 
                                pos: int, dsb_site: int, strand: str) -> bool:
        """
        Validate sequence lengths and print detailed information if invalid.
        
        Args:
            upstream: Upstream sequence
            downstream: Downstream sequence
            pos: PAM site position
            dsb_site: DSB site position
            strand: Strand orientation
            
        Returns:
            bool: True if sequence lengths are valid
        """
        if len(upstream) < self.search_range or len(downstream) < self.search_range:
            print(f"\nSkipping PAM site at position {pos} (insufficient sequence length)")
            print(f"  - Upstream length: {len(upstream)} (required: {self.search_range})")
            print(f"  - Downstream length: {len(downstream)} (required: {self.search_range})")
            print(f"  - DSB site: {dsb_site}, Strand: {strand}")
            return False
        return True

    def process_microhomology(self, mh: Dict, upstream: Seq, downstream: Seq) -> Optional[Dict]:
        """Process and validate microhomology pattern."""
        if mh['del_length'] < 0:
            return None

        # 距離の計算と検証
        relative_upstream_end = mh['upstream_end']
        relative_downstream_start = mh['downstream_start']
        absolute_upstream_end = relative_upstream_end
        absolute_downstream_start = relative_downstream_start + len(downstream)
        mh_distance = absolute_downstream_start - absolute_upstream_end - 1
        if ((self.min_mh_distance is not None and mh_distance < self.min_mh_distance) or
            (self.max_mh_distance is not None and mh_distance > self.max_mh_distance)):
            return None

        # パターンスコアの計算
        pattern_score = self.calculate_pattern_score(
            mh['del_length'],
            mh['gc_content'],
            mh['length']
        )

        # シーケンスの可視化
        sequence_vis = self.create_sequence_visualization(
            upstream=upstream,
            downstream=downstream,
            mh_data=mh
        )

        return {
            'pattern_score': pattern_score,
            'mh_distance': mh_distance,
            'visualization': sequence_vis,
            'frameshift': mh['del_length'] % 3 != 0
        }

    def create_sequence_visualization(self, upstream: Seq, downstream: Seq,
                                mh_data: Dict) -> Dict[str, str]:
        """
        Create sequence visualization with microhomology markers.
        """
        upstream_str = str(upstream)
        downstream_str = str(downstream)
        wt_seq = upstream_str + downstream_str
        dsb_pos = len(upstream_str)

        # 1-based to 0-based conversion
        mh_start_up = mh_data['upstream_start'] - 1
        mh_end_up = mh_data['upstream_end']
        mh_start_down = mh_data['downstream_start'] - 1
        mh_end_down = mh_data['downstream_end']

        # マーカー付きシーケンスの作成
        wt_with_markers = (
            # Upstream region with microhomology
            wt_seq[:mh_start_up] +
            "[" + wt_seq[mh_start_up:mh_end_up] + "]" +
            wt_seq[mh_end_up:dsb_pos] +
            # DSB marker
            '*' +
            # Downstream region with microhomology
            downstream_str[:mh_start_down] +
            "[" + downstream_str[mh_start_down:mh_end_down] + "]" +
            downstream_str[mh_end_down:]
        )

        # 削除配列（マイクロホモロジーは1回だけ含める）
        del_seq = (
            wt_seq[:mh_end_up] +
            downstream_str[mh_end_down:]
        )

        # 削除の可視化
        del_seq_viz = (
            wt_seq[:mh_end_up] +
            "-" * mh_data['del_length'] +
            downstream_str[mh_end_down:]
        )

        return {
            'wt_seq': wt_with_markers,
            'del_seq': del_seq,
            'del_seq_viz': del_seq_viz
        }