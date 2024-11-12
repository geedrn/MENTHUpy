#!/usr/bin/env python3
"""
MENTHU (Microhomology-mediated End joining Target HSting Utility)
Command-line interface for CRISPR target site analysis.
"""

import argparse
from pathlib import Path
import logging
import sys

# メインの機能をインポート
from .calculators.menthu_calculator import calculate_menthu_score
from .handlers.genbank_handler import retrieve_genbank_file, trim_genbank_record
from .handlers.cds_handler import extract_cds_ranges, add_cds_information

# ロギング設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('menthu.log')
    ]
)
logger = logging.getLogger(__name__)

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='MENTHU: Microhomology-mediated End joining Target HSting Utility',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--genbank_accession',
        type=str,
        help='GenBank accession number'
    )
    input_group.add_argument(
        '--genbank_file',
        type=str,
        help='Path to local GenBank file'
    )
    
    parser.add_argument(
        '--pam',
        type=str,
        default='NGG',
        help='PAM sequence'
    )
    parser.add_argument(
        '--min_microhomology_length',
        type=int,
        default=3,
        help='Minimum microhomology length'
    )
    parser.add_argument(
        '--start',
        type=int,
        help='Start coordinate (optional)'
    )
    parser.add_argument(
        '--end',
        type=int,
        help='End coordinate (optional)'
    )
    parser.add_argument(
        '--output_prefix',
        '-o',
        type=str,
        default='MENTHU_output',
        help='Prefix for output files'
    )
    parser.add_argument(
        '--search_range',
        type=int,
        default=50,
        help='Search range around DSB site'
    )
    parser.add_argument(
        '--min_mh_distance',
        type=int,
        default=0,
        help='Minimum distance between microhomology sequences'
    )
    parser.add_argument(
        '--max_mh_distance',
        type=int,
        help='Maximum distance between microhomology sequences (optional)'
    )
    parser.add_argument(
        '--min_menthu_score',
        type=float,
        help='Minimum MENTHU score (optional)' 
    )
    
    return parser.parse_args()

def run_menthu_analysis(args: argparse.Namespace) -> bool:
    """
    Run the MENTHU analysis pipeline.
    """
    try:
        # GenBankデータの取得と処理
        genbank_record = retrieve_genbank_file(
            accession_id=args.genbank_accession,
            file_path=args.genbank_file
        )
        
        # CDS情報の取得
        cds_ranges = extract_cds_ranges(genbank_record)
        
        trimmed_record = trim_genbank_record(
            genbank_record,
            start_coord=args.start,
            end_coord=args.end
        )
        
        # MENTHU解析の実行
        results = calculate_menthu_score(
            trimmed_record=trimmed_record,
            cds_ranges=cds_ranges,  # CDS情報を渡す
            pam=args.pam,
            min_microhomology_length=args.min_microhomology_length,
            start_coord=args.start or 1,
            output_prefix=args.output_prefix,
            search_range=args.search_range,
            min_mh_distance=args.min_mh_distance,
            max_mh_distance=args.max_mh_distance,
            weight=20.0
        )
        
        if results.empty:
            logger.warning("No results found")
            return False
        
        return True
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}", exc_info=True)
        return False

def main() -> int:
    """Main entry point for the CLI."""
    args = parse_arguments()
    success = run_menthu_analysis(args)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())