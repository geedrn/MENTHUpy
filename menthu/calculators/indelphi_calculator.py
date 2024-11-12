# menthu/calculators/indelphi_calculator.py

import os
from typing import Dict, Tuple
import pandas as pd
import sys
import importlib.util

def import_indelphi():
    """Import inDelphi module from the correct location."""
    # プロジェクトルートの絶対パスを取得
    current_dir = os.path.dirname(os.path.abspath(__file__))
    indelphi_path = os.path.join(current_dir, 'inDelphi-model', 'inDelphi.py')
    
    print(f"Attempting to import inDelphi from: {indelphi_path}")  # デバッグ用
    
    if not os.path.exists(indelphi_path):
        raise ImportError(f"inDelphi.py not found at {indelphi_path}")

    spec = importlib.util.spec_from_file_location("inDelphi", indelphi_path)
    indelphi = importlib.util.module_from_spec(spec)
    sys.modules["inDelphi"] = indelphi  # グローバルにモジュールを登録
    spec.loader.exec_module(indelphi)

    # モデルパスを設定
    model_dir = os.path.join(current_dir, 'inDelphi-model', 'model-sklearn-0.20.0')
    if not hasattr(indelphi, 'MODEL_DIR'):
        setattr(indelphi, 'MODEL_DIR', model_dir)
    
    print(f"Model directory set to: {model_dir}")  # デバッグ用
    return indelphi

# Import inDelphi
try:
    inDelphi = import_indelphi()
except Exception as e:
    print(f"Error importing inDelphi: {str(e)}")
    raise

def get_indelphi_predictions(upstream: str, downstream: str) -> Dict[str, list]:
    """
    Get InDelphi predictions for given sequences.
    
    Args:
        upstream: 50bp upstream sequence
        downstream: 50bp downstream sequence
    
    Returns:
        Dictionary containing prediction results
    """
    try:
        # モデルの初期化
        inDelphi.init_model(celltype='mESC')
        
        # 配列の結合とカットサイトの設定
        sequence = upstream + downstream
        cutsite = len(upstream)
        
        # 予測の実行
        pred_df, stats = inDelphi.predict(sequence, cutsite)
        pred_df = inDelphi.add_genotype_column(pred_df, stats)
        
        # 結果のソートと上位3件の抽出
        pred_df_sorted = pred_df.sort_values(by='Predicted frequency', ascending=False)
        
        return {
            "Top3 patterns": pred_df_sorted.head(3)['Category'].tolist(),
            "Top3 frequencies": pred_df_sorted.head(3)['Predicted frequency'].tolist(),
            "Top3 categories": pred_df_sorted.head(3)['Category'].tolist(),
            "Top3 positions": pred_df_sorted.head(3)['Genotype position'].tolist(),
            "Top3 lengths": pred_df_sorted.head(3)['Length'].tolist(),
            "Top3 genotypes": pred_df_sorted.head(3)['Genotype'].tolist(),
            "stats": stats
        }
        
    except Exception as e:
        print(f"Error in InDelphi prediction: {str(e)}")
        import traceback
        traceback.print_exc()
        return None