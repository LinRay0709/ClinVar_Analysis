# verify_fold_distribution.py
# 檢查 5 fold 索引檔內的良惡性比例分佈
import os
import sys
import pandas as pd

# 加入 src 目錄到路徑
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg

def check_fold_distribution():
    print("=" * 60)
    print("5-Fold 良惡性比例分佈檢查")
    print(f"索引檔目錄: {cfg.SPLIT_INFO_DIR}")
    print("=" * 60)
    
    results = []
    
    for i in range(1, 6):
        filename = f"training_{i}.csv"
        filepath = os.path.join(cfg.SPLIT_INFO_DIR, filename)
        
        if not os.path.exists(filepath):
            print(f"[警告] 找不到檔案: {filepath}")
            continue
        
        df = pd.read_csv(filepath)
        total = len(df)
        pos_count = (df['label'] == 1).sum()  # 惡性
        neg_count = (df['label'] == 0).sum()  # 良性
        pos_ratio = pos_count / total * 100 if total > 0 else 0
        neg_ratio = neg_count / total * 100 if total > 0 else 0
        
        results.append({
            'fold': i,
            'total': total,
            'pos': pos_count,
            'neg': neg_count,
            'pos_ratio': pos_ratio,
            'neg_ratio': neg_ratio
        })
        
        print(f"\nFold {i}: {filename}")
        print(f"  總筆數: {total}")
        print(f"  惡性 (label=1): {pos_count} ({pos_ratio:.2f}%)")
        print(f"  良性 (label=0): {neg_count} ({neg_ratio:.2f}%)")
    
    # 統計摘要
    if results:
        print("\n" + "=" * 60)
        print("摘要統計")
        print("=" * 60)
        
        all_pos_ratios = [r['pos_ratio'] for r in results]
        all_neg_ratios = [r['neg_ratio'] for r in results]
        
        print(f"惡性比例範圍: {min(all_pos_ratios):.2f}% ~ {max(all_pos_ratios):.2f}%")
        print(f"良性比例範圍: {min(all_neg_ratios):.2f}% ~ {max(all_neg_ratios):.2f}%")
        print(f"比例差異 (惡性): {max(all_pos_ratios) - min(all_pos_ratios):.2f}%")
        
        if max(all_pos_ratios) - min(all_pos_ratios) < 1.0:
            print("\n✅ 各 Fold 的良惡性比例一致！")
        else:
            print("\n⚠️ 各 Fold 的良惡性比例存在差異")
    
    return results

if __name__ == "__main__":
    check_fold_distribution()
