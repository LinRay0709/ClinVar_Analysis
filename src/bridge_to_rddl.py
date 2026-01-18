# bridge_to_rddl.py
import os
import numpy as np
import config as cfg
import bridge_utils as utils

def main():
    print("="*50)
    print(f"啟動 RDDL 資料轉換橋梁 (模組化版)")
    print(f"任務名稱: {cfg.TASK_NAME}")
    print("="*50)

    # 1. 建立目錄結構
    for d in [cfg.POS_DIR, cfg.NEG_DIR, cfg.SPLIT_INFO_DIR]:
        os.makedirs(d, exist_ok=True)
        print(f"檢查目錄: {d}")

    # 2. 處理 Test Set
    test_df = utils.process_dataset('test.csv', 'Test Set')
    if test_df is not None:
        out_path = os.path.join(cfg.SPLIT_INFO_DIR, 'test.csv')
        test_df.to_csv(out_path, index=False, header=False)
        print(f"-> 已建立測試清單: {out_path}")

    # 3. 處理 Ensemble Set
    ensem_df = utils.process_dataset('ensem.csv', 'Ensemble Set')
    if ensem_df is not None:
        out_path = os.path.join(cfg.SPLIT_INFO_DIR, 'ensemble.csv')
        ensem_df.to_csv(out_path, index=False, header=False)
        print(f"-> 已建立驗證清單: {out_path}")

    # 4. 處理 Training Set 並切分 5-Fold
    train_df = utils.process_dataset('train.csv', 'Training Set')
    if train_df is not None:
        print("正在將 Training Set 切分為 5 份 (for RDDL 5-Fold CV)...")
        shuffled_train = train_df.sample(frac=1, random_state=42)
        folds = np.array_split(shuffled_train, 5)
        
        for i, fold_df in enumerate(folds):
            out_path = os.path.join(cfg.SPLIT_INFO_DIR, f'training_{i+1}.csv')
            fold_df.to_csv(out_path, index=False, header=False)
            print(f"  -> Fold {i+1}: 已建立")

    print("\n✅ 所有轉換作業完成！")

if __name__ == "__main__":
    main()