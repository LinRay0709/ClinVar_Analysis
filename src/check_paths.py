# src/check_paths.py
import config  # 這會匯入你剛剛寫的 config.py

print("="*30)
print("路徑檢查工具")
print("="*30)
print(f"專案根目錄 (BASE_DIR):\n  {config.BASE_DIR}")
print("-" * 20)
print(f"BED 檔案路徑:\n  {config.BED_FILE}")
print("-" * 20)
print(f"VCF 檔案路徑:\n  {config.VCF_FILE}")
print("-" * 20)
print(f"輸出檔案路徑:\n  {config.OUTPUT_FILE}")
print("="*30)

# 檢查檔案是否真的存在
import os
print("檔案存在檢查:")
print(f"BED 存在?  {os.path.exists(config.BED_FILE)}")
print(f"VCF 存在?  {os.path.exists(config.VCF_FILE)}")