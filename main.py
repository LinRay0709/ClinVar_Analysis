# main.py (在專案根目錄)

from src import processing

def main():
    try:
        processing.process_bed_regions()
    except KeyboardInterrupt:
        print("\n使用者中斷執行 (KeyboardInterrupt)")
        sys.exit(0)
    except Exception as e:
        print(f"\n發生未預期的錯誤: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()