# ClinVar RDDL 訓練診斷報告

> 日期：2026-02-09
> 分析對象：OS1-1, OS1-2 方法，1 fold 訓練結果

---

## 1. 資料概況

| 項目 | 數值 |
|------|------|
| Positive (pathogenic) 樣本總數 | 423 |
| Negative (benign) 樣本總數 | 56,231 |
| 正負比例 | 1:133（極度不平衡） |
| 每 fold 正樣本數 | ~67-68 |
| 每 fold 負樣本數 | ~8,997 |
| 5-fold CV 訓練集 (4 folds) | pos ~270 / neg ~35,988 |
| 驗證集 (1 fold) | pos 68 / neg 8,997 |
| 輸入形狀 | (1024, 8) — 1024bp 序列 × 8 channels |

## 2. 目前 Hyperparameters（OS1-1 與 OS1-2 完全相同）

```json
{
  "batch_size": 64,
  "dropout_rate": 0.3,
  "epochs": 100,
  "initial_lr": 0.001,
  "decay_rate": 0.995,
  "bottom_lr": 0.00001
}
```

## 3. 訓練結果

### 3.1 指標總覽

| 方法 | Train F1 | Val F1（原始） | Val F1（resampled） | Val AUC | Val Precision | Val Recall |
|------|----------|----------------|----------------------|---------|---------------|------------|
| OS1-1 (1:1 oversampling) | 0.928 | 0.052 | 0.565 | 0.724 | 0.028 | 0.462 |
| OS1-2 (1:2 oversampling) | 0.813 | 0.040 | 0.437 | 0.677 | 0.021 | 0.296 |

### 3.2 Prediction 層面（最終模型對驗證集的預測）

| 方法 | AUC | F1 | Accuracy | 驗證集 pos | 驗證集 neg |
|------|-----|-----|----------|-----------|-----------|
| OS1-1 | 0.724 | 0.046 | 0.857 | 68 | 8,997 |
| OS1-2 | 0.677 | 0.044 | 0.890 | 68 | 8,997 |

### 3.3 Learning Curve 觀察

- **OS1-1**: Training F1 從 0.53 持續上升至 0.93；Val F1 始終在 0.02~0.06 之間徘徊，無改善趨勢。
- **OS1-2**: Training F1 從 0.00 上升至 0.81；Val F1 同樣始終在 0.00~0.06 之間，無改善。
- 兩者的 Sampled Validation F1 分別為 0.565 和 0.437，顯示模型在平衡的驗證集上有一定區分能力，但在原始分佈下完全失效。

## 4. 問題診斷

### 4.1 嚴重 Overfitting

Training F1 高達 92%，真實 Validation F1 僅 ~5%，差距極大。

**根本原因**：OS1-1 將 ~270 個正樣本 oversample 到 ~36,000 個，每個正樣本被重複約 133 次。模型直接記住了這些樣本而非學到泛化特徵。

### 4.2 極端類別不平衡

1:133 的比例使得：
- Oversampling 方法需要極大量的重複，導致嚴重 overfitting
- 驗證集僅有 68 個正樣本，precision 極低（~0.02），意味著模型做出的正預測中絕大多數是錯的
- 高 accuracy（85-89%）是假象，全部預測為負也能達到 ~99% accuracy

### 4.3 模型架構過於簡單

```
Conv1D(64, k=3) → GlobalAvgPool → Dense(64) → Dropout(0.3) → Dense(2, softmax)
```

僅一層 Conv1D（kernel_size=3），對 1024bp 長序列的特徵提取能力不足，無法捕捉中長距離的序列模式。

### 4.4 Hyperparameter 問題

- `epochs=100`：從 learning curve 看，~20 epoch 後 validation 已停滯，後 80 epochs 只加劇 overfitting
- `dropout_rate=0.3`：正則化強度不足以對抗嚴重 overfitting
- `decay_rate=0.995`：衰減過慢，LR 在後期仍維持較高值

## 5. 調整建議

### 5.1 最優先：改用 Cost-Sensitive 方法

建議優先嘗試不做 resampling 的方法，避免大量重複正樣本導致的 overfitting：

1. **CW（Class Weight）**：用 class weight 加權 loss，不改變資料分佈
2. **FL-gamma-2（Focal Loss）**：自動降低 easy negative 的權重，聚焦 hard examples
3. **MFE（Mean False Error）**：專門為不平衡問題設計

### 5.2 Hyperparameter 調整

| 參數 | 原值 | 建議值 | 理由 |
|------|------|--------|------|
| `batch_size` | 64 | 32 | 較小 batch 增加梯度噪音，有助正則化；非 resampling 方法下更易採樣到正樣本 |
| `dropout_rate` | 0.3 | 0.5 | 加強正則化以對抗 overfitting |
| `epochs` | 100 | 50 | Validation 約 20 epoch 後停滯，減少無效訓練 |
| `initial_lr` | 0.001 | 0.0005 | 較低 LR 配合較少 epochs，學習更穩定 |
| `decay_rate` | 0.995 | 0.98 | 更快衰減，避免後期 overfitting |
| `bottom_lr` | 0.00001 | 0.00001 | 維持不變 |

### 5.3 模型架構建議（未來可修改 USER_model.py）

建議改為多層 CNN + BatchNorm：

```python
x = Conv1D(32, 7, activation='relu', padding='same')(inputs)
x = BatchNormalization()(x)
x = MaxPooling1D(4)(x)
x = Conv1D(64, 5, activation='relu', padding='same')(x)
x = BatchNormalization()(x)
x = GlobalAveragePooling1D()(x)
x = Dense(32, activation='relu')(x)
x = Dropout(dropout_rate)(x)
outputs = Dense(2, activation='softmax')(x)
```

### 5.4 建議實驗順序

1. 先用 **CW** 方法 + 新 hyperparameters 跑 1 fold 測試
2. 再用 **FL-gamma-2** + 同樣 hyperparameters 比較
3. 如果效果仍不理想，再修改模型架構
4. 最後回頭跑 OS1-1 做對照組

---

## 附錄：LR Schedule 細節

```
Epoch 1-20 (前 20%): cosine warmup，LR = initial_lr * cos(1 - epoch/20)
Epoch 21-90 (20%-90%): exponential decay，LR = initial_lr * decay_rate^(epoch-20)
Epoch 91-100 (後 10%): 固定為 bottom_lr
```
