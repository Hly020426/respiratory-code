import os
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score

# =========================
# Dataset
# =========================
class MotionDataset(Dataset):
    def __init__(self, X, y, ids):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32)
        self.ids = ids

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.ids[idx]


def reshape_input(batch_x):
    B = batch_x.shape[0]
    return batch_x.view(B, 24, 8)

class BetterTransformer(nn.Module):
    def __init__(self, input_dim=8, d_model=16, nhead=2, num_layers=1):
        super().__init__()

        self.proj = nn.Linear(input_dim, d_model)
        self.pos_embedding = nn.Parameter(torch.randn(1, 24, d_model))

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            batch_first=True
        )

        self.transformer = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_layers
        )

        self.attn_pool = nn.Sequential(
            nn.Linear(d_model, 1)
        )

        self.fc = nn.Sequential(
            nn.Linear(d_model,8),
            nn.ReLU(),
            nn.BatchNorm1d(8),
            nn.Linear(8, 1)
        )

    def forward(self, x):
        x = self.proj(x) + self.pos_embedding
        x = self.transformer(x)

        weights = torch.softmax(self.attn_pool(x), dim=1)
        x = torch.sum(weights * x, dim=1)

        return self.fc(x)

def test_model(disease, df, model_path):

    print(f"\n====== {disease} ======")

    feature_cols = [str(i) for i in range(192)]

    X = df[feature_cols].values
    y = df[disease].values
    ids = df["Participant ID"].values

    loader = DataLoader(
        MotionDataset(X, y, ids),
        batch_size=128,
        shuffle=False
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = BetterTransformer().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()

    all_probs, all_labels, all_ids = [], [], []

    with torch.no_grad():
        for batch_x, batch_y, batch_ids in loader:
            batch_x = batch_x.to(device)

            x = reshape_input(batch_x)
            outputs = model(x)

            probs = torch.sigmoid(outputs).squeeze(1).cpu().numpy()

            all_probs.extend(probs)
            all_labels.extend(batch_y.numpy())
            all_ids.extend(batch_ids.numpy())

    auc = roc_auc_score(all_labels, all_probs)
    print(f"✅ AUC: {auc:.4f}")

    print(f"score: {np.min(all_probs):.4f} ~ {np.max(all_probs):.4f}")

    os.makedirs("result_test", exist_ok=True)

    result_df = pd.DataFrame({
        "id": all_ids,
        "label": all_labels,
        "score": all_probs
    })

    save_path = f"result_all/Fall/{disease}.csv"
    result_df.to_csv(save_path, index=False)

    print(f"result: {save_path}")
'''
# =========================
# ceshi
# =========================
if __name__ == "__main__":

    target_diseases = ['F05','F32']

    wales_df = pd.read_csv("icd_wales.csv")
    scotland_df = pd.read_csv("icd_scotland.csv")
    label_df = pd.read_csv("all_icd_f_group_1000_10y_wear.csv")

    test_motion = pd.concat([wales_df, scotland_df])
    t_df = pd.merge(test_motion, label_df, on="Participant ID")

    feature_cols = [str(i) for i in range(192)]

    for disease in target_diseases:

        print(f"\n====== {disease} ======")

        test_df_d = t_df[t_df[disease] != 2].dropna(subset=feature_cols + [disease])

        model_path = f"model/final/transformer_{disease}.pth"

        test_model(disease, test_df_d, model_path)

'''
# =========================
# All
# =========================
if __name__ == "__main__":

    target_diseases = ['F00', 'F01', 'F17', 'F40']

    all_df = pd.read_csv("all_movement_features.csv")
    label_df = pd.read_csv("all_icd_f_group_1000_10y_wear.csv")

    t_df = pd.merge(all_df, label_df, on="Participant ID")

    feature_cols = [str(i) for i in range(192)]

    for disease in target_diseases:

        print(f"\n======  {disease} ======")

        test_df_d = t_df[t_df[disease] != 2].dropna(subset=feature_cols + [disease])

        model_path = f"model/final/transformer_{disease}.pth"

        test_model(disease, test_df_d, model_path)
