import pandas as pd
import numpy as np

import os

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

import joblib

# 定义模型（增加模块控制参数）
class LSTM_SelfAttention_Model(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, static_feature_size, 
                 lstm_dropout_rate=0.5, attention_dropout_rate=0.3,
                 use_time_aware=True, use_dnn=True):
        super(LSTM_SelfAttention_Model, self).__init__()
        # 控制参数
        self.use_time_aware = use_time_aware
        self.use_dnn = use_dnn
        
        # LSTM层
        self.lstm_workday = nn.ModuleList([nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True) for _ in range(input_size)])
        self.lstm_restday = nn.ModuleList([nn.LSTM(input_size=1, hidden_size=hidden_size, batch_first=True) for _ in range(input_size)])
        
        # Dropout层
        self.lstm_dropout = nn.Dropout(lstm_dropout_rate)
        self.attention_dropout = nn.Dropout(attention_dropout_rate)
        
        # 自注意力层
        self.self_attention_workday = nn.ModuleList([nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2) for _ in range(input_size)])
        self.self_attention_restday = nn.ModuleList([nn.MultiheadAttention(embed_dim=hidden_size, num_heads=2) for _ in range(input_size)])
        
        # 时间感知模块（可选）
        if self.use_time_aware:
            self.time_aware_module = nn.Sequential(
                nn.Linear(1, 16),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(16, 8),
                nn.ReLU(),
                nn.Dropout(0.5)
            )
        
        # DNN编码器（可选）
        if self.use_dnn:
            self.dnn = nn.Sequential(
                nn.Linear(static_feature_size, 64),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(64, 32),
                nn.ReLU(),
                nn.Dropout(0.5)
            )
        
        # 计算最终全连接层的输入维度
        fc_input_size = hidden_size * 2 * input_size  # 基础LSTM特征
        
        if self.use_dnn:
            fc_input_size += 32  # 添加DNN特征维度
        
        if self.use_time_aware:
            fc_input_size += 8  # 添加时间感知特征维度
            
        # 最终全连接层
        self.fc = nn.Linear(fc_input_size, output_size)

    def self_attention(self, lstm_out, self_attention_layer):
        lstm_out = lstm_out.permute(1, 0, 2)
        attn_output, attn_weights = self_attention_layer(lstm_out, lstm_out, lstm_out)
        attn_output = attn_output.permute(1, 0, 2)
        weighted_sum = torch.mean(attn_output, dim=1)
        return weighted_sum, attn_weights
    
    def forward(self, workday_input, restday_input, static_input, time_gap):
        # 存储注意力权重
        workday_attn_weights_list = []
        restday_attn_weights_list = []

        # 处理工作日分支
        workday_attended_list = []
        for i in range(workday_input.size(2)):
            workday_lstm_out, _ = self.lstm_workday[i](workday_input[:, :, i].unsqueeze(2))
            workday_lstm_out = self.lstm_dropout(workday_lstm_out)
            workday_attended, attn_weights = self.self_attention(workday_lstm_out, self.self_attention_workday[i])
            workday_attended = self.attention_dropout(workday_attended)
            workday_attended_list.append(workday_attended)
            workday_attn_weights_list.append(attn_weights)

        # 处理休息日分支
        restday_attended_list = []
        for i in range(restday_input.size(2)):
            restday_lstm_out, _ = self.lstm_restday[i](restday_input[:, :, i].unsqueeze(2))
            restday_lstm_out = self.lstm_dropout(restday_lstm_out)
            restday_attended, attn_weights = self.self_attention(restday_lstm_out, self.self_attention_restday[i])
            restday_attended = self.attention_dropout(restday_attended)
            restday_attended_list.append(restday_attended)
            restday_attn_weights_list.append(attn_weights)

        # 拼接LSTM特征
        combined_features = torch.cat(workday_attended_list + restday_attended_list, dim=1)
        combined_features = self.lstm_dropout(combined_features)

        # 处理静态特征（如果启用DNN）
        if self.use_dnn:
            static_features = self.dnn(static_input)
            combined_features = torch.cat((combined_features, static_features), dim=1)

        # 处理时间差特征（如果启用时间感知）
        if self.use_time_aware:
            time_gap_features = self.time_aware_module(time_gap)
            combined_features = torch.cat((combined_features, time_gap_features), dim=1)

        # 最终预测
        output = self.fc(combined_features)

        return output, workday_attn_weights_list, restday_attn_weights_list
    
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# 读数据
data = pd.read_excel("Move.xlsx")

ids = data["Participant ID"].values
X = data.iloc[:, 0:192].values  # (N, 192)
print(X.shape)
X = X.reshape(len(X), 48, 4)
workday = X[:, :24, :]
restday = X[:, 24:, :]

# dummy 特征（与你训练时 use_dnn=False / use_time_aware=False 对齐）
static_feat = np.zeros((len(X), 1))
time_gap_feat = np.zeros((len(X), 1))

# 转 tensor
workday_tensor = torch.tensor(workday, dtype=torch.float32).to(device)
restday_tensor = torch.tensor(restday, dtype=torch.float32).to(device)
static_tensor = torch.tensor(static_feat, dtype=torch.float32).to(device)
time_gap_tensor = torch.tensor(time_gap_feat, dtype=torch.float32).to(device)

result_df = pd.DataFrame({"Participant ID": ids})

MODEL_DIR = "r2jn"

from torch.utils.data import TensorDataset, DataLoader

BATCH_SIZE = 32   # 如果还炸，降到 16 / 8
dataset = TensorDataset(
    workday_tensor,
    restday_tensor,
    static_tensor,
    time_gap_tensor
)

loader = DataLoader(
    dataset,
    batch_size=BATCH_SIZE,
    shuffle=False
)


for file in os.listdir(MODEL_DIR):
    if not file.endswith(".joblib"):
        continue

    disease = file[0:3]
    model_path = os.path.join(MODEL_DIR, file)

    print(f"🔍 推理疾病: {disease}")

    '''
    checkpoint = torch.load(
        model_path,
        map_location=device,
        weights_only=False
    )


    model = LSTM_SelfAttention_Model(
        input_size=4,
        hidden_size=32,
        output_size=1,
        static_feature_size=static_tensor.shape[1],
        use_time_aware=False,
        use_dnn=False
    ).to(device)

    model.load_state_dict(checkpoint["model_state_dict"])
    '''

    model = joblib.load(model_path)
    model = model.to(device)
    model.eval()

    all_probs = []

    model.eval()
    with torch.no_grad():
        for workday_b, restday_b, static_b, time_gap_b in loader:
            logits, _, _ = model(
                workday_b,
                restday_b,
                static_b,
                time_gap_b
            )
            probs = torch.sigmoid(logits)
            all_probs.append(probs.cpu())

    probs = torch.cat(all_probs).numpy().squeeze()

    # 加一列
    result_df[disease] = probs

result_df.to_csv("scores5.csv", index=False)
print("✅ 已生成 all_disease_risk_scores.csv")