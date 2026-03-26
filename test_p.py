import pandas as pd
import numpy as np

import os

diagnoses = pd.read_csv('data/test_p.csv')
movement_features = pd.read_csv('data/all_movement_features.csv')
validation = pd.read_csv('data/participant_country_classification.csv')
baseline_features =  pd.read_csv('data/Baseline_characteristics_6.csv') 
static_features = pd.read_csv('data/all_static_features.csv')
chronic_diseases = pd.read_csv('data/chronic_diseases_before.csv') 

MODEL_SAVE_DIR = 'result3'

# 分组
validation = validation[['Participant ID', 'Country']]

# 创建test列，将国家映射为数字编码
country_code_mapping = {
    'England': 1,
    'Scotland': 2,
    'Wales': 3
}

validation['test'] = validation['Country'].map(country_code_mapping)

# 获取movement_features中的Participant ID列表
valid_ids = movement_features['Participant ID'].unique()

# 筛选diagnoses数据，只保留存在于movement_features中的ID
diagnoses = diagnoses[diagnoses['Participant ID'].isin(valid_ids)].copy()

# 合并 `validation` 和 `diagnoses`，确保使用 `Participant ID` 进行合并，并保持 'right' 合并方式
diagnoses = pd.merge(validation[["Participant ID", 'test']], diagnoses, on='Participant ID', how='right')

# 合并 `validation` 和 `movement_features`，确保使用 `Participant ID` 进行合并，并保持 'right' 合并方式
movement_features = pd.merge(validation[["Participant ID", 'test']], movement_features, on='Participant ID', how='right')

# 合并 `validation` 和 `movement_features`，确保使用 `Participant ID` 进行合并，并保持 'right' 合并方式
static_features = pd.merge(chronic_diseases, static_features, on='Participant ID', how='right')

# 合并 `validation` 和 `movement_features`，确保使用 `Participant ID` 进行合并，并保持 'right' 合并方式
static_features = pd.merge(validation[["Participant ID", 'test']], static_features, on='Participant ID', how='right')

# 按 `Participant ID` 排序两者，确保排序一致
diagnoses = diagnoses.sort_values(by='Participant ID').reset_index(drop=True)
movement_features = movement_features.sort_values(by='Participant ID').reset_index(drop=True)
static_features = static_features.sort_values(by='Participant ID').reset_index(drop=True)

static_features = static_features.drop(columns=['Start time of wear'])

# 筛选movement_features中第三列后不全为NaN的ID
movement_non_nan_ids = movement_features.loc[
    ~movement_features.iloc[:, 2:].isna().all(axis=1), 
    'Participant ID'
]

# 筛选static_features中第三列后不全为NaN的ID
static_non_nan_ids = static_features.loc[
    ~static_features.iloc[:, 2:].isna().all(axis=1), 
    'Participant ID'
]

# 获取两个数据集中都有效的ID
valid_ids = set(movement_non_nan_ids) & set(static_non_nan_ids)

# 筛选两个数据集
movement_features = movement_features[
    movement_features['Participant ID'].isin(valid_ids)
]

static_features = static_features[
    static_features['Participant ID'].isin(valid_ids)
]

# 获取movement_features中的Participant ID列表
valid_ids = movement_features['Participant ID'].unique()

# 筛选diagnoses数据，只保留存在于movement_features中的ID
diagnoses = diagnoses[diagnoses['Participant ID'].isin(valid_ids)].copy()

import pandas as pd
from sklearn.impute import KNNImputer

# 假设 test_movement_features 已经是一个 pandas DataFrame
# 从第三列开始进行 KNN 插补
movement_features_to_impute = movement_features.iloc[:, 2:]  # 选择第三列及之后的列
static_features_to_impute = static_features.iloc[:, 2:]  # 选择第三列及之后的列

# 初始化 KNN 插补器，k 设置为 5（可以根据需要调整）
knn_imputer = KNNImputer(n_neighbors=5)

# 进行插补，返回插补后的数据
movement_features_imputed = movement_features.copy()
static_features_imputed = static_features.copy()
movement_features.iloc[:, 2:] = knn_imputer.fit_transform(movement_features_to_impute)
# static_features.iloc[:, 2:] = knn_imputer.fit_transform(static_features_to_impute)

time_gap = static_features[["Participant ID", 'test', 'Time difference_x']]
static_features = static_features.drop(columns=['Time difference_x'])

train_diagnoses = diagnoses[diagnoses['test'] == 1]
test_diagnoses = diagnoses[diagnoses['test'].isin([2, 3])]

# 假设 movement_features 是 memmap 对象
train_movement_features = movement_features.iloc[:, 2:][movement_features['test'] == 1].values.reshape(-1, 48, 4)
test_movement_features = movement_features.iloc[:, 2:][movement_features['test'].isin([2, 3])].values.reshape(-1, 48, 4)

train_static_features = static_features.iloc[:, 2:][static_features['test'] == 1]
test_static_features = static_features.iloc[:, 2:][static_features['test'].isin([2, 3])]

train_time_gap = time_gap.iloc[:, 2:][time_gap['test'] == 1]
test_time_gap = time_gap.iloc[:, 2:][time_gap['test'].isin([2, 3])]

train_static_features = train_static_features.to_numpy()
test_static_features = test_static_features.to_numpy()

train_time_gap = train_time_gap.to_numpy()
test_time_gap = test_time_gap.to_numpy()

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed
import numpy as np
import pandas as pd

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


# 修改后的训练函数，增加模块控制参数
def train_and_evaluate_model_lstm_multihead_attention(disease, use_time_aware=True, use_dnn=True):
    # 获取当前疾病的标签
    y = train_diagnoses[disease]
    valid_indices = y.isin([0, 1])
    y = y[valid_indices]

    # 获取特征数据
    X = train_movement_features[valid_indices]
    static_features = train_static_features[valid_indices] if use_dnn else torch.zeros((len(y), 1))  # 如果不使用DNN，传入空特征
    time_gap = train_time_gap[valid_indices] if use_time_aware else torch.zeros((len(y), 1))  # 如果不使用时间感知，传入0
    
    # 分割数据
    workday_data = X[:, :24, :]
    restday_data = X[:, 24:, :]

    # 转换为张量
    workday_tensor = torch.tensor(workday_data, dtype=torch.float32)
    restday_tensor = torch.tensor(restday_data, dtype=torch.float32)
    static_tensor = torch.tensor(static_features, dtype=torch.float32)
    time_gap_tensor = torch.tensor(time_gap, dtype=torch.float32)
    y_tensor = torch.tensor(y.values, dtype=torch.float32).unsqueeze(1)
    
    # 划分训练验证集
    workday_train, workday_val, restday_train, restday_val, static_train, static_val, time_gap_train, time_gap_val, y_train, y_val = train_test_split(
        workday_tensor, restday_tensor, static_tensor, time_gap_tensor, y_tensor, 
        test_size=0.2, random_state=42, stratify=y_tensor.numpy()
    )

    # 创建数据加载器
    train_dataset = TensorDataset(workday_train, restday_train, static_train, time_gap_train, y_train)
    train_loader = DataLoader(train_dataset, batch_size=128, shuffle=True)

    # 初始化模型（传入模块控制参数）
    model = LSTM_SelfAttention_Model(
        input_size=4, 
        hidden_size=32, 
        output_size=1, 
        static_feature_size=static_features.shape[1],  # 动态获取静态特征维度
        lstm_dropout_rate=0.5,
        attention_dropout_rate=0.3,
        use_time_aware=use_time_aware,
        use_dnn=use_dnn
    )
    
    # 训练配置
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # 训练循环
    model.train()
    for epoch in range(20):
        epoch_loss = 0
        for workday_inputs, restday_inputs, static_inputs, time_gap_inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs, _, _ = model(workday_inputs, restday_inputs, static_inputs, time_gap_inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        print(f"Epoch {epoch + 1}, Loss: {epoch_loss / len(train_loader)}")

    # 评估模型
    model.eval()
    with torch.no_grad():
        y_pred_proba, workday_attn_weights, restday_attn_weights = model(workday_val, restday_val, static_val, time_gap_val)
        y_pred_proba = y_pred_proba.squeeze()
        y_pred = (torch.sigmoid(y_pred_proba) > 0.5).float()
        
        auc = roc_auc_score(y_val.numpy(), torch.sigmoid(y_pred_proba).numpy())
        
        cm = confusion_matrix(y_val.numpy(), y_pred.numpy())
        TP = cm[1, 1]
        FN = cm[1, 0]
        sensitivity = TP / (TP + FN) if (TP + FN) != 0 else 0

    # ========= 10. 保存模型 =========
    os.makedirs(MODEL_SAVE_DIR, exist_ok=True)

    model_path = os.path.join(
        MODEL_SAVE_DIR,
        f"{disease}_auc_{auc:.4f}.pth"
    )

    torch.save({
        "model_state_dict": model.state_dict(),
        "auc": auc,
        "sensitivity": sensitivity,
        "disease": disease,
        "use_time_aware": use_time_aware,
        "use_dnn": use_dnn
    }, model_path)

    print(f"💾 模型已保存: {model_path}")

    return disease, auc, sensitivity, model, workday_attn_weights, restday_attn_weights


# 主执行函数（增加模块控制参数）
def run_training(disease_columns, use_time_aware=True, use_dnn=True, n_jobs=60):
    # 存储结果
    auc_scores = {}
    sensitivities = {}
    models = {}
    attention_weights = {}
    
    # 并行训练
    results = Parallel(n_jobs=n_jobs)(
        delayed(train_and_evaluate_model_lstm_multihead_attention)(
            disease, use_time_aware, use_dnn
        ) for disease in disease_columns
    )

    # 收集结果
    for disease, auc, sensitivity, model, workday_attn_weights, restday_attn_weights in results:
        auc_scores[disease] = auc
        sensitivities[disease] = sensitivity
        models[disease] = model
        attention_weights[disease] = {
            "workday_attn_weights": workday_attn_weights,
            "restday_attn_weights": restday_attn_weights
        }
        print(f"Model for {disease}: AUC = {auc:.4f}, Sensitivity = {sensitivity:.4f}")

    # 打印汇总结果
    print("\nFinal Results:")
    for disease in disease_columns:
        print(f"{disease}: AUC = {auc_scores[disease]:.4f}, Sensitivity = {sensitivities[disease]:.4f}")
    
    return models, attention_weights, auc_scores, sensitivities


# 使用示例
if __name__ == "__main__":
    # 假设已经定义了 train_diagnoses, train_movement_features, train_static_features, train_time_gap
    
    # 获取疾病列
    disease_columns = train_diagnoses.columns[2:]
    
    # 运行训练（可以自由选择启用哪些模块）
    models, attention_weights, auc_scores, sensitivities = run_training(
        disease_columns=disease_columns,
        use_time_aware=False,  # 启用时间感知模块
        use_dnn=False,        # 启用DNN模块
        n_jobs=60           # 并行任务数
    )