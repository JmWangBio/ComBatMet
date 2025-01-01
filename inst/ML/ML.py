
## Author: Junmin Wang
## Date: December 31st, 2024

# Import libraries
import numpy as np
import pandas as pd
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import random

# Read the methylation data from RDS files
NT_raw_dat = pyreadr.read_r("/path/to/NT_raw_dat.rds")
TP_raw_dat = pyreadr.read_r("/path/to/TP_raw_dat.rds")
NT_adj_dat = pyreadr.read_r("/path/to/NT_ComBat_met_dat.rds")
TP_adj_dat = pyreadr.read_r("/path/to/TP_ComBat_met_dat.rds")

# Extract the data matrices from the RDS files
normal_raw_data = NT_raw_dat[None]
cancer_raw_data = TP_raw_dat[None]
normal_adj_data = NT_adj_dat[None]
cancer_adj_data = TP_adj_dat[None]

# Transpose the matrices so rows are samples and columns are features
normal_raw_data = normal_raw_data.transpose()
cancer_raw_data = cancer_raw_data.transpose()
normal_adj_data = normal_adj_data.transpose()
cancer_adj_data = cancer_adj_data.transpose()

# Combine normal and cancer data into a single DataFrame
all_raw_data = pd.concat([normal_raw_data, cancer_raw_data], axis = 0)
all_adj_data = pd.concat([normal_adj_data, cancer_adj_data], axis = 0)
print(f"Combined data shape: {all_raw_data.shape}")

# Load metadata
metadata = pd.read_csv("/path/to/metadata.csv")
print(metadata.head())

# Ensure the sample IDs in the metadata match the columns in the data
assert set(metadata["barcode"]) == set(all_raw_data.index), "Sample IDs do not match raw data!"
assert set(metadata["barcode"]) == set(all_adj_data.index), "Sample IDs do not match adjusted data!"

# Replace 'NT' with 'normal' and 'TP' with 'cancer' in the metadata DataFrame
metadata['group'] = metadata['group'].replace({'NT': 'normal', 'TP': 'cancer'})

# Merge metadata with the data to assign labels
all_raw_data['group'] = metadata.set_index('barcode')['group']
all_adj_data['group'] = metadata.set_index('barcode')['group']
print(all_raw_data['group'].value_counts())
print(all_adj_data['group'].value_counts())

# Prepare features and labels
# Remove columns with any NAs in the data matrix
X = all_raw_data.drop(columns = ['group']).dropna(axis=1).values
X_adjusted = all_adj_data.drop(columns = ['group']).dropna(axis=1).values
y = (all_raw_data['group'] == 'cancer').astype(int).values

# Number of random features to select
num_random_features = 3

# Number of iterations to repeat the experiment
num_iterations = 50

# Set the random seed
random.seed(100)

# Function to build a simple neural network
def build_model(input_dim):
    model = Sequential([
        Dense(16, activation='relu', input_dim=input_dim),
        Dense(4, activation='relu'),
        Dense(1, activation='sigmoid')  # Binary classification
    ])
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    return model

# Pre-adjustment experiment
accuracies_pre_adjustment = []
for _ in range(num_iterations):
    # Randomly select 3 features
    random_features  = random.sample(range(X.shape[1]), num_random_features)
    X_random = X[:, random_features]
    
    # Split into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X_random, y, test_size=0.2, stratify=y, random_state=100)

    # Scale the data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train the model
    model = build_model(input_dim=num_random_features)
    model.fit(X_train_scaled, y_train, epochs=10, batch_size=16, verbose=0)
    
    # Evaluate the model
    y_pred = (model.predict(X_test_scaled) > 0.5).astype(int)
    accuracy = accuracy_score(y_test, y_pred)
    accuracies_pre_adjustment.append(accuracy)

# Post-adjustment experiment
accuracies_post_adjustment = []
for _ in range(num_iterations):
    # Randomly select 3 features
    random_features  = random.sample(range(X_adjusted.shape[1]), num_random_features)
    X_random = X_adjusted[:, random_features]
    
    # Split into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X_random, y, test_size=0.2, stratify=y, random_state=100)

    # Scale the data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train the model
    model = build_model(input_dim=num_random_features)
    model.fit(X_train_scaled, y_train, epochs=10, batch_size=16, verbose=0)
    
    # Evaluate the model
    y_pred = (model.predict(X_test_scaled) > 0.5).astype(int)
    accuracy = accuracy_score(y_test, y_pred)
    accuracies_post_adjustment.append(accuracy)

print(np.mean(accuracies_pre_adjustment), np.std(accuracies_pre_adjustment))
print(np.mean(accuracies_post_adjustment), np.std(accuracies_post_adjustment))

# Create a DataFrame with the accuracies
df_accuracies = pd.DataFrame({
    'Pre-Adjustment Accuracy': accuracies_pre_adjustment,
    'Post-Adjustment Accuracy': accuracies_post_adjustment
})

# Save the DataFrame to a CSV file
df_accuracies.to_csv('/path/to/accuracies_comparison.csv', index=False)
