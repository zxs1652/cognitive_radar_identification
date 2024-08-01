import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from scipy.io import loadmat


# Define the autoencoder model
# Define the CNN-based autoencoder model
class CNN_Autoencoder(nn.Module):
    def __init__(self):
        super(CNN_Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Conv1d(2, 16, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(16, 8, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(2)
        )
        self.decoder = nn.Sequential(
            nn.ConvTranspose1d(8, 16, kernel_size=2, stride=2),
            nn.ReLU(),
            nn.ConvTranspose1d(16, 2, kernel_size=2, stride=2),
            nn.Sigmoid()
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded


# # Generate sample data (300 pairs of 1x300 input and 1x300 output)
# np.random.seed(0)
# train_data_input = np.random.rand(1, 1, 300)  # 1x300 input features
# train_data_output = np.random.rand(1, 1, 300)  # 1x300 output features
#
# # Combine input and output for training the autoencoder
# train_data = np.concatenate((train_data_input, train_data_output), axis=1)  # Combined 2x300 data
# train_data = torch.tensor(train_data, dtype=torch.float32)
#
# # Generate new data (2x100)
# new_data_input = np.random.rand(1, 1, 100)
# new_data_output = np.random.rand(1, 1, 100)
# new_data = np.concatenate((new_data_input, new_data_output), axis=1)  # Combined 2x100 new data
# new_data = torch.tensor(new_data, dtype=torch.float32)


# load data from .mat
prx1_data = loadmat('prx1_800_rng3.mat')
prx4_data = loadmat('prx1_800_rng5.mat')
prx6_data = loadmat('prx3_800_rng5.mat')
rho123_data = loadmat('rho_800_rng3.mat')
rho456_data = loadmat('rho_800_rng5.mat')
tarpos_data = loadmat('tarpos_800_rng3.mat')

# reformat the training data
i_dataset = 0
received_power_1 = prx1_data['prx1'][i_dataset][0]
rho123 = rho123_data['rho'][i_dataset][0]
rho123_normalized = rho123 ** 4
rho123_total = np.sum(rho123_normalized, axis=0)
rho123_normalized /= rho123_total

transmitted_power_1 = received_power_1 * (rho123 ** 2)
max_transmitted_power_1 = np.max(transmitted_power_1)
transmitted_power_1 /= max_transmitted_power_1

tar_pos = tarpos_data['target_pos'][i_dataset][0].T

# reformat the test data
received_power_4 = prx4_data['prx1'][i_dataset][0]
rho456 = rho456_data['rho'][i_dataset][0]
rho456_normalized = rho456 ** 4
rho456_total = np.sum(rho456_normalized, axis=0)
rho456_normalized /= rho456_total

transmitted_power_4 = received_power_4 * (rho456 ** 2)
max_transmitted_power_4 = np.max(transmitted_power_4)
transmitted_power_4 /= max_transmitted_power_4

receiver_power_6 = prx6_data['prx3'][i_dataset][0]
transmitted_power_6 = receiver_power_6 * (rho456 ** 2)
max_transmitted_power_6 = np.max(transmitted_power_6)
transmitted_power_6 /= max_transmitted_power_6

# Prepare the training data
iRadar = 0
train_data_input = rho123_normalized[iRadar].reshape(1, 1, -1)
train_data_output = transmitted_power_1[iRadar].reshape(1, 1, -1)
train_data = np.concatenate((train_data_input, train_data_output), axis=1)  # Combined 2x300 data
train_data = torch.tensor(train_data, dtype=torch.float32)

# Prepare the test data
new_data_input = rho456_normalized[iRadar].reshape(1, 1, -1)
new_data_output = transmitted_power_6[iRadar].reshape(1, 1, -1)
new_data = np.concatenate((new_data_input, new_data_output), axis=1)  # Combined 2x100 new data
new_data = torch.tensor(new_data, dtype=torch.float32)

# Initialize the model, loss function and optimizer
model = CNN_Autoencoder()
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Training the autoencoder
num_epochs = 50
batch_size = 16

train_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, shuffle=True)

for epoch in range(num_epochs):
    for data in train_loader:
        # Forward pass
        output = model(data)
        loss = criterion(output, data)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item():.4f}')

# Pad new data to match the training data size (2x300)
padding = (0, train_data.shape[2] - new_data.shape[2])
padded_new_data = torch.nn.functional.pad(new_data, padding, "constant", 0)

# Calculate reconstruction error for new data
model.eval()
with torch.no_grad():
    reconstructions = model(padded_new_data)
    unpadded_reconstructions = reconstructions[:, :, :new_data.shape[2]]  # Unpad the output to original size
    reconstruction_error = torch.sqrt(torch.mean((new_data - unpadded_reconstructions) ** 2, dim=[1, 2]))

# Determine a threshold for reconstruction error based on training data
with torch.no_grad():
    train_reconstructions = model(train_data)
    train_error = torch.sqrt(torch.mean((train_data - train_reconstructions) ** 2, dim=[1, 2]))
    train_error_std = torch.std(train_data - train_reconstructions, dim=[1, 2])
    threshold = train_error  # + 1.96 * train_error_std

# Check if new data is within the threshold
new_data_within_distribution = reconstruction_error < threshold

# Output results
print("Reconstruction Error for New Data:", reconstruction_error)
print("Threshold:", threshold)
print("New Data Within Distribution:", new_data_within_distribution)
