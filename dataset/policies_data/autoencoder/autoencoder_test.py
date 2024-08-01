import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# Generate sample data (300 pairs of 1x300 input and 1x300 output)
np.random.seed(0)
train_data_input = np.random.rand(300, 1, 300)  # 1x300 input features
train_data_output = np.random.rand(300, 1, 300)  # 1x300 output features

# Combine input and output for training the autoencoder
train_data = np.concatenate((train_data_input, train_data_output), axis=1)  # Combined 2x300 data
train_data = torch.tensor(train_data, dtype=torch.float32)

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

    print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# Generate new data (2x100)
new_data_input = np.random.rand(50, 1, 100)
new_data_output = np.random.rand(50, 1, 100)
new_data = np.concatenate((new_data_input, new_data_output), axis=1)  # Combined 2x100 new data
new_data = torch.tensor(new_data, dtype=torch.float32)

# Pad new data to match the training data size (2x300)
padding = (0, train_data.shape[2] - new_data.shape[2])
padded_new_data = torch.nn.functional.pad(new_data, padding, "constant", 0)

# Calculate reconstruction error for new data
model.eval()
with torch.no_grad():
    reconstructions = model(padded_new_data)
    unpadded_reconstructions = reconstructions[:, :, :new_data.shape[2]]  # Unpad the output to original size
    reconstruction_error = torch.mean((new_data - unpadded_reconstructions) ** 2, dim=[1, 2])

# Determine a threshold for reconstruction error based on training data
with torch.no_grad():
    train_reconstructions = model(train_data)
    train_error = torch.mean((train_data - train_reconstructions) ** 2, dim=[1, 2])
    threshold = torch.mean(train_error) + 2 * torch.std(train_error)

# Check if new data is within the threshold
new_data_within_distribution = reconstruction_error < threshold

# Output results
print("Reconstruction Error for New Data:", reconstruction_error)
print("Threshold:", threshold)
print("New Data Within Distribution:", new_data_within_distribution)
