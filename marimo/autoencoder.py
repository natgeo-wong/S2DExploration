import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # PyTorch Testing on ERA5 Data

    Just some pytorch tests, can we run an ANN or something to predict:

    $$x(t-n:t) \rightarrow x(t+1)$$
    """)
    return


@app.cell(hide_code=True)
def _():
    import os
    import marimo as mo
    import datetime as dt
    import netCDF4 as nc
    import numpy as np
    import ultraplot as uplt

    import torch
    import torch.nn as nn
    import torch.optim as optim

    import config
    return config, mo, nc, nn, np, optim, os, torch, uplt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Test Loading some Data
    """)
    return


@app.cell
def _():
    geoID = "BNF"
    return (geoID,)


@app.cell(hide_code=True)
def _(config, os):
    fol = os.path.join(config.datadir,"era5hr","climatology")
    return (fol,)


@app.cell(hide_code=True)
def _(fol, geoID, nc, os):
    ds = nc.Dataset(os.path.join(fol,
        geoID + "-t-19800101-20241231.nc"
    ))
    dtv = ds["valid_time"][:]
    pre = ds["pressures"][:]
    t = ds["t"][:,:]
    ds.close()
    ds = nc.Dataset(os.path.join(fol,
        geoID + "-w-19800101-20241231.nc"
    ))
    w = ds["w"][:,:]
    ds.close()
    ds = nc.Dataset(os.path.join(fol,
        geoID + "-q-19800101-20241231.nc"
    ))
    q = ds["q"][:,:]
    ds.close()
    return dtv, pre, q, t, w


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try Some PyTorch Things?
    """)
    return


@app.cell
def _(Autoencoder, nn):
    class SampleAE(nn.Module):
        def __init__(self, input_dim, latent_dim):
            super(Autoencoder, self).__init__()

            # Encoder: Takes the input (8 features) and compresses it to the latent space (3 features)
            self.encoder = nn.Sequential(
                nn.Linear(input_dim, 5),
                nn.ReLU(),
                nn.Linear(5, latent_dim) # The bottleneck layer
            )

            # Decoder: Takes the latent code (3 features) and reconstructs the original input (8 features)
            self.decoder = nn.Sequential(
                nn.Linear(latent_dim, 5),
                nn.ReLU(),
                nn.Linear(5, input_dim) # Output layer matches input dimension
            )

        def forward(self, x):
            # 1. Encoding: Convert input into a compressed latent representation
            encoded = self.encoder(x)
            # 2. Decoding: Reconstruct the input from the latent representation
            decoded = self.decoder(encoded)
            return decoded, encoded # Return both the reconstruction and the latent code
    return (SampleAE,)


@app.cell
def _(Autoencoder, input_dim, latent_dim, nn):
    class SampleEncoder(nn.Module):
        def __init__(self, input_dim, latent_dim):
            super(Autoencoder, self).__init__()

            # Encoder: Takes the input (8 features) and compresses it to the latent space (3 features)
            self.encoder = nn.Sequential(
                nn.Linear(input_dim, 5),
                nn.ReLU(),
                nn.Linear(5, latent_dim) # The bottleneck layer
            )

            # Decoder: Takes the latent code (3 features) and reconstructs the original input (8 features)
            self.decoder = nn.Sequential(
                nn.Linear(latent_dim, 5),
                nn.ReLU(),
                nn.Linear(5, input_dim) # Output layer matches input dimension
            )

        def forward(self, x):
            nn.Linear(input_dim, 5),
            nn.ReLU(),
            nn.Linear(5, latent_dim)
    return


@app.cell
def _():
    nt = 10000
    return (nt,)


@app.cell
def _(np, nt, q, t, torch, w):
    # --- 2. Create Dummy Data (A simple XOR-like classification dataset) ---
    # This data is NOT linearly separable, requiring a non-linear ANN to solve.
    X = torch.tensor(
        # ((w - np.mean(w,axis=0)) / np.std(w))[0:nt,:],
        np.concatenate((
            (t - np.mean(t,axis=0)) / np.std(t),
            (q - np.mean(q,axis=0)) / np.std(q),
        #     (w - np.mean(w,axis=0)) / np.std(w)
        ),axis=1)[0:nt,:],
        dtype=torch.float32
    )
    # Target labels (Classes 0 or 1)
    y = torch.tensor((w[1:(nt+1),:]- np.mean(w,axis=0)) / np.std(w), dtype=torch.float32)
    return X, y


@app.cell
def _(torch):
    # --- 2. Create Dummy Data (Unsupervised: Target y = Input X) ---
    # We use 8 features to demonstrate the dimensionality reduction to the latent space (3 features).
    INPUT_SIZE = 8
    LATENT_SIZE = 3 # The bottleneck dimension
    N_SAMPLES = 10
    EPOCHS = 1000

    # Create random input data
    torch.manual_seed(42) # For reproducible results
    return EPOCHS, INPUT_SIZE, LATENT_SIZE, N_SAMPLES


@app.cell
def _(INPUT_SIZE, N_SAMPLES, torch):
    X2 = torch.randn(N_SAMPLES, INPUT_SIZE, dtype=torch.float32)
    return (X2,)


@app.cell
def _(X2):
    X2
    return


@app.cell
def _(INPUT_SIZE, LATENT_SIZE, SampleAE):
    # Instantiate the model
    model = SampleAE(INPUT_SIZE, LATENT_SIZE)
    return (model,)


@app.cell
def _(model, nn, optim):
    # --- 3. Define Loss Function and Optimizer ---
    # BCEWithLogitsLoss is perfect for binary classification: it combines Sigmoid + Binary Cross Entropy
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.01)
    return criterion, optimizer


@app.cell
def _(EPOCHS, X, criterion, model, optimizer, torch, y):
    for epoch in range(EPOCHS):
        # 1. Forward Pass
        logits = model(X) # Logits are raw, unnormalized outputs

        # 2. Compute Loss
        loss = criterion(logits, y)

        # 3. Backward Pass & Optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # Log progress
        if (epoch + 1) % 10 == 0:
            # Calculate training accuracy
            predictions = torch.round(torch.sigmoid(logits)) # Convert logits to 0 or 1
            correct = (predictions == y).sum().item()
            accuracy = correct / y.size(0) * 100

            print(f"Epoch [{epoch+1}/{EPOCHS}], Loss: {loss.item():.4f}, Accuracy: {accuracy:.1f}%")
    return


@app.cell
def _(dtv, np, nt, pre, uplt, w, wtest2):
    uplt.close(); f2,a2 = uplt.subplots(nrows=3,aspect=3)

    a2[0].contourf(dtv[nt:(nt+1000)],pre,wtest2[0:1000,:].T,levels=np.arange(-.5,.5,0.1),extend="both")
    p2_1 = a2[0].panel("r",width=0.2)
    p2_1.plot(np.mean(w[(nt+1):,:],axis=0),pre)
    p2_1.plot(np.mean(wtest2,axis=0),pre)
    p2_2 = a2[0].panel("l",width=0.2)
    p2_2.plot(np.mean(wtest2,axis=0) - np.mean(w[(nt+1):,:],axis=0),pre)

    a2[1].contourf(dtv[nt:(nt+1000)],pre,w[nt:(nt+1000),:].T,levels=np.arange(-.5,.5,0.1),extend="both")

    a2[2].contourf(dtv[nt:(nt+1000)],pre,wtest2[0:1000,:].T - w[nt:(nt+1000),:].T,levels=np.arange(-.5,.5,0.1),extend="both")

    for ax2 in a2:
        ax2.format(ylim=(1000,1))
    p2_1.format(xlim=(-0.05,0.05))

    f2
    return


if __name__ == "__main__":
    app.run()
