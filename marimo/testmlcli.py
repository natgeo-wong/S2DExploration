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
def _(dtv, pre, q, t, uplt, w):
    uplt.close(); fig,axs = uplt.subplots(nrows=3,aspect=2)

    axs[0].contourf(dtv[0:500],pre,w[0:500,:].T)
    axs[1].contourf(dtv[0:500],pre,t[0:500,:].T)
    axs[2].contourf(dtv[0:500],pre,q[0:500,:].T)

    for ax in axs:
        ax.format(yscale="log",ylim=(1000,1))

    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try Some PyTorch Things?
    """)
    return


@app.cell
def _(nn):
    class TestTorch(nn.Module):
        def __init__(self, input_size, hidden_size, output_size):
            super(TestTorch, self).__init__()
            # 1st Layer: Maps inputs to a hidden layer of size 'hidden_size'
            self.fc1 = nn.Linear(input_size, hidden_size)
            # Non-linear activation function (essential for an ANN!)
            self.relu = nn.ReLU()
            # 2nd Layer (Output): Maps the hidden layer to the final output (1 logit for binary classification)
            self.fc2 = nn.Linear(hidden_size, output_size)

        def forward(self, x):
            # Pass through the first layer, apply ReLU, and then pass to the output layer
            out = self.fc1(x)
            out = self.relu(out)
            out = self.fc2(out)
            return out
    return (TestTorch,)


@app.cell
def _(np, q, t, torch, w):
    # --- 2. Create Dummy Data (A simple XOR-like classification dataset) ---
    # This data is NOT linearly separable, requiring a non-linear ANN to solve.
    X = torch.tensor(
        np.concatenate((t,q,w),axis=1)[0:-1,:],
        dtype=torch.float32
    )
    # Target labels (Classes 0 or 1)
    y = torch.tensor(w[1:,:], dtype=torch.float32)
    return X, y


@app.cell
def _():
    # Parameters
    INPUT_SIZE = 37 * 3
    HIDDEN_SIZE = 10
    OUTPUT_SIZE = 37
    EPOCHS = 500
    return EPOCHS, HIDDEN_SIZE, INPUT_SIZE, OUTPUT_SIZE


@app.cell
def _(HIDDEN_SIZE, INPUT_SIZE, OUTPUT_SIZE, TestTorch):
    # Instantiate the model
    model = TestTorch(INPUT_SIZE, HIDDEN_SIZE, OUTPUT_SIZE)
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
def _(X, model):
    wtest = model(X)
    return (wtest,)


@app.cell
def _(np, wtest):
    np.size(wtest)
    return


@app.cell(disabled=True)
def _(X, model, torch, y):
    # Put model in evaluation mode (important for production/testing)
    model.eval() 
    with torch.no_grad(): # Disable gradient calculation for efficiency
        final_logits = model(X)
        final_predictions = torch.round(torch.sigmoid(final_logits))

        print(f"Input X:\n{X}")
        print(f"True Labels (y):\n{y.T}")
        print(f"Predicted Labels:\n{final_predictions.T}")

        # Test new input
        test_input = torch.tensor([[0.5, 0.5]], dtype=torch.float32)
        test_logit = model(test_input).item()
        test_prob = torch.sigmoid(torch.tensor(test_logit)).item()
        test_class = round(test_prob)

        print("\n--- Prediction for [0.5, 0.5] ---")
        print(f"Logit: {test_logit:.2f}")
        print(f"Probability of Class 1: {test_prob:.2f}")
        print(f"Predicted Class: {test_class}")
    return


if __name__ == "__main__":
    app.run()
