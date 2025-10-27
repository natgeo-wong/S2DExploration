import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # PyTorch Testing on ERA5 Data

    Just some pytorch tests, can we run an ANN or something to predict:

    $$x(t-n:t) \rightarrow x(t+1)$$
    """
    )
    return


@app.cell(hide_code=True)
def _():
    import os
    import marimo as mo
    import datetime as dt
    import netCDF4 as nc
    import ultraplot as uplt

    import torch
    import torch.nn as nn
    import torch.optim as optim

    import config
    return config, mo, nc, nn, os, uplt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Test Loading some Data""")
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
    mo.md(r"""## Try Some PyTorch Things?""")
    return


@app.cell
def _(nn):
    class ANNClassifier(nn.Module):
        def __init__(self, input_size, hidden_size, output_size):
            super(ANNClassifier, self).__init__()
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
    return


if __name__ == "__main__":
    app.run()
