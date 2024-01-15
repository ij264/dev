import pygmt
grid = pygmt.datasets.load_earth_relief(
    resolution="30s", region=[144.5, 145.5, 13, 14.5], registration="gridline"
)

fig = pygmt.Figure()
fig.grdimage(grid=grid, frame="a", projection="M10c", cmap="oleron")
fig.grdcontour(grid=grid, interval=500, annotation=1000)
fig.coast(shorelines="2p", land="lightgray")
fig.colorbar(frame=["a1000", "x+lElevation", "y+lm"])
fig.show()
