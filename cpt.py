import pygmt

pygmt.makecpt(
    cmap="polar",
    series=[-200, 200, 50],
    continuous=False,
    reverse=False,
    output='polar_geoid.cpt'
)