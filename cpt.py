import pygmt

# define figure configuration
pygmt.config(COLOR_BACKGROUND=31.875/31.875/255, COLOR_FOREGROUND=255/31.875/31.875)

pygmt.makecpt(
    cmap="polar",
    series=[-100, 100, 25],
    continuous=False,
    reverse=False,
    output='polar_geoid.cpt'
)