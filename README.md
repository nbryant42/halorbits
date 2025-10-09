# halorbits

This Python project implements interactive 3D visualizations of:

- The proposed [Near-Rectilinear Halo Orbit](https://en.wikipedia.org/wiki/Near-rectilinear_halo_orbit) for the
[Lunar Gateway](https://en.wikipedia.org/wiki/Lunar_Gateway).
- As-flown trajectory of the [Genesis probe](https://en.wikipedia.org/wiki/Genesis_(spacecraft))
- Planned trajectory of the [Lucy probe](https://en.wikipedia.org/wiki/Lucy_(spacecraft))
- Planned trajectory of the [JUICE probe](https://en.wikipedia.org/wiki/Jupiter_Icy_Moons_Explorer)
- As-flown trajectory of the [James Webb Space Telescope](https://en.wikipedia.org/wiki/James_Webb_Space_Telescope)

3D pan/tilt/zoom is supported via [Plotly](https://plotly.com/python/).

## Preview:

![Preview](docs/nrho_preview.png)

## Interactive versions:

* [Lunar Gateway](https://nbryant42.github.io/halorbits/nrho_plot.html)
(This works best in a desktop browser, in a maximized window, on a 16:9 monitor.)
* [Genesis probe](https://nbryant42.github.io/halorbits/genesis_halo_plot.html)
* [Lucy probe](https://nbryant42.github.io/halorbits/lucy_plot.html)
* Juice probe:
  - [Cruise phase](https://nbryant42.github.io/halorbits/juice_cruise.html)
  - [Jupiter phase](https://nbryant42.github.io/halorbits/juice_jupiter.html)
  - [Ganymede phase](https://nbryant42.github.io/halorbits/juice_ganymede.html)
* [James Webb Space Telescope](https://nbryant42.github.io/halorbits/jwst_plot.html)

## Local dev notes

SPICE kernels that fit within Github's 100M size limit have been committed.
Others will have to be downloaded from NASA or ESA.

Steps to use the scripts:

1. Install dependencies (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

2. Run the script (plot will open in browser):

```bash
python <script>.py
```

To write static HTML and, for the Lunar Gateway, the PNG preview (requires `kaleido`):

```bash
python <script>.py --static
```