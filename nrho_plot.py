import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

moon_radius = 1737
earth_radius = 6371

def make_sphere(center, radius, color, name):
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    xs = radius * np.cos(u) * np.sin(v) + center[0]
    ys = radius * np.sin(u) * np.sin(v) + center[1]
    zs = radius * np.cos(v) + center[2]
    return go.Surface(x=xs, y=ys, z=zs,
                      colorscale=[[0, color],[1, color]],
                      showscale=False, name=name, opacity=0.8)

def get_embr_transform(et):
    """
    Compute the rotation matrix from J2000 to the EMBR frame at epoch `et`.

    EMBR definition:
      +X = unit vector from EMB (Earth-Moon Barycenter) to MOON (radial)
      +Z = unit vector normal to orbital plane and angular momentum (r × v)
      +Y = Z × X

    Parameters
    ----------
    et : float
        Ephemeris time (seconds past J2000 TDB)

    Returns
    -------
    xfm : (3,3) ndarray
        Rotation matrix that transforms J2000 vectors into EMBR coordinates.
    """

    # State of Moon wrt Earth-Moon barycenter in J2000
    st_moon, _ = sp.spkezr("MOON", et, "J2000", "NONE", "EMB")
    r = st_moon[:3]
    v = st_moon[3:]

    # Normalize helper
    def vhat(vec):
        return vec / np.linalg.norm(vec)

    # Basis vectors
    x_hat = vhat(r)           # radial (EMB → Moon)
    z_hat = vhat(np.cross(r, v))  # cross product guarantees normal to both X and velocity
    y_hat = vhat(np.cross(z_hat, x_hat))  # completes right-handed system

    # Assemble rotation matrix: rows = EMBR basis in J2000
    xfm = np.vstack([x_hat, y_hat, z_hat])

    return xfm

sp.furnsh("naif0012.tls")
sp.furnsh("pck00010.tpc")
sp.furnsh("de432s.bsp")
sp.furnsh("receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp")

ids = sp.spkobj("receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp")
print("Bodies in kernel:", ids)

cover = sp.spkcov("receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp", ids[0])
start, end = sp.wnfetd(cover, 0)
print("Coverage (ET):", start, "to", end)
print("Coverage (UTC):", sp.et2utc(start, "C", 0), "to", sp.et2utc(end, "C", 0))

# ----- sample ~one NRHO orbit from the first coverage window -----
cover = sp.spkcov("receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp", -60000)
et_start, _ = sp.wnfetd(cover, 0)

period = 6.7 * 86400.0   # hard-coded NRHO period
npts   = 400             # enough for a smooth line
ets    = np.linspace(et_start, et_start + period, npts)
et_mid = 0.5 * (ets[0] + ets[-1])

print("Sampling NRHO trajectory:")
print("  Start UTC:", sp.et2utc(ets[0], "C", 0))
print("  End UTC:  ", sp.et2utc(ets[-1], "C", 0))
print(f"  Samples:   {npts}")

# --------- build Moon-centric EMBR coordinates ---------
x, y, z = [], [], []
for et in ets:
    # J2000 states wrt EMB (Earth-Moon barycenter)
    st_gw, _    = sp.spkezr("-60000", et, "J2000", "NONE", "EMB")
    st_moon, _  = sp.spkezr("MOON",   et, "J2000", "NONE", "EMB")

    # EMBR transform
    xfm = get_embr_transform(et)

    # Positions in EMBR
    r_gw_embr   = xfm @ st_gw[:3]
    r_moon_embr = xfm @ st_moon[:3]     # ~[R, ~0, ~0]

    # Moon-centric (subtract Moon position)
    r_mc = r_gw_embr - r_moon_embr

    x.append(r_mc[0]); y.append(r_mc[1]); z.append(r_mc[2])

# --------- static references for Earth and L2 (Moon-centric) ---------
# Use the mid-epoch to place a representative Earth position
st_earth_mid, _ = sp.spkezr("EARTH", et_mid, "J2000", "NONE", "EMB")
xfm_mid         = get_embr_transform(et_mid)
r_earth_embr    = xfm_mid @ st_earth_mid[:3]
st_moon_mid, _  = sp.spkezr("MOON",  et_mid, "J2000", "NONE", "EMB")
r_moon_embr_mid = xfm_mid @ st_moon_mid[:3]

print("Angular momentum Z sign:", np.cross(st_moon_mid[:3], st_moon_mid[3:])[2])

# Moon-centric Earth vector (draw Earth once, static)
earth_center_mc = r_earth_embr - r_moon_embr_mid

# L2 marker: fixed CR3BP estimate from the Moon (along +X in Moon-centric EMBR)
mu      = 0.012150585609624
R_moon  = np.linalg.norm(r_moon_embr_mid)        # EMB→Moon distance at mid epoch
r_L2    = R_moon * (mu/3.0)**(1.0/3.0)           # Moon→L2 distance
l2_mc   = (r_L2, 0.0, 0.0)

# --------- PLOT ---------
fig = go.Figure()

# Real NRHO trajectory (Moon-centric EMBR)
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z, mode='lines',
    line=dict(color='green', width=4),
    name="NRHO (Moon-centric EMBR)"
))

# Moon at origin (static sphere)
fig.add_trace(make_sphere((0.0, 0.0, 0.0), moon_radius, 'gray', "Moon"))

# Earth as a static sphere at its mid-epoch Moon-centric position
fig.add_trace(make_sphere(tuple(earth_center_mc), earth_radius, 'blue', "Earth"))

# L2 as a fixed reference dot along +X
fig.add_trace(go.Scatter3d(
    x=[l2_mc[0]], y=[l2_mc[1]], z=[l2_mc[2]],
    mode='markers+text',
    text=["L2 (approx)"], textposition="top center",
    marker=dict(size=6, color='red'),
    name="L2 (approx)"
))

# +Z arrow (north) in Moon-centric EMBR coords
fig.add_trace(go.Scatter3d(
    x=[0, 0], y=[0, 0], z=[0, 5000],
    mode="lines+text",
    text=["", "North (+Z)"], textposition="top center",
    line=dict(width=6),
    name="North"
))

fig.update_layout(
    scene=dict(aspectmode='data',
               xaxis=dict(range=[100000, -410000], title="x (km) [Moon-centric EMBR]"),
               yaxis_title="y (km)",
               zaxis_title="z (km)"),
    scene_camera=dict(
        eye=dict(x=0., y=3.0, z=0.5)
    ),
    title="Lunar Gateway NRHO — Moon-centric, EMBR (Earth-Moon Barycentric Rotating frame of reference) axes (L2 marker approx)"
)

# ---------------------------------------------
# Output control
if "--static" in sys.argv:
    os.makedirs("docs", exist_ok=True)
    html_path = os.path.join("docs", "nrho_plot.html")
    png_path = os.path.join("docs", "nrho_preview.png")

    # Write self-contained HTML
    fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
    print(f"Wrote interactive HTML: {html_path}")

    # Write PNG preview (requires kaleido)
    try:
        fig.write_image(png_path, width=1280, height=600)
        print(f"Wrote static preview PNG: {png_path}")
    except Exception as e:
        print("PNG export failed (do you have 'kaleido' installed?):", e)

else:
    # Default: launch one-shot server and open in browser
    fig.show()