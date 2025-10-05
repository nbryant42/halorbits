import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

moon_radius = 1737
earth_radius = 6371

# km^3/s^2
gm_sun = 1.32712440018e11
gm_earth = 3.986004418e5

def make_sphere(center, radius, color, name):
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    xs = radius * np.cos(u) * np.sin(v) + center[0]
    ys = radius * np.sin(u) * np.sin(v) + center[1]
    zs = radius * np.cos(v) + center[2]
    return go.Surface(x=xs, y=ys, z=zs,
                      colorscale=[[0, color],[1, color]],
                      showscale=False, name=name, opacity=0.8)


def load_kernels(kernels):
    for k in kernels:
        sp.furnsh(k)


def sample_spacecraft(sc_id, et_start, period, npts):
    """
    Sample the spacecraft and return coordinates in the reference GSE (Geocentric Solar Ecliptic) rotating-frame
    """
    targ = str(sc_id)
    ets = np.linspace(et_start, et_start + period, npts)
    x, y, z = [], [], []
    for et in ets:
        st_sc, _ = sp.spkezr(targ, et, "GSE", "NONE", "EARTH")
        r_plot = st_sc[:3]
        x.append(r_plot[0]); y.append(r_plot[1]); z.append(r_plot[2])

    return (x, y, z)


def main():
    # Use the DE405 ephemeris (version used by Genesis) and the provided Genesis SPK
    bsp_name = "gns_010811_041125_101231.bsp"
    kernels = [
        "naif0012.tls",
        "pck00010.tpc",
        "de405s.bsp",
        bsp_name,
        "gse.tk"
    ]

    load_kernels(kernels)

    sc_id = -47

    # Get coverage window from the genesis BSP
    cover = sp.spkcov(bsp_name, sc_id)
    et_start, et_end = sp.wnfetd(cover, 0)
    print("Coverage start (UTC):", sp.et2utc(et_start, 'C', 0))
    print("Coverage end   (UTC):", sp.et2utc(et_end, 'C', 0))

    # Primary mission ended at Earth return on 2004 SEP 08. BSP contains additional trajectory for the spacecraft bus
    # after it dropped the sample, swung past Earth, made a loop near L1, and finally fell out of the Earth's influence
    # into a heliocentric orbit. We can plot about as far as April-May 2005 before it becomes irrelevant.
    #et_end = sp.utc2et("2004 SEP 08 16:00:00")
    et_end = sp.utc2et("2005 MAY 01 00:00:00")

    # Estimate sampling rate based on ~5 halo cycles in 29.3 months (approx mission halo cycles)
    # At minimum, we want to sample ~400 points per cycle to get a smooth line, but a little more resolution near the
    # Earth flyby is nice, so we'll oversample a bit.
    seconds_per_month = 30.436875 * 86400.0
    total_cycles = (et_end - et_start) * 5.0 / (29.3 * seconds_per_month)
    total_pts = int(np.ceil(total_cycles * 400 * 64))

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Estimated cycles: {total_cycles:.2f}, total samples: {total_pts}")

    # Define mid-epoch as midpoint of coverage window
    et_mid = 0.5 * (et_start + et_end)

    (x, y, z) = sample_spacecraft(sc_id, et_start, (et_end - et_start), total_pts)

    # Mid-epoch static references: use the et_mid computed above
    st_moon_mid, _ = sp.spkezr("MOON", et_mid, "GSE", "NONE", "EARTH")

    # Earth at origin
    r_earth_plot = np.array([0.0, 0.0, 0.0])
    # compute static moon at mid-epoch
    r_moon_plot = st_moon_mid[:3]

    # --- sample one representative lunar orbital cycle around mission midpoint ---
    # Sample slightly more than one cycle so the plotted lunar trace closes cleanly,
    # with a small overlap to show precession.
    lunar_period_days = 31.0
    lunar_period = lunar_period_days * 86400.0
    n_lunar = 400 # set sampling resolution

    et_moon_start = et_mid - 0.5 * lunar_period
    et_moon_end = et_mid + 0.5 * lunar_period
    ets_moon = np.linspace(et_moon_start, et_moon_end, n_lunar)
    moon_x, moon_y, moon_z = [], [], []
    for etm in ets_moon:
        st_m, _ = sp.spkezr("MOON", etm, "GSE", "NONE", "EARTH")
        r_moon_i = st_m[:3]
        moon_x.append(r_moon_i[0])
        moon_y.append(r_moon_i[1])
        moon_z.append(r_moon_i[2])

    # Hill sphere approximation for Sun-Earth L1 and L2 points
    st_sun_mid, _ = sp.spkezr("SUN", et_mid, "GSE", "NONE", "EARTH")
    R_se = np.linalg.norm(st_sun_mid[:3])
    d_SE_L = R_se * (gm_earth / (3.0 * (gm_sun + gm_earth))) ** (1.0 / 3.0)

    # L1/L2 in the plotted (Earth-centered, rotated) frame: along +X from Earth
    l1 = np.array([-d_SE_L, 0.0, 0.0])
    l2 = np.array([d_SE_L, 0.0, 0.0])

    # Build plot
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='orange', width=4),
                               name="Genesis probe (GSE)"))
    # Moon and Earth spheres at their translated positions
    fig.add_trace(make_sphere(tuple(r_moon_plot), moon_radius, 'gray', "Moon"))
    fig.add_trace(make_sphere(tuple(r_earth_plot), earth_radius, 'blue', "Earth"))

    # Markers for Sun-Earth L1 and L2
    fig.add_trace(go.Scatter3d(x=[l1[0]], y=[l1[1]], z=[l1[2]], mode='markers+text', text=["L1 (approx)"],
                               textposition="top center", marker=dict(size=6, color='green'), name="L1 (approx)"))
    fig.add_trace(go.Scatter3d(x=[l2[0]], y=[l2[1]], z=[l2[2]], mode='markers+text', text=["L2 (approx)"],
                               textposition="top center", marker=dict(size=6, color='red'), name="L2 (approx)"))
    # Lunar orbital cycle trace (representative)
    fig.add_trace(go.Scatter3d(x=moon_x, y=moon_y, z=moon_z, mode='lines', line=dict(color='gray', width=2),
                               name='Moon orbit (1 cycle)'))

    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(range=[1750000, -2000000], title="x (km) [GSE rotating-frame, -X points toward Sun]"),
            yaxis=dict(range=[1250000, -1250000], title="y (km)"),
            zaxis_title="z (km)"),
        title="Genesis probe (including post-mission trajectory)")

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "genesis_halo_plot.html")
        png_path = os.path.join("docs", "genesis_halo_preview.png")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
        try:
            fig.write_image(png_path, width=1280, height=600)
            print(f"Wrote static preview PNG: {png_path}")
        except Exception as e:
            print("PNG export failed (do you have 'kaleido' installed?):", e)
    else:
        fig.show()


if __name__ == "__main__":
    main()
