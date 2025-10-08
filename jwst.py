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
    # sample epochs
    ets = np.linspace(et_start, et_start + period, npts)
    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    for et in ets:
        st_sc, _ = sp.spkezr(targ, et, "GSE", "NONE", "EARTH")
        r_plot = st_sc[:3]
        v_plot = st_sc[3:6]
        x.append(r_plot[0]); y.append(r_plot[1]); z.append(r_plot[2])
        vx.append(v_plot[0]); vy.append(v_plot[1]); vz.append(v_plot[2])

    # return positions, epochs, and velocities (all lists)
    return (x, y, z, ets, vx, vy, vz)


def main():
    # Use the DE440 ephemeris and the provided JWST SPK
    bsp_name = "jwst_rec.bsp"
    kernels = [
        "naif0012.tls",
        "pck00010.tpc",
        "de440s.bsp",
        bsp_name,
        "gse.tk"
    ]

    load_kernels(kernels)

    sc_id = -170

    # Get coverage window from the BSP
    cover = sp.spkcov(bsp_name, sc_id)
    et_start, et_end = sp.wnfetd(cover, 0)
    print("Coverage start (UTC):", sp.et2utc(et_start, 'C', 0))
    print("Coverage end   (UTC):", sp.et2utc(et_end, 'C', 0))

    # One sample per day
    total_pts = int(np.ceil((et_end - et_start) / 86400.0))

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Total samples: {total_pts}")

    # Define mid-epoch as midpoint of coverage window
    et_mid = 0.5 * (et_start + et_end)

    (x, y, z, ets, vx, vy, vz) = sample_spacecraft(sc_id, et_start, (et_end - et_start), total_pts)

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

    # Compute an estimated L2 distance for each sampled epoch to capture variation in Earth-Sun distance.
    # Hill-sphere-based approximation is evaluated at each sample.
    R_se_list = []
    ets = np.linspace(et_start, et_end, total_pts)
    for et in ets:
        st_sun, _ = sp.spkezr("SUN", et, "GSE", "NONE", "EARTH")
        R_se = np.linalg.norm(st_sun[:3])
        R_se_list.append(R_se)

    R_se_arr = np.array(R_se_list)
    # Hill-sphere approx for L1/L2 distance from Earth along Sun-Earth line
    d_SE_L_arr = R_se_arr * (gm_earth / (3.0 * (gm_sun + gm_earth))) ** (1.0 / 3.0)

    # L2 positions (in plotted Earth-centered GSE frame): along +X from Earth
    l2_x = d_SE_L_arr
    l2_y = np.zeros_like(l2_x)
    l2_z = np.zeros_like(l2_x)

    # Compute mean L2 location to use as the single marker
    l2_mean = np.array([np.mean(l2_x), 0.0, 0.0])

    # Build plot
    fig = go.Figure()

    # Prepare customdata for hover: timestamp (UTC string) and velocity magnitude (km/s)
    # Convert ET epochs to UTC strings for hover display
    utc_times = [sp.et2utc(et, 'C', 0) for et in ets]
    # Compute speed relative to Earth in km/s (GSE frame) as norm of velocity vector
    speeds = [np.linalg.norm([vx[i], vy[i], vz[i]]) for i in range(len(vx))]

    # customdata columns: utc string and speed (float).
    # Use a list of tuples (utc_string, speed) so we keep the string and numeric types separate.
    customdata = list(zip(utc_times, speeds))

    hovertemplate = (
        "UTC: %{customdata[0]}<br>"
        "Speed (in rotating frame, km/s): %{customdata[1]:.3f}<br>"
        "x: %{x:.1f} km<br>"
        "y: %{y:.1f} km<br>"
        "z: %{z:.1f} km<br>"
        "<extra></extra>"
    )

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='orange', width=4),
                               name="James Webb Space Telescope (GSE)",
                               customdata=customdata,
                               hovertemplate=hovertemplate))
    # Moon and Earth spheres at their translated positions
    fig.add_trace(make_sphere(tuple(r_moon_plot), moon_radius, 'gray', "Moon"))
    fig.add_trace(make_sphere(tuple(r_earth_plot), earth_radius, 'blue', "Earth"))

    # Marker for Sun-Earth L2
    l2_color = 'red'

    fig.add_trace(go.Scatter3d(x=l2_x, y=l2_y, z=l2_z, mode='lines',
                               line=dict(color=l2_color, width=1),
                               name='L2 (estimated range)'))

    # Plot mean L2 location as a single marker
    fig.add_trace(go.Scatter3d(x=[l2_mean[0]], y=[l2_mean[1]], z=[l2_mean[2]], mode='markers+text', text=["L2 (mean)"],
                               textposition="top center", marker=dict(size=6, color=l2_color), name="L2 (mean)"))
    # Lunar orbital cycle trace (representative)
    fig.add_trace(go.Scatter3d(x=moon_x, y=moon_y, z=moon_z, mode='lines', line=dict(color='gray', width=2),
                               name='Moon orbit (1 cycle)'))

    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(title="x (km) [GSE rotating-frame, -X points toward Sun]"),
            yaxis=dict(title="y (km)"),
            zaxis_title="z (km)"),
        title="James Webb Space Telescope (as flown through {et_end})".format(et_end=sp.et2utc(et_end, 'C', 0)))

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "jwst_plot.html")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
