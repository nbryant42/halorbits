import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

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


def sample_object(sc_id, et_start, period, npts):
    """
    Sample the object and return coordinates in the HJB (Heliocentric Jupiter Barycenter) rotating-frame
    defined in `hjb.tk`.
    """
    targ = str(sc_id)
    # sample epochs
    ets = np.linspace(et_start, et_start + period, npts)
    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    for et in ets:
        st_sc, _ = sp.spkezr(targ, et, "HJB", "NONE", "SUN")
        r_plot = st_sc[:3]
        v_plot = st_sc[3:6]
        x.append(r_plot[0]); y.append(r_plot[1]); z.append(r_plot[2])
        vx.append(v_plot[0]); vy.append(v_plot[1]); vz.append(v_plot[2])

    # return positions, epochs, and velocities (all lists)
    return (x, y, z, ets, vx, vy, vz)


def compute_l4_l5_from_jupiter(jx, jy, jz, deg=60.0):
    """
    Given arrays (or lists) of Jupiter barycenter positions in the HJB frame (relative to Sun at origin),
    compute the corresponding L4 (leading, +60 deg) and L5 (trailing, -60 deg) Lagrange points
    by rotating the Sun->Jupiter vector about the +Z axis by +-60 degrees.

    Returns: (l4_x, l4_y, l4_z), (l5_x, l5_y, l5_z) as lists.
    """
    theta = np.deg2rad(deg)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    l4_x, l4_y, l4_z = [], [], []
    l5_x, l5_y, l5_z = [], [], []

    for xi, yi, zi in zip(jx, jy, jz):
        # rotate (xi, yi) by +theta for L4
        x4 = xi * cos_t - yi * sin_t
        y4 = xi * sin_t + yi * cos_t
        z4 = zi

        # rotate by -theta for L5
        x5 = xi * cos_t + yi * sin_t
        y5 = -xi * sin_t + yi * cos_t
        z5 = zi

        l4_x.append(x4); l4_y.append(y4); l4_z.append(z4)
        l5_x.append(x5); l5_y.append(y5); l5_z.append(z5)

    return (l4_x, l4_y, l4_z), (l5_x, l5_y, l5_z)


def main():
    # Use the DE440 ephemeris (version used by Lucy) and the provided Lucy SPK
    bsp_name = "lcy_211016_330402_240718_v1.bsp"
    kernels = [
        "naif0012.tls",
        "pck00010.tpc",
        "de440s.bsp",
        bsp_name,
        "hjb.tk"
    ]

    load_kernels(kernels)

    sc_id = -49

    # Get coverage window from the Lucy BSP
    cover = sp.spkcov(bsp_name, sc_id)
    et_start, et_end = sp.wnfetd(cover, 0)
    print("Coverage start (UTC):", sp.et2utc(et_start, 'C', 0))
    print("Coverage end   (UTC):", sp.et2utc(et_end, 'C', 0))

    # Plot one sample per day.  Lucy has a ~12 year mission, so this is about 4000 points.
    seconds_per_day = 86400.0
    total_days = (et_end - et_start) / seconds_per_day
    total_pts = int(np.ceil(total_days))

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Total samples: {total_pts}")

    x, y, z, ets, vx, vy, vz = sample_object(sc_id, et_start, (et_end - et_start), total_pts)
    earth_x, earth_y, earth_z, _, _, _, _ = sample_object("EARTH", et_start, (et_end - et_start), total_pts)
    jupiter_x, jupiter_y, jupiter_z, _, _, _, _ = sample_object("JUPITER_BARYCENTER", et_start, (et_end - et_start), total_pts)

    # Mission target bodies (from BSP 'Bodies:' section) to plot
    target_defs = [
        ("PATROCLUS BARYCENTER", "20000617", 'purple'),
        ("DONALDJOHANSON", "20052246", 'teal'),
        ("LEUCUS", "20011351", 'brown'),
        ("20152830", "20152830", 'crimson'),
        ("POLYMELE", "20015094", 'gold'),
        ("EURYBATES", "920003548", 'cyan'),
        ("ORUS", "20021900", 'lime')
    ]

    target_samples = {}
    for label, ident, color in target_defs:
        tx, ty, tz, _, _, _, _ = sample_object(ident, et_start, (et_end - et_start), total_pts)
        target_samples[label] = (tx, ty, tz, color)

    # Compute L4/L5 estimated positions from Jupiter barycenter samples
    # (computed by rotating the Sun->Jupiter vector by +-60 deg)
    (l4_x, l4_y, l4_z), (l5_x, l5_y, l5_z) = compute_l4_l5_from_jupiter(jupiter_x, jupiter_y, jupiter_z)

    # Build plot
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode='markers+text', text=["Sun"],
                               textposition="top center", marker=dict(size=6, color='yellow'), name="Sun"))

    # Prepare customdata for hover: timestamp (UTC string) and velocity magnitude (km/s)
    # Convert ET epochs to UTC strings for hover display
    utc_times = [sp.et2utc(et, 'C', 0) for et in ets]
    # Compute speed in km/s as norm of velocity vector
    speeds = [np.linalg.norm([vx[i], vy[i], vz[i]]) for i in range(len(vx))]

    # customdata columns: utc string and speed (float).
    # Use a list of tuples (utc_string, speed) so we keep the string and numeric types separate.
    customdata = list(zip(utc_times, speeds))

    hovertemplate = (
        "UTC: %{customdata[0]}<br>"
        "Speed (in rotating frame, km/s): %{customdata[1]:.3f}<br>"
        "x: %{x:.7s} km<br>"
        "y: %{y:.7s} km<br>"
        "z: %{z:.7s} km<br>"
        "<extra></extra>"
    )

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='magenta', width=4),
                               name="Lucy probe (GSE)", customdata=customdata, hovertemplate=hovertemplate))

    fig.add_trace(go.Scatter3d(x=earth_x, y=earth_y, z=earth_z, mode='lines', line=dict(color='gray', width=2),
                               name='Earth orbit'))
    
    fig.add_trace(go.Scatter3d(x=jupiter_x, y=jupiter_y, z=jupiter_z, mode='lines',
                                 line=dict(color='orange', width=2, dash='dash'),
                                 name='Jupiter system barycenter'))

    # Mark the mean (center) of the Jupiter barycenter motion range
    mean_jup = (np.mean(jupiter_x), np.mean(jupiter_y), np.mean(jupiter_z))
    fig.add_trace(go.Scatter3d(x=[mean_jup[0]], y=[mean_jup[1]], z=[mean_jup[2]], mode='markers+text',
                                marker=dict(size=6, color='orange'), text=['Jupiter system barycenter'],
                                textposition='top center', name='Jupiter system barycenter, mean'))
    # Mark representative points (mean positions) for L4/L5 and label them
    mean_l4 = (np.mean(l4_x), np.mean(l4_y), np.mean(l4_z))
    mean_l5 = (np.mean(l5_x), np.mean(l5_y), np.mean(l5_z))
    fig.add_trace(go.Scatter3d(x=[mean_l4[0]], y=[mean_l4[1]], z=[mean_l4[2]], mode='markers+text',
                                marker=dict(size=6, color='green'), text=['L4'], textposition='top center', name='L4'))
    fig.add_trace(go.Scatter3d(x=[mean_l5[0]], y=[mean_l5[1]], z=[mean_l5[2]], mode='markers+text',
                                marker=dict(size=6, color='blue'), text=['L5'], textposition='top center', name='L5'))
    fig.add_trace(go.Scatter3d(x=l4_x, y=l4_y, z=l4_z, mode='lines', line=dict(color='green', width=2, dash='dash'),
                            name='L4 (estimated)'))
    fig.add_trace(go.Scatter3d(x=l5_x, y=l5_y, z=l5_z, mode='lines', line=dict(color='blue', width=2, dash='dash'),
                            name='L5 (estimated)'))

    # Plot mission target orbits
    for label, (tx, ty, tz, color) in target_samples.items():
        fig.add_trace(go.Scatter3d(x=tx, y=ty, z=tz, mode='lines', line=dict(color=color, width=2),
                                   name=f'{label} orbit'))

    # Map identifiers to color for intercept markers
    ident_to_color = {ident: color for (_, ident, color) in target_defs}

    # Intercept dates (UTC) from Wikipedia table — plot each object's position at the encounter epoch
    intercepts = [
        ("EARTH", "2022-10-16T00:00:00", "Earth 2022-10-16"),
        ("20152830", "2023-11-01T00:00:00", "152830 Dinkinesh 2023-11-01"),
        ("EARTH", "2024-12-13T00:00:00", "Earth 2024-12-13"),
        ("20052246", "2025-04-20T00:00:00", "52246 Donaldjohanson 2025-04-20"),
        ("920003548", "2027-08-12T00:00:00", "Eurybates 2027-08-12"),
        ("20015094", "2027-09-15T00:00:00", "Polymele 2027-09-15"),
        ("20011351", "2028-04-18T00:00:00", "Leucus 2028-04-18"),
        ("20021900", "2028-11-11T00:00:00", "Orus 2028-11-11"),
        ("EARTH", "2030-12-26T00:00:00", "Earth 2030-12-26"),
        ("20000617", "2033-03-02T00:00:00", "Patroclus–Menoetius 2033-03-02")
    ]

    for ident, utc_str, label in intercepts:
        et = sp.str2et(utc_str)
        st, _ = sp.spkezr(str(ident), et, "HJB", "NONE", "SUN")
        r = st[:3]
        col = ident_to_color.get(str(ident), 'black')
        fig.add_trace(go.Scatter3d(x=[r[0]], y=[r[1]], z=[r[2]], mode='markers+text',
                                    marker=dict(size=6, color=col, symbol='circle'),
                                    text=[label], textposition='top center',
                                    name=f'{label}'))

    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(title="x (km) [Distance from Sun]"),
            yaxis=dict(title="y (km)"),
            zaxis_title="z (km)"),
        title="Lucy spacecraft trajectory (Jupiter rotating frame)")

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "lucy_plot.html")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
