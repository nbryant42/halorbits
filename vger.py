import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

# Voyager 1 trajectory based on voyager_1.ST+1991_a54418u.merged.bsp
# Voyager 2 trajectory based on voyager_2.ST+1992_m05208u.merged.bsp
#
# This is plotted as the Sun-referenced position in a custom frozen frame, "VGER"
# which has SATURN BARYCENTER on the X axis at the time of closest Voyager 1 approach.

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


def sample_object(sc_id, et_array):
    """Sample object positions and compute speed magnitude and UTC strings.

    et_array must be iterable of ET seconds. Returns (x,y,z,speeds,utcs)
    where utcs is a list of UTC strings suitable for hover display.
    """
    targ = str(sc_id)
    x, y, z, speeds, utcs = [], [], [], [], []
    for et in et_array:
        st_sc, _ = sp.spkezr(targ, et, "VGER", "NONE", "SUN")
        r = st_sc[:3]
        if r[1] <= 6e9: # filter out extreme Y values for plotting
            v = st_sc[3:6]
            speed = float(np.linalg.norm(v))
            x.append(r[0]); y.append(r[1]); z.append(r[2])
            speeds.append(speed)
            utcs.append(sp.et2utc(et, 'C', 0))
    return x, y, z, speeds, utcs


def main():
    bsp_name = "voyager_1.ST+1991_a54418u.merged.bsp"
    kernels = [
        "naif0012.tls",
        "voyager_2.ST+1992_m05208u.merged.bsp",
        "vger.tk",
        bsp_name
    ]

    load_kernels(kernels)

    sc_id = -31

    # Get coverage window from the BSP
    cover = sp.spkcov(bsp_name, sc_id)
    et_start, et_end = sp.wnfetd(cover, 0)
    print("Coverage start (UTC):", sp.et2utc(et_start, 'C', 0))
    print("Coverage end   (UTC):", sp.et2utc(et_end, 'C', 0))

    # We want 1 sample per day; enough to render smooth curves.
    ets = np.arange(et_start, et_end, 86400.0)
    total_pts = len(ets)

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Total samples: {total_pts}")

    (x, y, z, speeds, utcs) = sample_object(sc_id, ets)
    (v2_x, v2_y, v2_z, v2_speeds, v2_utcs) = sample_object(-32, ets)
    (v_x, v_y, v_z, _, _) = sample_object("VENUS", ets)
    (e_x, e_y, e_z, _, _) = sample_object("EARTH", ets)
    (j_x, j_y, j_z, _, _) = sample_object("JUPITER BARYCENTER", ets)
    (s_x, s_y, s_z, _, _) = sample_object("SATURN BARYCENTER", ets)
    (me_x, me_y, me_z, _, _) = sample_object("MERCURY", ets)
    (ma_x, ma_y, ma_z, _, _) = sample_object("MARS", ets)
    (u_x, u_y, u_z, _, _) = sample_object("URANUS BARYCENTER", ets)
    (n_x, n_y, n_z, _, _) = sample_object("NEPTUNE BARYCENTER", ets)
    (p_x, p_y, p_z, _, _) = sample_object("PLUTO BARYCENTER", ets)

    # Build plot
    fig = go.Figure()

    hovertemplate = (
        "UTC: %{customdata[0]}<br>"
        "Speed (km/s): %{customdata[1]:.3f}<br>"
        "x: %{x:.7s} km<br>"
        "y: %{y:.7s} km<br>"
        "z: %{z:.7s} km<br>"
        "<extra></extra>"
    )

    # Attach customdata (UTC and speed) and a hovertemplate for Voyager
    voyager_custom = np.column_stack([utcs, speeds]) if len(utcs) == len(speeds) else None
    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='magenta', width=1),
                               name="Voyager 1",
                               customdata=voyager_custom,
                               hovertemplate=hovertemplate))
    v2_custom = np.column_stack([v2_utcs, v2_speeds]) if len(v2_utcs) == len(v2_speeds) else None
    fig.add_trace(go.Scatter3d(x=v2_x, y=v2_y, z=v2_z, mode='lines', line=dict(color='cyan', width=1),
                               name="Voyager 2",
                               customdata=v2_custom,
                               hovertemplate=hovertemplate))
    fig.add_trace(go.Scatter3d(x=v_x, y=v_y, z=v_z, mode='lines', line=dict(color='gold', width=1),
                               name="Venus"))
    fig.add_trace(go.Scatter3d(x=e_x, y=e_y, z=e_z, mode='lines', line=dict(color='blue', width=1),
                               name="Earth"))
    fig.add_trace(go.Scatter3d(x=j_x, y=j_y, z=j_z, mode='lines', line=dict(color='orange', width=1),
                               name="Jupiter"))
    fig.add_trace(go.Scatter3d(x=s_x, y=s_y, z=s_z, mode='lines', line=dict(color='brown', width=1),
                               name="Saturn"))
    fig.add_trace(go.Scatter3d(x=me_x, y=me_y, z=me_z, mode='lines', line=dict(color='gray', width=1),
                               name="Mercury"))
    fig.add_trace(go.Scatter3d(x=ma_x, y=ma_y, z=ma_z, mode='lines', line=dict(color='red', width=1),
                               name="Mars"))
    fig.add_trace(go.Scatter3d(x=u_x, y=u_y, z=u_z, mode='lines', line=dict(color='teal', width=1),
                               name="Uranus"))
    fig.add_trace(go.Scatter3d(x=n_x, y=n_y, z=n_z, mode='lines', line=dict(color='royalblue', width=1),
                               name="Neptune"))
    fig.add_trace(go.Scatter3d(x=p_x, y=p_y, z=p_z, mode='lines', line=dict(color='purple', width=1),
                               name="Pluto"))
    # Sun at origin
    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode='markers+text',
                               marker=dict(size=6, color='yellow'), text=['Sun'],
                               textposition='top center', name='Sun'))


    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(title="x (km) [Distance from Sun]", autorange="reversed"),
            yaxis=dict(title="y (km)", autorange="reversed"),
            zaxis_title="z (km)"),
        title="Voyagers (Saturn on X axis at Voyager 1's closest approach)",)

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "vger.html")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
