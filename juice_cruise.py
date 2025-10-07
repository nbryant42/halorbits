import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

# Plot planned trajectory for the cruise portion of the JUICE mission, ending at Jupiter orbit insertion (July 2031)
# as per the so-called Dataset v460 -- 150lb CReMA 5.1 Option A3.2 scenario, as of 20250918_001
#
# ```
# Scenario 2: modifications with respect to original Option A3 correspond
#    to the inclusion of the Ganymede phase with the following characteristics:
#   
#     - GCO-500 (N,P,Q)=(56,5,19) extended by 30 days to recover lost science at
#       the beginning of the phase due to the superior conjunction (165 days
#       from GCO500 start till last GCO5--/GCO200 manoeuvre).
#   
#     - Beta angle profile adjusted at the beginning/end of GCO500 to reduce
#       delta-v
#   
#     - GEO of 150 days with 12 hours period
#   
#     - GCO200 (N, P Q): (65, -1, 4).
# ```
#
# This is plotted as the Sun-referenced position vector in a custom frozen frame, "JUICE"

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
    Sample the object and return state relative to Sun in the JUICE frame
    """
    targ = str(sc_id)
    ets = np.linspace(et_start, et_start + period, npts)
    x, y, z = [], [], []
    for et in ets:
        st_sc, _ = sp.spkezr(targ, et, "JUICE", "NONE", "SUN")
        r_plot = st_sc[:3]
        x.append(r_plot[0]); y.append(r_plot[1]); z.append(r_plot[2])

    return (x, y, z)


def main():
    # Use the DE432 ephemeris (version used by JUICE) and the provided JUICE SPK
    bsp_name = "juice_crema_5_1_150lb_23_1_a3_2_v01.bsp"
    kernels = [
        "naif0012.tls",
        "pck00010.tpc",
        "de432s.bsp",
        "jup365_19900101_20500101.bsp",
        "juice.tk",
        bsp_name
    ]

    load_kernels(kernels)

    sc_id = -28

    # Get coverage window from the JUICE BSP
    cover = sp.spkcov(bsp_name, sc_id)
    et_start, et_end = sp.wnfetd(cover, 0)
    print("Coverage start (UTC):", sp.et2utc(et_start, 'C', 0))
    print("Coverage end   (UTC):", sp.et2utc(et_end, 'C', 0))

    # Plot until just after the first Jupiter swingby. After this point,
    # it makes more sense to look at relative motion; see juice_jupiter.py and juice_ganymede.py
    et_end = sp.utc2et("2031 AUG 01 00:00:00")

    # We want 1 sample per day; enough to render smooth curves
    total_pts = int((et_end - et_start) / 86400.0)

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Total samples: {total_pts}")

    (x, y, z) = sample_object(sc_id, et_start, (et_end - et_start), total_pts)
    (v_x, v_y, v_z) = sample_object("VENUS", et_start, (et_end - et_start), total_pts)
    (e_x, e_y, e_z) = sample_object("EARTH", et_start, (et_end - et_start), total_pts)
    (j_x, j_y, j_z) = sample_object("JUPITER", et_start, (et_end - et_start), total_pts)

    # Build plot
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='magenta', width=1),
                               name="JUICE probe"))
    fig.add_trace(go.Scatter3d(x=v_x, y=v_y, z=v_z, mode='lines', line=dict(color='gold', width=1),
                               name="Venus"))
    fig.add_trace(go.Scatter3d(x=e_x, y=e_y, z=e_z, mode='lines', line=dict(color='blue', width=1),
                               name="Earth"))
    fig.add_trace(go.Scatter3d(x=j_x, y=j_y, z=j_z, mode='lines', line=dict(color='orange', width=1),
                               name="Jupiter"))
    # Sun at origin
    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode='markers+text',
                               marker=dict(size=6, color='yellow'), text=['Sun'],
                               textposition='top center', name='Sun'))


    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(title="x (km) [Distance from Sun]"),
            yaxis=dict(title="y (km)"),
            zaxis_title="z (km)"),
        scene_camera=dict(eye=dict(x=3.0, y=3.0, z=3.0)),
        title="JUICE probe",)

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "juice_cruise.html")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
