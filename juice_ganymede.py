import spiceypy as sp
import numpy as np
import plotly.graph_objects as go
import sys
import os

# Plot trajectory for the final portion of the JUICE mission, starting just before the approach to Jupiter-Ganymede L2
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
# This is plotted as the JUICE-Ganymede relative position vector in the J2000 frame,
# in order to remove the effect of Ganymede's rotation beneath the probe.

ganymede_radius = 2634

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
    Sample the object and return state relative to Ganymede in the JUICE frame
    """
    targ = str(sc_id)
    ets = np.linspace(et_start, et_start + period, npts)
    x, y, z = [], [], []
    for et in ets:
        st_sc, _ = sp.spkezr(targ, et, "JUICE", "NONE", "GANYMEDE")
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

    # Plot only the portion starting just before the approach to Jupiter-Ganymede L2
    et_start = sp.utc2et("2034 OCT 28 00:00:00")

    # We want 256 samples per day; enough to render smooth curves
    total_pts = int(256.0 * (et_end - et_start) / 86400.0)

    print(f"Sampling mission from {sp.et2utc(et_start,'C',0)} to {sp.et2utc(et_end,'C',0)}")
    print(f"Total samples: {total_pts}")

    (x, y, z) = sample_object(sc_id, et_start, (et_end - et_start), total_pts)

    # Build plot
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='magenta', width=1),
                               name="JUICE probe"))
    # Ganymede sphere
    fig.add_trace(make_sphere(tuple(np.array([0.0, 0.0, 0.0])), ganymede_radius, 'gray', "Ganymede"))

    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(title="x (km) [Distance from Ganymede]"),
            yaxis=dict(title="y (km)"),
            zaxis_title="z (km)"),
        title="JUICE probe (Position relative to Ganymede)",)

    if "--static" in sys.argv:
        os.makedirs("docs", exist_ok=True)
        html_path = os.path.join("docs", "juice_ganymede.html")
        fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
        print(f"Wrote interactive HTML: {html_path}")
    else:
        fig.show()


if __name__ == "__main__":
    main()
