import spiceypy as sp
import numpy as np

# Load kernels
sp.kclear()
sp.reset()
sp.furnsh("naif0012.tls")
sp.furnsh("pck00010.tpc")
sp.furnsh("de432s.bsp")

et = sp.utc2et("2025-01-01T00:00:00")

# 1. Get state vectors in J2000 relative to EMB
st_moon_j2000, _ = sp.spkezr("MOON", et, "J2000", "NONE", "EMB")
pos_moon_j2000 = st_moon_j2000[:3]
vel_moon_j2000 = st_moon_j2000[3:]

# 2. Build the EMBR basis vectors in J2000
# X-axis: Normalized position vector (radial)
x_axis = sp.vhat(pos_moon_j2000)

# Z-axis: Normalized orbital normal (angular momentum vector)
z_candidate = sp.vcrss(pos_moon_j2000, vel_moon_j2000)
z_axis = sp.vhat(z_candidate)

# Y-axis: Cross product of Z and X to complete the right-handed system
y_axis = sp.vcrss(z_axis, x_axis)

# 3. Create the manual J2000 -> EMBR rotation matrix
# The rows of this matrix are the EMBR basis vectors in J2000 coordinates.
xfm_manual = np.array([x_axis, y_axis, z_axis])

# 4. Transform the moon's state vector from J2000 to the manual EMBR frame
# We need to apply the rotation matrix to the position and velocity parts separately
r_embr_manual = xfm_manual.dot(pos_moon_j2000)
v_embr_manual = xfm_manual.dot(vel_moon_j2000)

# Calculate and print results
R_embr_manual = np.linalg.norm(r_embr_manual)
ang_embr_x_manual = np.degrees(np.arccos(r_embr_manual[0] / R_embr_manual))
ang_embr_y_manual = np.degrees(np.arccos(r_embr_manual[1] / R_embr_manual))
ang_embr_z_manual = np.degrees(np.arccos(r_embr_manual[2] / R_embr_manual))

print("Manual EMBR basis vectors in J2000 (rows of the rotation matrix):")
print(xfm_manual)

print("\nMoon in manual EMBR (km):", r_embr_manual)
print("Distance (km):", R_embr_manual)
print("Angle to X (deg):", ang_embr_x_manual)
print("Angle to Y (deg):", ang_embr_y_manual)
print("Angle to Z (deg):", ang_embr_z_manual)
