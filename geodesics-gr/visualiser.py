import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation

# ---------------- settings ----------------
MARKER_SIZE = 2.5      # <- red ball size (try 1.5 - 4)
FIGSIZE     = (16, 12) # bigger figure = smaller-looking balls relative to the frame
DPI         = 120      # 16x12 in * 120 = 1920x1440 video
FPS         = 120
SKIP        = 10
ROTATE      = False     # slow camera orbit
LIMIT       = None     # e.g. 12 to force a fixed view; None = automatic
# ------------------------------------------

data = pd.read_csv("trajectory.csv", sep=";")

r, theta, phi = (data[c].to_numpy() for c in ("r", "theta", "phi"))
data["x"] = r * np.sin(theta) * np.cos(phi)
data["y"] = r * np.sin(theta) * np.sin(phi)
data["z"] = r * np.cos(theta)

# Precompute coordinates per frame (much faster than filtering every frame)
frames_xyz = {f: g[["x", "y", "z"]].to_numpy() for f, g in data.groupby("frame")}
frames = np.array(sorted(frames_xyz))[::SKIP]

# Scale: ignore a few far-away outliers
M = 1.0
r_h = 2 * M
limit = LIMIT if LIMIT is not None else np.percentile(r, 99) * 1.1
limit = max(limit, 2 * r_h)

fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
ax = fig.add_subplot(projection="3d")
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)   # fill the whole frame

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)
ax.set_box_aspect((1, 1, 1), zoom=1.25)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

# Event horizon
u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, np.pi, 40)
ax.plot_surface(
    r_h * np.outer(np.cos(u), np.sin(v)),
    r_h * np.outer(np.sin(u), np.sin(v)),
    r_h * np.outer(np.ones_like(u), np.cos(v)),
    color="black", shade=False, linewidth=0, antialiased=False,
)

# Particles: small, no marker edge
particles, = ax.plot(
    [], [], [],
    linestyle="", marker="o", color="red",
    markersize=MARKER_SIZE, markeredgewidth=0,
)

def animate(i):
    xyz = frames_xyz[frames[i]]
    particles.set_data(xyz[:, 0], xyz[:, 1])
    particles.set_3d_properties(xyz[:, 2])
    if ROTATE:
        ax.view_init(elev=25, azim=30 + 0.15 * i)
    return particles,

ani = FuncAnimation(
    fig, animate,
    frames=len(frames),
    interval=1000 / FPS,
    blit=False,          # blit doesn't work reliably with 3D axes
)

ani.save("geodesics.mp4", writer="ffmpeg", fps=FPS, dpi=DPI, bitrate=8000)
plt.show()