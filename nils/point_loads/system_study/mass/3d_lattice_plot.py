# cube: lines = constant coordinate values (one colour per coord), faces culled by view.
import numpy as np, itertools, matplotlib.pyplot as plt

size = 1.0
linewidth = 1.5
cols = ("#871a5b","#206095","#a8bd3a")  # x-const, y-const, z-const

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(projection="3d")

t = np.array([0, size])                  # parameter for the varying coordinate on a line
samples = np.linspace(0, size, 4)

def view_vec(ax):
    az, el = np.deg2rad(ax.azim), np.deg2rad(ax.elev)
    return np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])

v = view_vec(ax)

# faces where each constant-coordinate (axis) should appear:
faces = {
    0: [(1, 0.0), (2, size)],   # x-const → front(y=0) & top(z=size)
    1: [(2, size), (0, size)],  # y-const → top(z=size) & right(x=size)
    2: [(1, 0.0), (0, size)]    # z-const → front(y=0) & right(x=size)
}

for axis, col in enumerate(cols):            # axis = coordinate which is constant (x,y,z)
    for fixed_coord, fixed_val in faces[axis]:
        # remaining coord that varies along this face line
        samp_coord = ({0,1,2} - {axis, fixed_coord}).pop()
        normal = (-1 if np.isclose(fixed_val, 0.0) else 1) * np.eye(3)[fixed_coord]
        if np.dot(normal, v) <= 0:           # face not visible
            continue
        for s in samples:                    # each constant-coordinate sample
            pts = [None, None, None]
            pts[axis] = np.full_like(t, s)           # constant coordinate
            pts[fixed_coord] = np.full_like(t, fixed_val)
            pts[samp_coord] = t                      # vary across face
            ln, = ax.plot(pts[0], pts[1], pts[2], linewidth=linewidth, color=col)
            try:
                ln.set_solid_capstyle('round'); ln.set_solid_joinstyle('round')
            except Exception:
                pass

try: ax.set_box_aspect((1,1,1))
except Exception: pass
ax.set_axis_off()
plt.savefig('3d_lattice.svg', format='svg')
plt.tight_layout(); plt.show()
 