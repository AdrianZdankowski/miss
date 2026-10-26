import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Data (change here for a different task variant)
# ---------------------------------------------------------------------------
G = 10.0           # gravitational acceleration [m/s^2]
M = 1.0            # mass [kg]
N_BOUNCES = 5      # number of simulated impacts

P0 = np.array([-5.0, 0.0, 200.0])   # x0, y0, z0
V0 = np.array([2.0, 4.0, 0.0])      # v0x, v0y, v0z
K_LOSS = 0.8                        # k = Ekin(after) / Ekin(before)


def f(x, y):
    """Surface z = f(x, y)."""
    return x**2 + y**2

def grad_f(x, y):
    """Partial derivatives [df/dx, df/dy] of the surface."""
    return np.array([2 * x, 2 * y])

T_MAX = 100.0      # longest flight time searched for an impact [s]
DT_SCAN = 1e-3     # step used to bracket the impact time [s]
T_EPS = 1e-9       # roots t <= T_EPS are rejected (t = 0 is the previous impact)


def position(p0, v0, t):
    """Eq. (2): position of the body after time t."""
    t = np.asarray(t, dtype=float)
    x = p0[0] + v0[0] * t
    y = p0[1] + v0[1] * t
    z = p0[2] + v0[2] * t - 0.5 * G * t**2
    return x, y, z


def velocity(v0, t):
    """Eq. (1): velocity of the body after time t."""
    return np.array([v0[0], v0[1], v0[2] - G * t])


def bisection(g, a, b, tol=1e-14, max_iter=200):
    """Root of g in [a, b], where g(a) and g(b) have opposite signs."""
    ga = g(a)
    for _ in range(max_iter):
        c = 0.5 * (a + b)
        gc = g(c)
        if gc == 0 or (b - a) / 2 < tol:
            return c
        if np.sign(gc) == np.sign(ga):
            a, ga = c, gc
        else:
            b = c
    return 0.5 * (a + b)


def impact_time(p0, v0):
    """First positive root of eq. (3): z(t) - f(x(t), y(t)) = 0.

    Returns None when there is no collision within T_MAX.
    """
    def g(t):
        x, y, z = position(p0, v0, t)
        return z - f(x, y)

    t = np.arange(T_EPS, T_MAX + DT_SCAN, DT_SCAN)
    gt = g(t)
    idx = np.nonzero(np.sign(gt[:-1]) != np.sign(gt[1:]))[0]
    if idx.size == 0:
        return None
    i = idx[0]
    return bisection(g, t[i], t[i + 1])


def unit_normal(x, y):
    """Eq. (7): n = N / ||N||, N = [-df/dx, -df/dy, 1]."""
    fx, fy = grad_f(x, y)
    N = np.array([-fx, -fy, 1.0])
    return N / np.linalg.norm(N)


def reflect(v, n, k):
    """Eq. (8): v_after = sqrt(k) * (v_before - 2 (v_before . n) n)."""
    return np.sqrt(k) * (v - 2 * np.dot(v, n) * n)


def simulate(p0, v0, k, n_bounces=N_BOUNCES):
    p, v = np.array(p0, float), np.array(v0, float)
    rows, trajectories = [], []
    for _ in range(n_bounces):
        tI = impact_time(p, v)
        if tI is None:
            print("No further collision - simulation stopped.")
            break

        t = np.linspace(0.0, tI, 200)
        trajectories.append(np.column_stack(position(p, v, t)))

        pI = np.array(position(p, v, tI), dtype=float)
        v_before = velocity(v, tI)
        e_kin = 0.5 * M * np.dot(v_before, v_before)
        e_pot = M * G * pI[2]
        rows.append([*pI, tI, e_kin, e_pot, e_kin + e_pot])

        v = reflect(v_before, unit_normal(pI[0], pI[1]), k)
        p = pI
    return np.array(rows), trajectories

HEADER = ["xI", "yI", "zI", "dt", "Ekin", "Epot", "Etotal"]

def print_table(title, rows):
    print(f"\n{title}")
    print("".join(f"{h:>12}" for h in ["no."] + HEADER))
    for i, r in enumerate(rows, 1):
        print(f"{i:>12}" + "".join(f"{val:12.4f}" for val in r))

def save_csv(filename, rows):
    np.savetxt(filename, rows, delimiter=",", fmt="%.4f",
               header=",".join(HEADER), comments="")

def plot(title, filename, p0, trajectories, rows):
    pts = np.vstack(trajectories)
    margin = 2.0
    lim = max(np.abs(pts[:, :2]).max(), abs(p0[0]), abs(p0[1])) + margin
    xs = np.linspace(-lim, lim, 60)
    X, Y = np.meshgrid(xs, xs)
    Z = f(X, Y)
    Z = np.where(Z <= pts[:, 2].max() * 1.5, Z, np.nan)   # clip tall surface

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2, color="gray",
                      linewidth=0.4, alpha=0.6)
    for i, tr in enumerate(trajectories, 1):
        ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], linewidth=2,
                label=f"trajectory {i}")
    ax.scatter(rows[:, 0], rows[:, 1], rows[:, 2], color="red", s=25,
               depthshade=False, label="impact points")
    ax.scatter(*p0, color="black", s=25, depthshade=False,
               label="start point")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(title)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    return fig

def main():
    cases = [
        ("Results (without losses - k = 1)", 1.0, "no_losses"),
        (f"Results (with losses - k = {K_LOSS})", K_LOSS, "with_losses"),
    ]
    for title, k, tag in cases:
        rows, trajectories = simulate(P0, V0, k)
        print_table(title, rows)
        save_csv(f"results_{tag}.csv", rows)
        plot(f"Body bouncing off z = x^2 + y^2, k = {k}",
             f"trajectories_{tag}.png", P0, trajectories, rows)
    plt.show()

if __name__ == "__main__":
    main()
