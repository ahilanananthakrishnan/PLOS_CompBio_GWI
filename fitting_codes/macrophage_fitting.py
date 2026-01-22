import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# ============================================================
# DATA
# Each row = [tPB, M(0), M(7), M(30)]
# ============================================================
M0 = 220.0  # used only to normalize experimental data

mac_data = np.array([
    [0.0,   M0, 171.0,  75.0],
    [3.5,   M0, 200.0, 107.5],
    [7.0,   M0, 333.0, 140.0],
], dtype=float)

tPB_sets = mac_data[:, 0]
t_exp = np.array([0.0, 7.0, 30.0], dtype=float)

# normalized observations m = M/M0
m_data = mac_data[:, 1:] / M0

# ============================================================
# NORMALIZED MODEL with d1_M FIXED TO 0
# Params to fit: a_m, d2_M   (d1_M = 0.0, not calibrated)
# ============================================================
d1_M_fixed = 0.0

def m_of_t(t, tPB, a_m, d2_M):
    """
    Normalized macrophage level m(t)=M(t)/M0 with m(0)=1.

      dm/dt = a_m              for t < tPB   (since d1_M = 0)
      dm/dt = -d2_M * m        for t >= tPB

    Analytic solution:
      For t < tPB:
         m(t) = 1 + a_m * t
      For t >= tPB:
         m(t) = (1 + a_m * tPB) * exp(-d2_M * (t - tPB))
    """
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)

    pre = t < tPB
    post = ~pre

    out[pre] = 1.0 + a_m * t[pre]

    m_tPB = 1.0 + a_m * tPB
    out[post] = m_tPB * np.exp(-d2_M * (t[post] - tPB))

    return out

# ============================================================
# GLOBAL FIT (all tPB sets share the same a_m and d2_M)
# ============================================================
t_all = np.tile(t_exp, len(tPB_sets))         # [0,7,30, 0,7,30, ...]
tPB_all = np.repeat(tPB_sets, len(t_exp))     # [tPB,tPB,tPB, ...]
m_all = m_data.reshape(-1)                    # flattened normalized observations

def residuals(p):
    a_m, d2_M = p
    if (a_m < 0) or (d2_M < 0):
        return 1e3 * np.ones_like(m_all)
    m_pred = np.array(
        [m_of_t(ti, tPB_i, a_m, d2_M) for ti, tPB_i in zip(t_all, tPB_all)],
        dtype=float
    ).reshape(-1)
    return (m_pred - m_all)

# Initial guess in normalized space: [a_m, d2_M]
p0 = np.array([0.06, 0.05], dtype=float)

res = least_squares(residuals, p0, max_nfev=200000)
a_m_fit, d2_M_fit = res.x

# ============================================================
# PRINT 
# ============================================================
print("*** FITTED NORMALIZED MACROPHAGE PARAMETERS ***")
print(f"a_m   (normalized production) [1/day] = {a_m_fit:.7g}")
print(f"d1_M  (decay for t < tPB)      [1/day] = {d1_M_fixed:.7g}")
print(f"d2_M  (decay for t >= tPB)     [1/day] = {d2_M_fit:.7g}")

# ============================================================
# (normalized plot with error bars)
# ============================================================
tPB_control = 1e-3
tPB_values_plot = np.array([tPB_control, 3.0, 7.0, 10.0], dtype=float)
t_fine = np.linspace(0.0, 30.0, 300)

# Experimental points (only for control and tPB=7)
t_data = np.array([7.0, 30.0], dtype=float)

M_var_obs_0 = np.array([171.0, 75.0], dtype=float)
M_var_obs_7 = np.array([333.0, 140.0], dtype=float)
m_var_obs_0 = M_var_obs_0 / M0
m_var_obs_7 = M_var_obs_7 / M0

M_std_0 = np.array([12.1, 42.2], dtype=float)
M_std_7 = np.array([22.02, 79.15], dtype=float)
m_std_0 = M_std_0 / M0
m_std_7 = M_std_7 / M0

inferno = plt.get_cmap('inferno')
color_control = 'k'
color_tpb3 = inferno(0.25)
color_tpb7 = inferno(0.55)
color_tpb10 = inferno(0.85)

def color_for_tpb(tPB):
    if tPB <= tPB_control:
        return color_control
    if np.isclose(tPB, 3.0):
        return color_tpb3
    if np.isclose(tPB, 7.0):
        return color_tpb7
    if np.isclose(tPB, 10.0):
        return color_tpb10
    return inferno(0.5)

def set_axis_font_times(ax, tick_size=10):
    for item in (ax.title, ax.xaxis.label, ax.yaxis.label):
        item.set_fontname('Times New Roman')
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontname('Times New Roman')
    ax.tick_params(labelsize=tick_size)

def axis_style(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.grid(False)
    set_axis_font_times(ax, 19)
    for spine in ax.spines.values():
        spine.set_linewidth(1.44)

fig, ax = plt.subplots(1, 1, figsize=(6.0, 5.0))

# Model curves (normalized)
for tPB in tPB_values_plot:
    col = color_for_tpb(tPB)
    m_curve = m_of_t(t_fine, tPB, a_m_fit, d2_M_fit)
    ax.plot(t_fine, m_curve, '-', lw=3.6, color=col)

# Experimental error bars (normalized)
ax.errorbar(t_data, m_var_obs_0, yerr=m_std_0,
            fmt='o', ms=12, mfc=(0.6, 1.0, 0.6), mec='k', mew=1.2,
            ecolor='k', elinewidth=1.44)

ax.errorbar(t_data, m_var_obs_7, yerr=m_std_7,
            fmt='o', ms=12, mfc=(0.6, 1.0, 0.6), mec='k', mew=1.2,
            ecolor='k', elinewidth=1.44)

ax.set_xlabel("Day", fontsize=19)
ax.set_ylabel(r"Normalized CD40$^{+}$ expression $M$", fontsize=19)
ax.set_title("Macrophage level over time", fontsize=22)

ax.set_box_aspect(1)
axis_style(ax, (0, 32), (0, 2.0))

plt.tight_layout()
plt.show()
