
import numpy as np
import pandas as pd
from scipy.optimize import minimize, least_squares
import matplotlib.pyplot as plt

tPB_control = 1e-3  # used everywhere to represent control

# ============================================================
# INITIAL GUESSES + BOUNDS
# ============================================================

# Gene parameters p = [r_drop, x_ss, k_rec]
#   r_drop : linear drop rate during exposure (t < t_PB)
#   x_ss   : steady-state fold-change for t >= t_PB
#   k_rec  : exponential recovery rate after exposure
p0_ch = np.array([0.10, 1.65, 0.10])
lb_ch = np.array([-1.00, -5.00, -2.00])
ub_ch = np.array([ 1.00,  5.00,  2.00])

p0_nos = np.array([0.08, 0.95, 0.50])
lb_nos = np.array([-1.00, -5.00, -2.00])
ub_nos = np.array([ 1.00,  5.00,  2.00])

# Active-stress parameters (NO q):
#   Excitatory theta_E = [alpha_act_ex, k_ex, gamma_ch]
#   Inhibitory theta_I = [beta_act_in,  k_in, gamma_nos]
th0E = np.array([1.80, 0.020, 0.005])
lbE  = np.array([-10.0, -1.000, -5.000])
ubE  = np.array([ 10.0,  1.000,  5.000])

th0I = np.array([3.00, 0.080, 0.015])
lbI  = np.array([-10.0, -1.000, -5.000])
ubI  = np.array([ 10.0,  1.000,  5.000])

use_bounded_lbfgsb = True

# ============================================================
# TIME / PB DURATIONS
# ============================================================
tPB_values_plot = np.array([tPB_control, 3.0, 7.0, 10.0], dtype=float)  # control first
t_gene   = np.arange(0.0, 30.1, 0.1)
t_values = np.linspace(0.0, 30.0, 300)

# ============================================================
# COLORS (inferno-based)
# ============================================================
inferno = plt.get_cmap('inferno')
color_control = 'k'
color_tpb3    = inferno(0.25)
color_tpb7    = inferno(0.55)
color_tpb10   = inferno(0.85)

def color_for_tpb(tPB):
    """Control (tPB≈0) is always black; others are inferno."""
    if tPB <= tPB_control:
        return color_control
    if np.isclose(tPB, 3.0):
        return color_tpb3
    if np.isclose(tPB, 7.0):
        return color_tpb7
    if np.isclose(tPB, 10.0):
        return color_tpb10
    return inferno(0.5)

# ============================================================
# EXPERIMENTAL GENE DATA
# ============================================================

# ChAT
chat_day7 = np.array([0.18, 0.22, 0.14, 0.6, 0.76, 0.64, 0.22, 0.2, 0.07])
chat_day30 = np.array([
    2.46, 2.76, 2.7, 1.76, 2.02, 1.79, 1.0, 1.06, 0.93,
    0.891, 1.106, 1.192, 1.543, 1.534, 1.159
])
chat_t  = np.array([7.0, 30.0])
chat_mu = np.array([chat_day7.mean(), chat_day30.mean()])
chat_sd = np.array([chat_day7.std(ddof=1), chat_day30.std(ddof=1)])

# Nos1
nos1_day7 = np.array([0.683,0.689,0.677,0.304,0.368,0.352,0.39,0.619,0.554])
nos1_day30 = np.array([
    0.63,0.621,0.674,0.649,0.618,0.633,0.363,0.365,0.373,
    0.895,1.208,1.059,0.585,0.514,0.575
])
nos1_t  = np.array([7.0, 30.0])
nos1_mu = np.array([nos1_day7.mean(), nos1_day30.mean()])
nos1_sd = np.array([nos1_day7.std(ddof=1), nos1_day30.std(ddof=1)])

# ============================================================
# GENE MODEL: piecewise (linear drop then exponential recovery)
# ============================================================

def piecewise_gene(t, tPB, r_drop, x_ss, k_rec):
    """
    Fold-change x(t; t_PB):

      if t < t_PB:
        x = 1 - r_drop * t

      if t >= t_PB:
        x = x_ss - (x_ss - x(t_PB)) * exp(-k_rec * (t - t_PB))
    """
    t = np.asarray(t, dtype=float)
    x = np.empty_like(t)

    pre  = t < tPB
    post = ~pre

    x[pre] = 1.0 - r_drop * t[pre]

    x_tPB = 1.0 - r_drop * tPB
    x[post] = x_ss - (x_ss - x_tPB) * np.exp(-k_rec * (t[post] - tPB))
    return x

def ch(t, tPB, p_ch):
    """Excitatory gene fold-change ch(t; t_PB)."""
    return piecewise_gene(t, tPB, p_ch[0], p_ch[1], p_ch[2])

def nos(t, tPB, p_nos):
    """Inhibitory gene fold-change nos(t; t_PB)."""
    return piecewise_gene(t, tPB, p_nos[0], p_nos[1], p_nos[2])

# ============================================================
# ============================================================

tpb0 = tPB_control
t_targets = np.array([7.0, 30.0], dtype=float)
base_targets = np.array([1.0, 1.0], dtype=float)
w0 = 1.0

def chat_obj(p):
    c7 = ch(t_targets, 7.0, p)
    c0 = ch(t_targets, tpb0, p)
    return np.sum((c7 - chat_mu)**2) + w0 * np.sum((c0 - base_targets)**2)

def nos_obj(p):
    n7 = nos(t_targets, 7.0, p)
    n0 = nos(t_targets, tpb0, p)
    return np.sum((n7 - nos1_mu)**2) + w0 * np.sum((n0 - base_targets)**2)

def bounded_minimize(obj, p0, lb, ub):
    res = minimize(
        obj, p0,
        bounds=list(zip(lb, ub)),
        method='L-BFGS-B',
        options={'maxfun': 2e5, 'maxiter': 1e5}
    )
    return res.x

if use_bounded_lbfgsb:
    p_ch  = bounded_minimize(chat_obj, p0_ch,  lb_ch,  ub_ch)
    p_nos = bounded_minimize(nos_obj,  p0_nos, lb_nos, ub_nos)
else:
    p_ch  = minimize(chat_obj, p0_ch).x
    p_nos = minimize(nos_obj,  p0_nos).x

# ============================================================
# STRESS DATA
# ============================================================

data_arr = np.array([
    [ 1,0.2,1e-4,7,   0.9402,0.56922],
    [ 1,0.2,1e-4,30,  0.4361,0.17835],
    [ 1,0.2,7,   7,   0.1438,0.058478],
    [ 1,0.2,7,   30,  0.8913,0.50597],
    [-1,0.2,1e-4,7,  -1.0426,0.28057],
    [-1,0.2,1e-4,30, -0.3091,0.17608],
    [-1,0.2,7,   7,  -0.0496,0.03136],
    [-1,0.2,7,   30, -0.1557,0.058446]
], dtype=float)

data = pd.DataFrame(
    data_arr,
    columns=['Etype', 'stretch', 'tPB', 'teval', 'stress', 'std']
)

# ============================================================
# MACROPHAGE MODEL (normalized)
# ============================================================

a_M  = 14.1612
d2_M = 0.0368
M0   = 220.0

def M_norm(t, tPB):
    """
    """
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)

    pre  = t < tPB
    post = ~pre

    out[pre]  = (a_M * t[pre] + M0) / M0
    out[post] = ((a_M * tPB + M0) / M0) * np.exp(-d2_M * (t[post] - tPB))
    return out

def intM_norm(teval, tPB):
    """
    """
    teval = np.atleast_1d(teval).astype(float)
    out = np.zeros_like(teval)
    for i, ti in enumerate(teval):
        if ti > 0.0:
            grid = np.linspace(0.0, ti, 800)
            out[i] = np.trapezoid(M_norm(grid, tPB), grid)
    return out

# ============================================================
# ACTIVE STRESS (Step-3 matching form; no q)
# ============================================================

def phi_ex(teval, k_ex):
    """Excitatory phenomenological factor phi_ex(teval)."""
    teval = np.asarray(teval, dtype=float)
    return 1.0 / (1.0 + np.exp(k_ex * teval))

def phi_in(teval, k_in):
    """Inhibitory phenomenological factor phi_in(teval)."""
    teval = np.asarray(teval, dtype=float)
    return 1.0 - 1.0 / (1.0 + np.exp(k_in * teval))

def chi_ch(teval, tPB, k_ex, gamma_ch, t_c):
    """
    Excitatory desensitization envelope excluding ch:
      chi_ch = phi_ex(teval) * exp(-gamma_ch * intM_norm(teval; t_PB) / t_c)
    """
    return phi_ex(teval, k_ex) * np.exp(-gamma_ch * intM_norm(teval, tPB) / t_c)

def chi_nos(teval, tPB, k_in, gamma_nos, t_c):
    """
    Inhibitory desensitization envelope excluding nos:
      chi_nos = phi_in(teval) * exp(-gamma_nos * intM_norm(teval; t_PB) / t_c)
    """
    return phi_in(teval, k_in) * np.exp(-gamma_nos * intM_norm(teval, tPB) / t_c)

def sigE(theta_E, tPB, teval, t_c=30.0):
    """
    Excitatory active stress sigma_ex^{act} (kPa):
      theta_E = [alpha_act_ex, k_ex, gamma_ch]
      eta_ex  = chi_ch * ch
      sigma_ex^{act} = alpha_act_ex * eta_ex
    """
    alpha_act_ex, k_ex, gamma_ch = theta_E
    teval = np.asarray(teval, dtype=float)

    ch_val = np.ones_like(teval) if tPB <= tPB_control else ch(teval, tPB, p_ch)
    eta_ex = chi_ch(teval, tPB, k_ex, gamma_ch, t_c) * ch_val
    return alpha_act_ex * eta_ex

def sigI(theta_I, tPB, teval, t_c=30.0):
    """
    Inhibitory active stress sigma_in^{act} (kPa):
      theta_I = [beta_act_in, k_in, gamma_nos]
      eta_in  = chi_nos * nos
      sigma_in^{act} = - beta_act_in * eta_in
    """
    beta_act_in, k_in, gamma_nos = theta_I
    teval = np.asarray(teval, dtype=float)

    nos_val = np.ones_like(teval) if tPB <= tPB_control else nos(teval, tPB, p_nos)
    eta_in  = chi_nos(teval, tPB, k_in, gamma_nos, t_c) * nos_val
    return -beta_act_in * eta_in

# ============================================================
# FIT theta_E and theta_I on t_PB in {control, 7}
# ============================================================

isCal = (data['tPB'] <= tPB_control) | (np.isclose(data['tPB'], 7.0))
Ecal = data[(data['Etype'] ==  1) & isCal].reset_index(drop=True)
Ical = data[(data['Etype'] == -1) & isCal].reset_index(drop=True)

def predE(tbl, theta):
    return np.array([float(sigE(theta, row.tPB, row.teval)) for row in tbl.itertuples(index=False)])

def predI(tbl, theta):
    return np.array([float(sigI(theta, row.tPB, row.teval)) for row in tbl.itertuples(index=False)])

def residE(theta):
    return Ecal['stress'].values - predE(Ecal, theta)

def residI(theta):
    return Ical['stress'].values - predI(Ical, theta)

theta_E = least_squares(residE, th0E, bounds=(lbE, ubE), max_nfev=200000, verbose=0).x
theta_I = least_squares(residI, th0I, bounds=(lbI, ubI), max_nfev=200000, verbose=0).x

# ============================================================
# REPORT
# ============================================================

print("\n======== FITTED GENE PARAMETERS ========")
print("ChAT  p_ch  = [r_drop, ch_ss, k_ch]  =", p_ch)
print("Nos1  p_nos = [r_drop, nos_ss, k_nos] =", p_nos)

print("\n======== FITTED ACTIVE-STRESS PARAMETERS ========")
print("theta_E = [alpha_act_ex, k_ex, gamma_ch]   =", theta_E)
print("theta_I = [beta_act_in,  k_in, gamma_nos]  =", theta_I)

# ============================================================
# PLOTTING UTILITIES
# ============================================================

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
    set_axis_font_times(ax, 10)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

# ============================================================
# PLOTS (NO chi FIGURE)
# ============================================================

# -------------------- FIGURE 1: Genes --------------------
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

# ChAT curves
for tPB in tPB_values_plot:
    col = color_for_tpb(tPB)
    if tPB <= tPB_control:
        ax1.plot(t_gene, np.ones_like(t_gene), '-', lw=2.2, color=col)
    else:
        ax1.plot(t_gene, ch(t_gene, tPB, p_ch), '-', lw=2.2, color=col)

ax1.errorbar(chat_t, chat_mu, yerr=chat_sd,
             fmt='o', ms=7, mfc='w', mec=color_tpb7, mew=1.8,
             ecolor=color_tpb7, elinewidth=1.2)
ax1.set_xlabel("Day", fontsize=11)
ax1.set_ylabel(r"ChAT fold change $ch$", fontsize=11)
ax1.set_title("A   Excitatory neuron gene expression", loc='left', fontsize=11)
axis_style(ax1, (0, 32), (0, 2.5))

# Nos1 curves
for tPB in tPB_values_plot:
    col = color_for_tpb(tPB)
    if tPB <= tPB_control:
        ax2.plot(t_gene, np.ones_like(t_gene), '-', lw=2.2, color=col)
    else:
        ax2.plot(t_gene, nos(t_gene, tPB, p_nos), '-', lw=2.2, color=col)

ax2.errorbar(nos1_t, nos1_mu, yerr=nos1_sd,
             fmt='o', ms=7, mfc='w', mec=color_tpb7, mew=1.8,
             ecolor=color_tpb7, elinewidth=1.2)
ax2.set_xlabel("Day", fontsize=11)
ax2.set_ylabel(r"Nos1 fold change $nos$", fontsize=11)
ax2.set_title("B   Inhibitory neuron gene expression", loc='left', fontsize=11)
axis_style(ax2, (0, 32), (0, 2.5))

fig1.tight_layout()

# -------------------- FIGURE 2: Active stress --------------------
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(7.0, 3.5))

# Excitatory curves
for tPB in tPB_values_plot:
    col = color_for_tpb(tPB)
    ax3.plot(t_values, sigE(theta_E, tPB, t_values), '-', lw=2.2, color=col)

maskE0 = (data['Etype'] == 1) & (data['tPB'] <= tPB_control)
maskE7 = (data['Etype'] == 1) & np.isclose(data['tPB'], 7.0)

if maskE0.any():
    ax3.errorbar(data.loc[maskE0, 'teval'], data.loc[maskE0, 'stress'],
                 yerr=data.loc[maskE0, 'std'],
                 fmt='o', ms=7, mfc='w', mec='k', mew=1.8,
                 ecolor='k', elinewidth=1.2)
if maskE7.any():
    ax3.errorbar(data.loc[maskE7, 'teval'], data.loc[maskE7, 'stress'],
                 yerr=data.loc[maskE7, 'std'],
                 fmt='o', ms=7, mfc='w', mec=color_tpb7, mew=1.8,
                 ecolor=color_tpb7, elinewidth=1.2)

ax3.axhline(0.0, linestyle='--', color='k', lw=1.2)
ax3.set_xlabel("Day", fontsize=11)
ax3.set_ylabel(r"Excitatory active stress $\sigma^{act}_{ex}$ (kPa)", fontsize=11)
ax3.set_title("C   Excitatory active stress adaptation", loc='left', fontsize=11)
axis_style(ax3, (0, 32), (-0.5, 2.0))

# Inhibitory curves
for tPB in tPB_values_plot:
    col = color_for_tpb(tPB)
    ax4.plot(t_values, sigI(theta_I, tPB, t_values), '-', lw=2.2, color=col)

maskI0 = (data['Etype'] == -1) & (data['tPB'] <= tPB_control)
maskI7 = (data['Etype'] == -1) & np.isclose(data['tPB'], 7.0)

if maskI0.any():
    ax4.errorbar(data.loc[maskI0, 'teval'], data.loc[maskI0, 'stress'],
                 yerr=data.loc[maskI0, 'std'],
                 fmt='o', ms=7, mfc='w', mec='k', mew=1.8,
                 ecolor='k', elinewidth=1.2)
if maskI7.any():
    ax4.errorbar(data.loc[maskI7, 'teval'], data.loc[maskI7, 'stress'],
                 yerr=data.loc[maskI7, 'std'],
                 fmt='o', ms=7, mfc='w', mec=color_tpb7, mew=1.8,
                 ecolor=color_tpb7, elinewidth=1.2)

ax4.axhline(0.0, linestyle='--', color='k', lw=1.2)
ax4.set_xlabel("Day", fontsize=11)
ax4.set_ylabel(r"Inhibitory active stress $\sigma^{act}_{in}$ (kPa)", fontsize=11)
ax4.set_title("D   Inhibitory active stress adaptation", loc='left', fontsize=11)
axis_style(ax4, (0, 32), (-2.0, 0.5))

fig2.tight_layout()

plt.show()

