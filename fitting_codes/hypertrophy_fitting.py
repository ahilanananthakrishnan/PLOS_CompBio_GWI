"""Hypertrophy fitting 
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# ============================================================
# EXPERIMENTAL TARGETS (ONLY VALUES USED IN FITTING)
# ============================================================

# --- IL-6 data (for t_PB = 7 days) ---
tPB_IL6 = 7.0
t_IL6_targets = np.array([7.0, 30.0], dtype=float)      # days
I_IL6_targets = np.array([7.3, 1.578], dtype=float)     # fold-change

# --- vartheta_max data (circular smooth muscle) ---
tPB_theta_max = np.array([0.0, 7.0], dtype=float)
theta_max_val = np.array([6.592, 7.029], dtype=float)

# --- Hypertrophy trajectories (growth multiplier vartheta^h) ---
t_growth = np.array([0.0, 7.0, 30.0], dtype=float)
theta_const_obs = np.array([1.0, 1.553, 6.592], dtype=float)  # healthy / control
theta_const_std = np.array([0.0, 0.429, 1.521], dtype=float)

theta_var_obs = np.array([1.0, 3.447, 7.029], dtype=float)    # PB, t_PB=7
theta_var_std = np.array([0.0, 0.811, 1.828], dtype=float)    # PB, t_PB=7


# ============================================================
# UTILITIES
# ============================================================

def set_axis_font_times(ax, tick_size=12):
    """Consistent manuscript-style axes."""
    for item in (ax.title, ax.xaxis.label, ax.yaxis.label):
        item.set_fontname('Times New Roman')
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontname('Times New Roman')
    ax.tick_params(labelsize=tick_size)


# ============================================================
# CYTOKINE MODEL: I(t)
# ============================================================

def il6_piecewise(t, tPB, a_cyto, d1_cyto, d2_cyto):
    """
    Cytokine fold-change model I(t) (piecewise, as used in the manuscript calibration).

    During PB exposure (t < tPB):
        dI/dt = a_cyto - d1_cyto * I
    After exposure (t >= tPB):
        dI/dt = - d2_cyto * I

    Initial condition:
        I(0) = 1
    """
    t = np.atleast_1d(t)
    I = np.zeros_like(t)

    # Control: no exposure
    if tPB <= 1e-3:
        return np.ones_like(t)

    # Special case: d1_cyto ~ 0 -> linear rise during exposure
    if abs(d1_cyto) < 1e-8:
        mask_rise = t < tPB
        mask_decay = ~mask_rise
        I[mask_rise] = 1.0 + a_cyto * t[mask_rise]
        I_tPB = 1.0 + a_cyto * tPB
        I[mask_decay] = I_tPB * np.exp(-d2_cyto * (t[mask_decay] - tPB))
    else:
        mask_rise = t < tPB
        mask_decay = ~mask_rise
        I[mask_rise] = a_cyto / d1_cyto + (1.0 - a_cyto / d1_cyto) * np.exp(-d1_cyto * t[mask_rise])
        I_tPB = a_cyto / d1_cyto + (1.0 - a_cyto / d1_cyto) * np.exp(-d1_cyto * tPB)
        I[mask_decay] = I_tPB * np.exp(-d2_cyto * (t[mask_decay] - tPB))

    return I


def fit_il6_parameters():
    """
    Fit a_cyto and d2_cyto using IL-6 data at (tPB=7, t=7 and 30),
    with d1_cyto fixed to 0 (as in your current setup).
    """
    tPB = tPB_IL6
    t_data = t_IL6_targets
    I_data = I_IL6_targets

    def residuals(params):
        a_cyto, d2_cyto = params
        d1_cyto = 0.0
        I_pred = il6_piecewise(t_data, tPB, a_cyto, d1_cyto, d2_cyto)
        return I_pred - I_data

    # Initial guesses
    a0 = (I_data[0] - 1.0) / t_data[0]
    d20 = 0.05

    res = least_squares(
        residuals,
        x0=[a0, d20],
        bounds=([0.0, 0.0], [10.0, 1.0])
    )

    a_fit, d2_fit = res.x
    d1_fit = 0.0

    print("\n*** IL-6 fit (t_PB = 7, targets at t = 7 & 30) ***")
    print(f"a_cyto  = {a_fit:.4f}")
    print(f"d1_cyto = {d1_fit:.4f} (fixed)")
    print(f"d2_cyto = {d2_fit:.4f}")
    print(f"IL-6 residual norm = {np.linalg.norm(res.fun):.4e}")

    return a_fit, d1_fit, d2_fit


# ============================================================
# vartheta_max(tPB) from time-integral of I(t) over exposure window
# ============================================================

def fit_theta_max(a_cyto, d1_cyto, d2_cyto):
    """
    Fit vartheta_max_0 and b_hyper from two points (tPB=0 and 7).

    For a given tPB, we compute:
        int_I = integral_0^{tPB} I(t) dt

    Then assume a linear model:
        vartheta_max(tPB) = vartheta_max_0 + b_hyper * int_I
    """
    integrals = []
    for tPB in tPB_theta_max:
        if tPB == 0.0:
            integrals.append(0.0)
        else:
            val, _ = quad(lambda τ: il6_piecewise(τ, tPB, a_cyto, d1_cyto, d2_cyto)[0], 0.0, tPB)
            integrals.append(val)
    integrals = np.array(integrals, dtype=float)

    # Fit slope/intercept (b_hyper, vartheta_max_0)
    b_hyper, theta0 = np.polyfit(integrals, theta_max_val, 1)

    print("\n*** vartheta_max fit (t_PB = 0, 7) ***")
    print(f"vartheta_max_0 = {theta0:.4f}")
    print(f"b_hyper        = {b_hyper:.4f}")

    return theta0, b_hyper


def theta_max_of_tPB(tPB, a_cyto, d1_cyto, d2_cyto, theta0, b_hyper):
    """
    Evaluate vartheta_max(tPB) using the linear model in the manuscript:
        vartheta_max(tPB) = vartheta_max_0 + b_hyper * integral_0^{tPB} I(t) dt
    """
    if tPB <= 1e-3:
        return theta0
    integral, _ = quad(lambda τ: il6_piecewise(τ, tPB, a_cyto, d1_cyto, d2_cyto)[0], 0.0, tPB)
    return theta0 + b_hyper * integral


# ============================================================
# HILL-TYPE HYPERTROPHY MODEL: vartheta^h(t)
# ============================================================

def theta_trajectory(params, tPB, theta_init, theta_max, a_cyto, d1_cyto, d2_cyto, t_eval):
    """
    Hypertrophy / growth multiplier evolution (Hill-type form):

        d(vartheta^h)/dt = (c_hyper * I(t)^n / (K_hyper^n + I(t)^n)) * (1 - vartheta^h / vartheta_max)

    where vartheta_max depends on tPB via theta_max_of_tPB().
    """
    c_hyper, K_hyper, n_hyper = params

    def rhs(t, theta):
        I_val = il6_piecewise(t, tPB, a_cyto, d1_cyto, d2_cyto)[0]
        return (c_hyper * I_val**n_hyper / (K_hyper**n_hyper + I_val**n_hyper)) * (1.0 - theta / theta_max)

    sol = solve_ivp(
        rhs,
        (t_eval[0], t_eval[-1]),
        [theta_init],
        t_eval=t_eval,
        max_step=0.1
    )
    return sol.y[0]


def fit_hill_parameters(a_cyto, d1_cyto, d2_cyto, theta0, b_hyper):
    """
    Fit (c_hyper, K_hyper, n_hyper) by matching both:
      - control hypertrophy trajectory (tPB ~ 0)
      - PB hypertrophy trajectory (tPB = 7)

    """
    tPB_const = 1e-3  # control (approx 0)
    tPB_var   = 7.0   # PB

    theta_max_const = theta_max_of_tPB(tPB_const, a_cyto, d1_cyto, d2_cyto, theta0, b_hyper)
    theta_max_var   = theta_max_of_tPB(tPB_var,   a_cyto, d1_cyto, d2_cyto, theta0, b_hyper)

    def residuals(p):
        th_const = theta_trajectory(
            p, tPB_const, theta_const_obs[0],
            theta_max_const, a_cyto, d1_cyto, d2_cyto, t_growth
        )
        th_var = theta_trajectory(
            p, tPB_var, theta_var_obs[0],
            theta_max_var, a_cyto, d1_cyto, d2_cyto, t_growth
        )
        return np.concatenate([th_const - theta_const_obs,
                               th_var   - theta_var_obs])

    res = least_squares(
        residuals,
        x0=[1.0, 1.0, 1.0],
        bounds=([0.0, 0.01, 0.1], [10.0, 10.0, 10.0])
    )

    c_fit, K_fit, n_fit = res.x

    print("\n*** Hill-type hypertrophy fit ***")
    print(f"c_hyper = {c_fit:.4f}, K_hyper = {K_fit:.4f}, n_hyper = {n_fit:.4f}")

    return (c_fit, K_fit, n_fit)


# ============================================================
# MAIN: BUILD FIGURE WITH PANELS A & B ONLY
# ============================================================

if __name__ == "__main__":
    # 1) Fit cytokine parameters (a_cyto, d1_cyto, d2_cyto)
    a_fit, d1_fit, d2_fit = fit_il6_parameters()

    # 2) Fit vartheta_max model (vartheta_max_0, b_hyper)
    theta0_fit, b_fit = fit_theta_max(a_fit, d1_fit, d2_fit)

    # 3) Fit hypertrophy kinetics (c_hyper, K_hyper, n_hyper)
    hill_params = fit_hill_parameters(a_fit, d1_fit, d2_fit, theta0_fit, b_fit)

    # Time grid for plotting (days)
    t_fine = np.linspace(0.0, 30.0, 400)

    # Inferno colors (continuous)
    inferno = plt.get_cmap('inferno')
    color_control = 'k'
    color_tpb3    = inferno(0.25)
    color_tpb7    = inferno(0.55)
    color_tpb10   = inferno(0.85)

    colors_IL6 = {
        3.0:  color_tpb3,
        7.0:  color_tpb7,
        10.0: color_tpb10
    }

    colors_growth = {
        1e-3: color_control,
        3.0:  color_tpb3,
        7.0:  color_tpb7,
        10.0: color_tpb10
    }

    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.5))

    # ---------------- A: Cytokine panel (I) ----------------
    # Control (no exposure): I(t)=1 baseline
    ax1.plot([0, 30], [1, 1], color_control, lw=2.2, label="Control")

    # t_PB = 3, 7, 10 curves (same fitted cytokine parameters, different exposure durations)
    tPB_list = [3.0, 7.0, 10.0]
    labels_IL6 = {
        3.0:  r"$t_{\mathrm{PB}} = 3$",
        7.0:  r"$t_{\mathrm{PB}} = 7$",
        10.0: r"$t_{\mathrm{PB}} = 10$"
    }
    for tPB in tPB_list:
        col = colors_IL6[tPB]
        I_curve = il6_piecewise(t_fine, tPB, a_fit, d1_fit, d2_fit)
        ax1.plot(t_fine, I_curve, '-', lw=2.2, color=col, label=labels_IL6[tPB])

    # Experimental IL-6 points (for t_PB = 7 only)
    ax1.plot(t_IL6_targets, I_IL6_targets, 'o', ms=7, mfc='w', mec=color_tpb7, mew=1.8, label="Experiment")

    ax1.set_xlabel("Day", fontsize=11)
    ax1.set_ylabel(r"Cytokine fold change $I$", fontsize=11)
    ax1.set_title("A   Cytokine level over time", loc='left', fontsize=11)
    ax1.set_xlim(0, 30)
    ax1.set_ylim(0, 10)
    ax1.grid(False)
    set_axis_font_times(ax1, 10)
    ax1.legend(fontsize=8, loc="upper right", frameon=False)

    # ---------------- B: Hypertrophy panel (vartheta^h) ----------------
    # Plot vartheta^h(t) for control and t_PB=3,7,10
    tPB_growth_list = [1e-3, 3.0, 7.0, 10.0]  # control, 3, 7, 10
    for tPB in tPB_growth_list:
        theta_max = theta_max_of_tPB(tPB, a_fit, d1_fit, d2_fit, theta0_fit, b_fit)
        theta_curve = theta_trajectory(
            hill_params, tPB, 1.0, theta_max, a_fit, d1_fit, d2_fit, t_fine
        )
        ax2.plot(t_fine, theta_curve, '-', lw=2.2, color=colors_growth[tPB])

    # Experimental points with error bars (healthy vs PB)
    ax2.errorbar(
        t_growth, theta_const_obs, yerr=theta_const_std,
        fmt='o', ms=6, color=color_control, ecolor=color_control,
        elinewidth=1.2, capsize=3, label="Control (exp.)"
    )
    ax2.errorbar(
        t_growth, theta_var_obs, yerr=theta_var_std,
        fmt='o', ms=6, color=color_tpb7, ecolor=color_tpb7,
        elinewidth=1.2, capsize=3, label=r"$t_{\mathrm{PB}}=7$ (exp.)"
    )

    ax2.set_xlabel("Day", fontsize=11)
    ax2.set_ylabel(r"Growth multiplier $\vartheta^{\mathrm{h}}$", fontsize=11)
    ax2.set_title("B   CSM hypertrophy evolution over time", loc='left', fontsize=11)
    ax2.set_xlim(0, 30)
    ax2.set_ylim(0, 10)
    ax2.grid(False)
    set_axis_font_times(ax2, 10)
    ax2.legend(fontsize=8, loc="upper left", frameon=False)

    fig.tight_layout()
    plt.show()
