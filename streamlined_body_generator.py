#!/usr/bin/env python3
"""
Exact-match implementation of your generate_streamlined_body_geometry
plus a finder for the beta_tail threshold that uses the exact same
grid / argmax logic as your code (no peak refinement).

Usage: run the script; it prints the threshold interval (step=0.01 deg by default)
and shows the geometry just below/above the threshold using preserve_global_max=True.
"""
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def generate_streamlined_body_geometry(
    R_le_over_c,
    beta_tail,
    Psi_zeta_max,
    zeta_max,
    zeta_te,
    dimensional_known_dict,
    preserve_global_max=False
):
    # Kulfan class function exponents for round-nose airfoil
    N_1, N_2 = 0.5, 1.0

    # shape-function endpoints (Kulfan Eqs. (4) and (5))
    S_le = np.sqrt(2.0 * R_le_over_c)
    S_te_orig = np.tan(np.deg2rad(beta_tail)) + zeta_te

    if not (0.0 < Psi_zeta_max < 1.0):
        raise ValueError("Psi_zeta_max must be between 0 and 1 (exclusive).")

    # class function C(Psi) and derivative
    def C(Psi):
        return (Psi**N_1) * ((1.0 - Psi)**N_2)
    def dC_dPsi(Psi):
        c = C(Psi)
        return c * (N_1 / Psi - N_2 / (1.0 - Psi))

    # Compute S at the specified Psi_zeta_max so zeta(Psi_zeta_max) == zeta_max
    C_at = C(Psi_zeta_max)
    S_at = (zeta_max - Psi_zeta_max * zeta_te) / C_at

    # Enforce stationary condition: d/dPsi [C*S + Psi*zeta_te] = 0 at Psi_zeta_max
    dC_at = dC_dPsi(Psi_zeta_max)
    S_prime_target = (- dC_at * S_at - zeta_te) / C_at

    # Helper: compute quadratic coefficients a,b,c passing through (x0,y0) with slope y0p and through (x1,y1)
    def quadratic_coeffs(x0, y0, y0p, x1, y1):
        dx = x1 - x0
        a = (y1 - y0 - y0p * dx) / (dx**2)
        b = y0p - 2.0 * a * x0
        c = y0 + a * x0**2 - y0p * x0
        return a, b, c

    # Build method that returns the zeta(Psi) for a given S_te
    def build_zeta_for_S_te(S_te, npts=201):
        a_f, b_f, c_f = quadratic_coeffs(Psi_zeta_max, S_at, S_prime_target, 0.0, S_le)
        a_a, b_a, c_a = quadratic_coeffs(Psi_zeta_max, S_at, S_prime_target, 1.0, S_te)
        Psi = np.linspace(0.0, 1.0, npts)
        zeta = np.empty_like(Psi)
        for i, p in enumerate(Psi):
            if p >= Psi_zeta_max:
                S_val = a_a*p**2 + b_a*p + c_a
            else:
                S_val = a_f*p**2 + b_f*p + c_f
            zeta[i] = C(p) * S_val + p * zeta_te
        return Psi, zeta

    # initial zeta with original S_te
    Psi, zeta = build_zeta_for_S_te(S_te_orig)
    idx_global = int(np.argmax(zeta))
    Psi_global = Psi[idx_global]    # NOTE: exact original behavior (no parabola refine)

    S_te_used = S_te_orig
    adjusted = False

    # if requested: force the global maximum to be at Psi_zeta_max by scaling S_te
    if preserve_global_max and abs(Psi_global - Psi_zeta_max) > 1e-8:
        # grid search on alpha in [0,1] where S_te = S_at + alpha*(S_te_orig - S_at)
        alphas = np.linspace(0.0, 1.0, 101)
        chosen_alpha = None
        for a in alphas:
            S_te_try = S_at + a * (S_te_orig - S_at)
            Psi_try, zeta_try = build_zeta_for_S_te(S_te_try)
            if Psi_try[int(np.argmax(zeta_try))] <= Psi_zeta_max + 1e-9:
                chosen_alpha = a
                break
        if chosen_alpha is None:
            S_te_used = S_at
        else:
            # refine alpha in a small interval around chosen_alpha
            low = max(0.0, chosen_alpha - 1.0/100)
            high = chosen_alpha
            for _ in range(20):
                a_mid = 0.5*(low + high)
                S_te_mid = S_at + a_mid*(S_te_orig - S_at)
                Psi_mid, zeta_mid = build_zeta_for_S_te(S_te_mid)
                if Psi_mid[int(np.argmax(zeta_mid))] <= Psi_zeta_max + 1e-10:
                    high = a_mid
                else:
                    low = a_mid
            S_te_used = S_at + high*(S_te_orig - S_at)
        # rebuild final zeta with chosen S_te_used
        Psi, zeta = build_zeta_for_S_te(S_te_used)
        adjusted = True

    # volume/length handling (keeps the same approach as before)
    l_over_d = 1.0 / (Psi_zeta_max * 2.0)
    Vol = None
    def volume_for_length(L):
        x = Psi * L
        y = zeta * L
        return np.pi * simpson(y**2, x=x)

    if list(dimensional_known_dict.keys())[0] == 'length':
        L = list(dimensional_known_dict.values())[0]
    else:
        target_vol = list(dimensional_known_dict.values())[0]
        x0 = (target_vol * 4.0 * l_over_d**2 / np.pi)**(1.0/3.0)
        res = minimize(lambda L_: abs(volume_for_length(L_[0]) - target_vol), x0=[x0], bounds=[(1e-6, None)], method='SLSQP')
        L = float(res.x[0])
        Vol = volume_for_length(L)

    idx_max = int(np.argmax(zeta))
    Psi_at_max = Psi[idx_max]

    return {
        'Psi': Psi, 'zeta': zeta,
        'L': L, 'Vol': Vol,
        'Psi_at_max': Psi_at_max, 'idx_max': idx_max,
        'S_te_used': S_te_used, 'adjusted_S_te': adjusted, 'S_te_original': S_te_orig
    }

# ------------------- Threshold finder that matches your exact argmax logic -------------------
def find_beta_threshold(R_le_over_c, Psi_zeta_max, zeta_max, zeta_te, step_deg=0.01, beta_min=0.0, beta_max=89.0):
    """
    Return (beta_last_ok, beta_first_bad) where beta_last_ok is the largest beta (multiple of step_deg)
    for which the discrete argmax index (on 1001 samples) remains at or before Psi_zeta_max,
    and beta_first_bad is the next sampled beta that pushes the argmax beyond Psi_zeta_max.
    This reproduces exactly the test in your original code:
        Psi_try[int(np.argmax(zeta_try))] <= Psi_zeta_max + 1e-9
    """
    # re-create same helper internals as in generate_streamlined_body_geometry
    N_1, N_2 = 0.5, 1.0
    def C(Psi): return (Psi**N_1) * ((1.0 - Psi)**N_2)
    def dC_dPsi(Psi):
        c = C(Psi)
        return c * (N_1 / Psi - N_2 / (1.0 - Psi))

    # same S_at and S' target
    C_at = C(Psi_zeta_max)
    S_at = (zeta_max - Psi_zeta_max * zeta_te) / C_at
    dC_at = dC_dPsi(Psi_zeta_max)
    S_prime_target = (- dC_at * S_at - zeta_te) / C_at

    def quadratic_coeffs(x0, y0, y0p, x1, y1):
        dx = x1 - x0
        a = (y1 - y0 - y0p * dx) / (dx**2)
        b = y0p - 2.0 * a * x0
        c = y0 + a * x0**2 - y0p * x0
        return a, b, c

    a_f, b_f, c_f = quadratic_coeffs(Psi_zeta_max, S_at, S_prime_target, 0.0, np.sqrt(2.0 * R_le_over_c))

    npts = 1001
    Psi_grid = np.linspace(0.0, 1.0, npts)
    idx_target = int(round(Psi_zeta_max * (npts - 1)))  # exactly the grid index used in your code

    betas = np.arange(beta_min, beta_max + 1e-12, step_deg)
    last_ok = None
    first_bad = None
    for beta in betas:
        S_te = np.tan(np.deg2rad(beta)) + zeta_te
        a_a, b_a, c_a = quadratic_coeffs(Psi_zeta_max, S_at, S_prime_target, 1.0, S_te)
        Svals = np.where(Psi_grid >= Psi_zeta_max, a_a*Psi_grid**2 + b_a*Psi_grid + c_a, a_f*Psi_grid**2 + b_f*Psi_grid + c_f)
        zeta = C(Psi_grid) * Svals + Psi_grid * zeta_te
        idx_global = int(np.argmax(zeta))
        if idx_global <= idx_target:
            last_ok = beta
        else:
            first_bad = beta
            break

    return last_ok, first_bad

# ------------------- Demo: compute threshold with step 0.01 deg and show below/above -------------------
if __name__ == "__main__":
    R_le_over_c = 0.06
    Psi_zeta_max = 0.35
    zeta_max = 1.0/5.0/2.0   # 0.1
    zeta_te = 0.0

    # find threshold with 0.01 degree resolution (you reported ~31.8..31.9)
    last_ok, first_bad = find_beta_threshold(R_le_over_c, Psi_zeta_max, zeta_max, zeta_te, step_deg=0.01)

    if last_ok is None:
        print("No beta >= 0 keeps the peak at or before Psi_zeta_max on the 1001-point grid.")
    else:
        print("Threshold interval (matching your code's argmax-on-1001-grid test):")
        print(f"  last beta that keeps peak at/before {Psi_zeta_max:.6f}  => {last_ok:.2f} deg")
        print(f"  first beta that moves peak past {Psi_zeta_max:.6f}    => {first_bad:.2f} deg")

        # show what your preserve_global_max does just below and just above
        beta_below = last_ok
        beta_above = first_bad

        res_below = generate_streamlined_body_geometry(R_le_over_c, beta_below, Psi_zeta_max, zeta_max, zeta_te, {'volume':6.0}, preserve_global_max=True)
        res_above = generate_streamlined_body_geometry(R_le_over_c, beta_above, Psi_zeta_max, zeta_max, zeta_te, {'volume':6.0}, preserve_global_max=True)

        print("\nPreserve behavior at/below/above threshold (using preserve_global_max=True):")
        print(f"  beta = {beta_below:.2f} -> adjusted S_te? {res_below['adjusted_S_te']}, S_te_used = {res_below['S_te_used']:.6g}, Psi_at_max = {res_below['Psi_at_max']:.6f}")
        print(f"  beta = {beta_above:.2f} -> adjusted S_te? {res_above['adjusted_S_te']}, S_te_used = {res_above['S_te_used']:.6g}, Psi_at_max = {res_above['Psi_at_max']:.6f}")

        # Plot below-threshold geometry (your plotting snippet)
        Psi = res_below['Psi']
        zeta = res_below['zeta']
        fig, ax = plt.subplots()
        ax.plot(Psi, zeta, color='blue', label=f'beta={beta_below:.2f} (below thresh)', linestyle='dashed')
        ax.set_xlabel('Chord (m)')
        ax.set_ylabel('Thickness (m)')
        ax.set_aspect('equal')
        ax.spines[['right', 'top']].set_visible(False)
        ax.tick_params(axis='x', which='both', top=False)
        ax.tick_params(axis='y', which='both', right=False)
        ax.tick_params(bottom=False)
        plt.show()

        # Plot above-threshold geometry
        Psi = res_above['Psi']
        zeta = res_above['zeta']
        fig, ax = plt.subplots()
        ax.plot(Psi, zeta, color='red', label=f'beta={beta_above:.2f} (above thresh)', linestyle='dashed')
        ax.set_xlabel('Chord (m)')
        ax.set_ylabel('Thickness (m)')
        ax.set_aspect('equal')
        ax.spines[['right', 'top']].set_visible(False)
        ax.tick_params(axis='x', which='both', top=False)
        ax.tick_params(axis='y', which='both', right=False)
        ax.tick_params(bottom=False)
        plt.show()
