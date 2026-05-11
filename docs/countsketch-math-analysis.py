"""
Monte Carlo simulations underlying docs/countsketch-math-analysis.md.

Reproduces:
1. Bias of raw squared sketched correlation r0^2
2. Bias of renormalized squared sketched cosine r'^2
3. Residual bias of the linear correction f(tilde r^2) = tilde r^2 * d/(d-2) - 1/(d-2)
4. Variance reduction from renormalization (raw / renormalized)
5. LD-score level bias and variance with realistic LD decay

Run with: /tmp/ldsc_py3_venv/bin/python3 docs/countsketch-math-analysis.py
(or any python with numpy).

The "active" estimator in ldsc-rs is the renormalized squared cosine with
the linear bias correction (compute.rs lines 799-818 + 1297-1306).
"""
import numpy as np


def make_pair(N, r, rng):
    """Generate two standardized columns of length N with population correlation r."""
    z1 = rng.standard_normal(N)
    z2 = rng.standard_normal(N)
    x = z1
    y = r * z1 + np.sqrt(max(0, 1 - r * r)) * z2
    x = (x - x.mean()) / x.std()
    y = (y - y.mean()) / y.std()
    return x, y


def countsketch(X, d, rng):
    """Apply CountSketch with sketch dim d. X shape (N, c) -> (d, c)."""
    N = X.shape[0]
    buckets = rng.integers(0, d, size=N)
    signs = rng.choice([-1, 1], size=N).astype(np.float64)
    out = np.zeros((d, X.shape[1]))
    np.add.at(out, buckets, signs[:, None] * X)
    return out


def simulation_1_per_pair_bias(N=2000, n_trials=30000):
    """§4: empirical bias of squared sketched estimators vs analytical predictions."""
    rng = np.random.default_rng(42)
    print("=" * 100)
    print("Simulation 1: per-pair bias of sketched squared estimators")
    print("=" * 100)
    h_renorm = "bias r'^2"
    print(f"{'r':>5} {'r2_emp':>8} {'d':>5} | "
          f"{'bias r0^2':>10} {'pred (1+r2)/d':>14} | "
          f"{h_renorm:>10} {'pred (1-r2)(1-2r2)/d':>22}")
    print("-" * 100)
    for r in [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]:
        x, y = make_pair(N, r, rng)
        r_emp = (x @ y) / N
        r2_emp = r_emp ** 2
        XY = np.stack([x, y], axis=1)
        for d in [50, 100, 200, 500]:
            bias_raw = 0.0
            bias_renorm = 0.0
            for _ in range(n_trials):
                t = countsketch(XY, d, rng)
                A = (t[:, 0] * t[:, 1]).sum()
                B = (t[:, 0] ** 2).sum()
                C = (t[:, 1] ** 2).sum()
                bias_raw += (A / N) ** 2 - r2_emp
                bias_renorm += A * A / (B * C) - r2_emp
            bias_raw /= n_trials
            bias_renorm /= n_trials
            pred_raw = (1 + r2_emp) / d
            pred_renorm = (1 - 3 * r2_emp + 2 * r2_emp ** 2) / d
            print(f"{r:>5.2f} {r2_emp:>8.4f} {d:>5} | "
                  f"{bias_raw:>10.5f} {pred_raw:>14.5f} | "
                  f"{bias_renorm:>10.5f} {pred_renorm:>22.5f}")
        print()


def simulation_2_correction_residual(N=2000, n_trials=30000):
    """§5: residual bias of the linear correction = r^2(2r^2-1)/(d-2)."""
    rng = np.random.default_rng(42)
    print("=" * 100)
    print("Simulation 2: residual bias of linear correction + variance reduction from renorm")
    print("=" * 100)
    print(f"{'r':>5} {'r2_emp':>8} {'d':>5} | "
          f"{'corrected':>10} {'r2_emp':>8} {'pred resid':>11} {'emp resid':>11} | "
          f"{'var ratio raw/renorm':>22}")
    print("-" * 100)
    for r in [0.0, 0.3, 0.5, 0.7, 0.9, 0.99]:
        x, y = make_pair(N, r, rng)
        r_emp = (x @ y) / N
        r2_emp = r_emp ** 2
        XY = np.stack([x, y], axis=1)
        for d in [50, 100, 200, 500]:
            corrected_sum = 0.0
            m_raw = 0.0
            m_renorm = 0.0
            v_raw = 0.0
            v_renorm = 0.0
            for _ in range(n_trials):
                t = countsketch(XY, d, rng)
                A = (t[:, 0] * t[:, 1]).sum()
                B = (t[:, 0] ** 2).sum()
                C = (t[:, 1] ** 2).sum()
                r0_sq = (A / N) ** 2
                rp_sq = A * A / (B * C)
                corr = rp_sq * d / (d - 2) - 1 / (d - 2)
                corrected_sum += corr
                m_raw += r0_sq
                m_renorm += rp_sq
                v_raw += r0_sq ** 2
                v_renorm += rp_sq ** 2
            corrected_sum /= n_trials
            m_raw /= n_trials
            m_renorm /= n_trials
            v_raw = v_raw / n_trials - m_raw ** 2
            v_renorm = v_renorm / n_trials - m_renorm ** 2
            emp_resid = corrected_sum - r2_emp
            pred_resid = r2_emp * (2 * r2_emp - 1) / (d - 2)
            ratio = v_raw / max(v_renorm, 1e-12)
            print(f"{r:>5.2f} {r2_emp:>8.4f} {d:>5} | "
                  f"{corrected_sum:>10.5f} {r2_emp:>8.5f} {pred_resid:>11.6f} {emp_resid:>11.6f} | "
                  f"{ratio:>22.2f}")
        print()


def simulation_3_ld_score(N=2000, M_window=200, n_trials=500):
    """§7: LD-score-level bias and variance, realistic LD decay."""
    rng = np.random.default_rng(42)
    print("=" * 100)
    print("Simulation 3: LD-score level (1 target + 200 neighbors with exp(-k/30) LD decay)")
    print("=" * 100)
    # Build the data: target SNP plus M neighbors with target correlation r_k = exp(-k/30)
    rs = np.exp(-np.arange(1, M_window + 1) / 30.0)
    z_target = rng.standard_normal(N)
    cols = [z_target]
    for r in rs:
        z = rng.standard_normal(N)
        cols.append(r * z_target + np.sqrt(max(0, 1 - r * r)) * z)
    X = np.stack(cols, axis=1)
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    target = X[:, 0]
    neighbors = X[:, 1:]
    r_emp = (target @ neighbors) / N
    l_exact = np.sum(r_emp ** 2)
    print(f"True empirical LD score (sum over {M_window} neighbors): {l_exact:.4f}")
    print()
    print(f"{'d':>5} | {'mean l_corr':>12} {'bias':>8} {'std l_corr':>12} {'CV':>7} | "
          f"{'bias uncorr':>13} {'bias raw':>10}")
    print("-" * 90)
    for d in [50, 100, 200, 500]:
        ls_corrected, ls_uncorrected, ls_raw = [], [], []
        for _ in range(n_trials):
            t = countsketch(X, d, rng)
            target_sketch = t[:, 0]
            neighbor_sketches = t[:, 1:]
            A = neighbor_sketches.T @ target_sketch
            B = (target_sketch ** 2).sum()
            C = (neighbor_sketches ** 2).sum(axis=0)
            rp_sq = A * A / (B * C)
            corr = rp_sq * d / (d - 2) - 1 / (d - 2)
            r0_sq = (A / N) ** 2
            ls_corrected.append(corr.sum())
            ls_uncorrected.append(rp_sq.sum())
            ls_raw.append(r0_sq.sum())
        ls_corrected = np.array(ls_corrected)
        ls_uncorrected = np.array(ls_uncorrected)
        ls_raw = np.array(ls_raw)
        bias_c = ls_corrected.mean() - l_exact
        bias_u = ls_uncorrected.mean() - l_exact
        bias_r = ls_raw.mean() - l_exact
        cv = ls_corrected.std() / ls_corrected.mean()
        print(f"{d:>5} | {ls_corrected.mean():>12.4f} {bias_c:>+8.4f} {ls_corrected.std():>12.4f} {cv:>7.4f} | "
              f"{bias_u:>+13.4f} {bias_r:>+10.4f}")


if __name__ == "__main__":
    simulation_1_per_pair_bias()
    print()
    simulation_2_correction_residual()
    print()
    simulation_3_ld_score()
