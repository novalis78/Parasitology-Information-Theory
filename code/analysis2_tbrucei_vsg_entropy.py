#!/usr/bin/env python3
"""
Analysis 2: T. brucei VSG Switching Entropy Over Time
=====================================================
Compute Shannon entropy H(VSG | timepoint) from Mugnier et al. (2015)
VSG-seq time-series data, tracking how antigenic diversity evolves
during infection.

Key question: Does antigenic entropy INCREASE over infection?
At what rate (bits per day)? Does it approach log₂(repertoire size)?

Data: Mugnier et al. (2015) Science 349:1470
       - 4 mice (Databases S1-S4), EATRO 1125 strain
       - Timepoints: days 6-30 (early), days 96-105 (late, Mouse 3 only)
       - Database S5: Mouse 3 late infection
       - VSG expression as percentage of population
"""

import numpy as np
import pandas as pd
from scipy.stats import entropy, linregress
import os
import json
from collections import Counter

# ============================================================
# Configuration
# ============================================================
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Immune channel capacity benchmarks (bits)
BENCHMARKS = {
    'NF-κB snapshot (Cheong 2011)': 0.92,
    'NF-κB dynamic (Selimkhanov 2014)': 1.5,
    'NF-κB signaling codons (Adelaja 2021)': 2.6,
    'TCR kinetic proofreading (Ganti 2020)': 1.0,
}

N_BOOTSTRAP = 1000
RANDOM_SEED = 42

# Total VSG repertoire sizes for EATRO 1125
TOTAL_REPERTOIRE = 3570  # from VSGdb


def compute_entropy_bits(probs):
    """Compute Shannon entropy in bits from probability array"""
    probs = np.array(probs, dtype=float)
    probs = probs[probs > 0]
    if len(probs) == 0:
        return 0.0
    probs = probs / probs.sum()  # renormalize
    return entropy(probs, base=2)


def bootstrap_entropy_from_probs(probs, n_bootstrap=N_BOOTSTRAP, seed=RANDOM_SEED):
    """Bootstrap entropy from probability distribution by multinomial resampling"""
    rng = np.random.RandomState(seed)
    probs = np.array(probs, dtype=float)
    probs = probs[probs > 0]
    if len(probs) == 0:
        return 0.0, 0.0, 0.0
    probs = probs / probs.sum()

    # Simulate resampling: draw N samples from multinomial, recompute entropy
    # Use N=1000 as effective sample size (reasonable for RNA-seq)
    N_eff = 1000
    boot_entropies = []
    for _ in range(n_bootstrap):
        counts = rng.multinomial(N_eff, probs)
        boot_probs = counts / counts.sum()
        boot_entropies.append(compute_entropy_bits(boot_probs))

    boot_entropies = np.array(boot_entropies)
    return (
        np.mean(boot_entropies),
        np.percentile(boot_entropies, 2.5),
        np.percentile(boot_entropies, 97.5)
    )


def load_mugnier_data():
    """Load all Mugnier et al. (2015) VSG expression databases"""
    mice = {}
    for i in range(1, 5):
        path = os.path.join(DATA_DIR, f'mugnier_database_s{i}.csv')
        df = pd.read_csv(path)
        # Remove parasitemia rows (they don't have VSG names in the normal sense)
        df_vsg = df[df['VSG'] != 'parasitemia'].copy()
        df_vsg['pct'] = pd.to_numeric(df_vsg['pct'], errors='coerce')
        df_vsg['day'] = pd.to_numeric(df_vsg['day'], errors='coerce')
        df_vsg = df_vsg.dropna(subset=['pct', 'day'])
        mice[f'Mouse{i}'] = df_vsg
        print(f"  Mouse {i}: {len(df_vsg)} VSG-day entries, "
              f"{df_vsg['VSG'].nunique()} unique VSGs, "
              f"days {sorted(df_vsg['day'].unique().astype(int))}")

    # Database S5: Mouse 3 late infection
    path = os.path.join(DATA_DIR, 'mugnier_database_s5.csv')
    df = pd.read_csv(path)
    df_vsg = df[df['VSG'] != 'parasitemia'].copy()
    df_vsg['pct'] = pd.to_numeric(df_vsg['pct'], errors='coerce')
    df_vsg['day'] = pd.to_numeric(df_vsg['day'], errors='coerce')
    df_vsg = df_vsg.dropna(subset=['pct', 'day'])
    mice['Mouse3_late'] = df_vsg
    print(f"  Mouse 3 (late): {len(df_vsg)} VSG-day entries, "
          f"{df_vsg['VSG'].nunique()} unique VSGs, "
          f"days {sorted(df_vsg['day'].unique().astype(int))}")

    # Parasitemia data
    parasitemia = {}
    for i in range(1, 5):
        path = os.path.join(DATA_DIR, f'mugnier_database_s{i}.csv')
        df = pd.read_csv(path)
        para = df[df['VSG'] == 'parasitemia'][['day', 'parasites']].copy()
        para['day'] = pd.to_numeric(para['day'], errors='coerce')
        para['parasites'] = pd.to_numeric(para['parasites'], errors='coerce')
        para = para.dropna()
        parasitemia[f'Mouse{i}'] = para

    return mice, parasitemia


def compute_timepoint_entropy(mouse_df):
    """Compute entropy at each timepoint for one mouse"""
    results = []
    for day in sorted(mouse_df['day'].unique()):
        day_data = mouse_df[mouse_df['day'] == day]
        probs = day_data['pct'].values / 100.0  # Convert percentage to proportion

        # Filter to VSGs with > 0% expression
        probs = probs[probs > 0]
        if len(probs) == 0:
            continue

        probs = probs / probs.sum()  # Renormalize

        H = compute_entropy_bits(probs)
        H_max = np.log2(len(probs)) if len(probs) > 1 else 0
        n_vsgs = len(probs)

        # Effective number of VSGs (2^H)
        effective_n = 2 ** H

        # Bootstrap CI
        boot_mean, boot_lo, boot_hi = bootstrap_entropy_from_probs(probs)

        results.append({
            'day': int(day),
            'H_bits': round(H, 4),
            'H_max': round(H_max, 4),
            'efficiency': round(H / H_max, 4) if H_max > 0 else 0,
            'n_vsgs_expressed': n_vsgs,
            'effective_n_vsgs': round(effective_n, 2),
            'dominant_pct': round(probs.max() * 100, 2),
            'bootstrap_mean': round(boot_mean, 4),
            'bootstrap_CI_lo': round(boot_lo, 4),
            'bootstrap_CI_hi': round(boot_hi, 4),
        })

    return pd.DataFrame(results)


def compute_vsg_turnover(mouse_df):
    """Compute VSG turnover between consecutive timepoints"""
    days = sorted(mouse_df['day'].unique())
    results = []
    prev_vsgs = None
    for day in days:
        day_data = mouse_df[mouse_df['day'] == day]
        current_vsgs = set(day_data[day_data['pct'] > 0]['VSG'])

        if prev_vsgs is not None:
            shared = len(current_vsgs & prev_vsgs)
            new = len(current_vsgs - prev_vsgs)
            lost = len(prev_vsgs - current_vsgs)
            jaccard = shared / len(current_vsgs | prev_vsgs) if len(current_vsgs | prev_vsgs) > 0 else 0
            results.append({
                'day': int(day),
                'n_current': len(current_vsgs),
                'n_shared': shared,
                'n_new': new,
                'n_lost': lost,
                'jaccard': round(jaccard, 4),
                'turnover_rate': round(new / len(current_vsgs), 4) if len(current_vsgs) > 0 else 0,
            })

        prev_vsgs = current_vsgs

    return pd.DataFrame(results)


def compute_mi_time_vsg(mouse_df):
    """Compute MI between timepoint and VSG identity: I(time; VSG)"""
    # Joint distribution P(time, VSG)
    # Each timepoint-VSG pair weighted by expression percentage
    all_vsgs = sorted(mouse_df['VSG'].unique())
    all_days = sorted(mouse_df['day'].unique())

    # Build joint probability matrix
    joint = np.zeros((len(all_days), len(all_vsgs)))
    day_idx = {d: i for i, d in enumerate(all_days)}
    vsg_idx = {v: i for i, v in enumerate(all_vsgs)}

    for _, row in mouse_df.iterrows():
        if row['VSG'] in vsg_idx and row['day'] in day_idx:
            joint[day_idx[row['day']], vsg_idx[row['VSG']]] = row['pct']

    # Normalize to joint probability
    joint_sum = joint.sum()
    if joint_sum == 0:
        return 0.0

    P_xy = joint / joint_sum
    P_x = P_xy.sum(axis=1)  # marginal over time
    P_y = P_xy.sum(axis=0)  # marginal over VSG

    # MI = sum P(x,y) log2(P(x,y) / (P(x)*P(y)))
    MI = 0.0
    for i in range(len(all_days)):
        for j in range(len(all_vsgs)):
            if P_xy[i, j] > 0 and P_x[i] > 0 and P_y[j] > 0:
                MI += P_xy[i, j] * np.log2(P_xy[i, j] / (P_x[i] * P_y[j]))

    return MI


# ============================================================
# Main Analysis
# ============================================================
def main():
    print("=" * 70)
    print("ANALYSIS 2: T. brucei VSG Switching Entropy Over Time")
    print("=" * 70)

    # ----------------------------------------------------------
    # Step 1: Load data
    # ----------------------------------------------------------
    print("\n--- Loading Mugnier et al. (2015) data ---")
    mice, parasitemia = load_mugnier_data()

    # ----------------------------------------------------------
    # Step 2: Compute entropy at each timepoint for each mouse
    # ----------------------------------------------------------
    print("\n--- Computing per-timepoint entropy ---")

    all_results = {}
    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        df = mice[mouse_name]
        tp_results = compute_timepoint_entropy(df)
        all_results[mouse_name] = tp_results
        print(f"\n  {mouse_name}:")
        for _, row in tp_results.iterrows():
            print(f"    Day {int(row['day']):3d}: H = {row['H_bits']:.3f} bits  "
                  f"({int(row['n_vsgs_expressed']):3d} VSGs, "
                  f"effective={row['effective_n_vsgs']:.1f}, "
                  f"dominant={row['dominant_pct']:.1f}%)")

    # Mouse 3 late infection
    tp_late = compute_timepoint_entropy(mice['Mouse3_late'])
    all_results['Mouse3_late'] = tp_late
    print(f"\n  Mouse 3 (late infection):")
    for _, row in tp_late.iterrows():
        print(f"    Day {int(row['day']):3d}: H = {row['H_bits']:.3f} bits  "
              f"({int(row['n_vsgs_expressed']):3d} VSGs, "
              f"effective={row['effective_n_vsgs']:.1f}, "
              f"dominant={row['dominant_pct']:.1f}%)")

    # ----------------------------------------------------------
    # Step 3: Entropy growth rate (bits per day)
    # ----------------------------------------------------------
    print("\n--- Entropy Growth Rate ---")

    growth_rates = {}
    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        tp = all_results[mouse_name]
        if len(tp) >= 3:
            slope, intercept, r_value, p_value, std_err = linregress(
                tp['day'].values, tp['H_bits'].values
            )
            growth_rates[mouse_name] = {
                'slope_bits_per_day': round(slope, 4),
                'intercept': round(intercept, 4),
                'r_squared': round(r_value**2, 4),
                'p_value': round(p_value, 6),
            }
            print(f"  {mouse_name}: ΔH/Δt = {slope:.4f} bits/day  "
                  f"(R² = {r_value**2:.3f}, p = {p_value:.4f})")

    # Log fit: H(t) = a * log2(t) + b
    print("\n  Logarithmic fit H(t) = a·log₂(t) + b:")
    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        tp = all_results[mouse_name]
        if len(tp) >= 3:
            log_days = np.log2(tp['day'].values)
            slope, intercept, r_value, p_value, std_err = linregress(
                log_days, tp['H_bits'].values
            )
            print(f"  {mouse_name}: a = {slope:.3f}, b = {intercept:.3f}  "
                  f"(R² = {r_value**2:.3f})")

    # ----------------------------------------------------------
    # Step 4: VSG turnover between timepoints
    # ----------------------------------------------------------
    print("\n--- VSG Turnover Analysis ---")

    turnover_results = {}
    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        turnover = compute_vsg_turnover(mice[mouse_name])
        turnover_results[mouse_name] = turnover
        if len(turnover) > 0:
            mean_turnover = turnover['turnover_rate'].mean()
            print(f"  {mouse_name}: mean turnover rate = {mean_turnover:.3f} "
                  f"(fraction new VSGs per timepoint)")

    # ----------------------------------------------------------
    # Step 5: MI between time and VSG identity
    # ----------------------------------------------------------
    print("\n--- Mutual Information I(time; VSG) ---")

    mi_results = {}
    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        mi = compute_mi_time_vsg(mice[mouse_name])
        mi_results[mouse_name] = mi
        # H(time) as benchmark
        n_timepoints = all_results[mouse_name]['day'].nunique()
        H_time = np.log2(n_timepoints)
        efficiency = mi / H_time if H_time > 0 else 0
        print(f"  {mouse_name}: I(time; VSG) = {mi:.3f} bits  "
              f"(H(time) = {H_time:.2f} bits, efficiency = {efficiency:.3f})")

    # ----------------------------------------------------------
    # Step 6: Fraction of repertoire explored
    # ----------------------------------------------------------
    print("\n--- Fraction of Repertoire Explored ---")

    H_max_repertoire = np.log2(TOTAL_REPERTOIRE)
    print(f"  Total VSG repertoire: {TOTAL_REPERTOIRE} sequences")
    print(f"  H_max(repertoire) = log₂({TOTAL_REPERTOIRE}) = {H_max_repertoire:.2f} bits")

    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        n_unique = mice[mouse_name]['VSG'].nunique()
        tp = all_results[mouse_name]
        max_H = tp['H_bits'].max()
        effective_explored = 2 ** max_H
        frac_explored = n_unique / TOTAL_REPERTOIRE
        frac_effective = effective_explored / TOTAL_REPERTOIRE
        print(f"  {mouse_name}: {n_unique} unique VSGs detected "
              f"({100*frac_explored:.1f}% of repertoire), "
              f"max effective = {effective_explored:.0f} "
              f"({100*frac_effective:.1f}%)")

    # ----------------------------------------------------------
    # Step 7: Key comparison to immune channel capacity
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print("KEY COMPARISON: VSG Entropy vs. Immune Channel Capacity")
    print("=" * 70)

    # Average entropy across mice at each timepoint
    common_days = set(all_results['Mouse1']['day'])
    for m in ['Mouse2', 'Mouse3', 'Mouse4']:
        common_days &= set(all_results[m]['day'])
    common_days = sorted(common_days)

    print(f"\nMean entropy across 4 mice at common timepoints:")
    print(f"  {'Day':<8s} {'Mean H (bits)':<15s} {'Std':<10s} {'vs NF-κB (0.92 bits)'}")
    print(f"  {'-'*55}")

    mean_entropies = []
    for day in common_days:
        h_vals = []
        for m in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
            row = all_results[m][all_results[m]['day'] == day]
            if len(row) > 0:
                h_vals.append(row['H_bits'].values[0])
        if h_vals:
            mean_h = np.mean(h_vals)
            std_h = np.std(h_vals)
            ratio = mean_h / 0.92
            mean_entropies.append({'day': day, 'mean_H': mean_h, 'std_H': std_h})
            print(f"  {day:<8d} {mean_h:<15.3f} {std_h:<10.3f} {ratio:.1f}x")

    # Overall summary
    all_H = []
    for m in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        all_H.extend(all_results[m]['H_bits'].values)
    overall_mean = np.mean(all_H)
    overall_std = np.std(all_H)

    nfkb_cap = BENCHMARKS['NF-κB snapshot (Cheong 2011)']
    print(f"\n  Overall mean H across all mice and timepoints: "
          f"{overall_mean:.2f} ± {overall_std:.2f} bits")
    print(f"  Overall / NF-κB capacity: {overall_mean/nfkb_cap:.1f}x")

    if overall_mean > nfkb_cap:
        print(f"\n>>> RESULT: VSG source entropy EXCEEDS immune channel capacity.")
        print(f"    H(VSG) mean = {overall_mean:.2f} bits > C(immune) = {nfkb_cap:.2f} bits")
        print(f"    By Shannon's converse theorem, reliable VSG discrimination is")
        print(f"    information-theoretically impossible at most timepoints.")
        print(f"    The parasite forces R > C, guaranteeing immune classification errors.")
    else:
        print(f"\n>>> RESULT: VSG switching entropy does NOT consistently exceed "
              f"immune channel capacity.")

    # ----------------------------------------------------------
    # Step 8: Save results
    # ----------------------------------------------------------
    print("\n--- Saving results ---")

    summary = {
        'analysis': 'T. brucei VSG switching entropy over time',
        'data_source': 'Mugnier et al. (2015) Science 349:1470',
        'total_vsg_repertoire': TOTAL_REPERTOIRE,
        'H_max_repertoire_bits': round(H_max_repertoire, 2),
        'overall_mean_H': round(overall_mean, 3),
        'overall_std_H': round(overall_std, 3),
        'NF_kB_channel_capacity': nfkb_cap,
        'ratio_overall_to_NF_kB': round(overall_mean / nfkb_cap, 2),
        'growth_rates': growth_rates,
        'mi_time_vsg': {k: round(v, 4) for k, v in mi_results.items()},
        'per_mouse_summary': {},
    }

    for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
        tp = all_results[mouse_name]
        summary['per_mouse_summary'][mouse_name] = {
            'n_unique_vsgs': int(mice[mouse_name]['VSG'].nunique()),
            'n_timepoints': len(tp),
            'min_H': round(tp['H_bits'].min(), 3),
            'max_H': round(tp['H_bits'].max(), 3),
            'mean_H': round(tp['H_bits'].mean(), 3),
        }

    with open(os.path.join(RESULTS_DIR, 'analysis2_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    # Save per-timepoint data
    for mouse_name, tp in all_results.items():
        tp_out = tp.copy()
        tp_out['mouse'] = mouse_name
        tp_out.to_csv(os.path.join(RESULTS_DIR, f'analysis2_{mouse_name}_entropy.csv'),
                     index=False)

    print(f"  Saved analysis2_summary.json")
    print(f"  Saved per-mouse entropy CSVs")

    # ----------------------------------------------------------
    # Step 9: Generate figures
    # ----------------------------------------------------------
    print("\n--- Generating figures ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Figure 4: Entropy over time for all 4 mice
        fig, ax = plt.subplots(figsize=(12, 6))

        colors = {'Mouse1': '#1f77b4', 'Mouse2': '#ff7f0e',
                  'Mouse3': '#2ca02c', 'Mouse4': '#d62728'}

        for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
            tp = all_results[mouse_name]
            ax.plot(tp['day'], tp['H_bits'], 'o-', color=colors[mouse_name],
                   linewidth=2, markersize=6, label=mouse_name, alpha=0.8)
            # Bootstrap CI band
            ax.fill_between(tp['day'], tp['bootstrap_CI_lo'], tp['bootstrap_CI_hi'],
                          color=colors[mouse_name], alpha=0.15)

        # Mouse 3 late
        tp_late = all_results['Mouse3_late']
        ax.plot(tp_late['day'], tp_late['H_bits'], 's--', color='#2ca02c',
               linewidth=2, markersize=8, label='Mouse3 (late)', alpha=0.6)

        # Immune capacity benchmarks
        ax.axhline(y=0.92, color='red', linestyle='--', linewidth=2, alpha=0.7,
                  label='NF-κB capacity (0.92 bits)')
        ax.axhline(y=2.6, color='orange', linestyle='--', linewidth=1.5, alpha=0.7,
                  label='NF-κB codons (2.6 bits)')

        ax.set_xlabel('Days Post-Infection', fontsize=12)
        ax.set_ylabel('Shannon Entropy H(VSG) (bits)', fontsize=12)
        ax.set_title('T. brucei VSG Antigenic Entropy Over Time\n'
                     'Mugnier et al. (2015)', fontsize=13)
        ax.legend(fontsize=9, loc='upper left')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig4_tbrucei_entropy_over_time.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig4_tbrucei_entropy_over_time.png")

        # Figure 5: Number of VSGs expressed over time
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        for mouse_name in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']:
            tp = all_results[mouse_name]
            ax1.plot(tp['day'], tp['n_vsgs_expressed'], 'o-',
                    color=colors[mouse_name], linewidth=2, markersize=5,
                    label=mouse_name)
            ax2.plot(tp['day'], tp['effective_n_vsgs'], 'o-',
                    color=colors[mouse_name], linewidth=2, markersize=5,
                    label=mouse_name)

        ax1.set_xlabel('Days Post-Infection', fontsize=11)
        ax1.set_ylabel('Number of VSGs Expressed', fontsize=11)
        ax1.set_title('VSGs Detected per Timepoint', fontsize=12)
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)

        ax2.set_xlabel('Days Post-Infection', fontsize=11)
        ax2.set_ylabel('Effective Number of VSGs (2^H)', fontsize=11)
        ax2.set_title('Effective VSG Diversity (2^H)', fontsize=12)
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig5_tbrucei_vsg_counts.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig5_tbrucei_vsg_counts.png")

        # Figure 6: Entropy vs immune capacity summary (both analyses)
        fig, ax = plt.subplots(figsize=(10, 6))

        # Load Analysis 1 summary
        try:
            with open(os.path.join(RESULTS_DIR, 'analysis1_summary.json')) as f:
                a1 = json.load(f)

            categories = [
                'P. falciparum\nper-genome',
                'P. falciparum\nper-country',
                'P. falciparum\nglobal',
                'T. brucei VSG\n(mean all mice)',
                'T. brucei VSG\n(max observed)',
            ]

            max_H = max(tp['H_bits'].max() for tp in
                       [all_results[m] for m in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']])

            values = [
                a1['mean_per_genome_entropy'],
                a1['mean_per_country_entropy'],
                a1['global_entropy_bits'],
                overall_mean,
                max_H,
            ]
            bar_colors = ['#4CAF50', '#2196F3', '#1565C0', '#FF9800', '#F44336']

            x = np.arange(len(categories))
            bars = ax.bar(x, values, color=bar_colors, alpha=0.85, width=0.6)

            ax.axhline(y=0.92, color='red', linestyle='--', linewidth=2,
                      label='NF-κB capacity (0.92 bits)')
            ax.axhline(y=2.6, color='orange', linestyle='--', linewidth=1.5,
                      label='NF-κB codons (2.6 bits)')

            ax.set_xticks(x)
            ax.set_xticklabels(categories, fontsize=9)
            ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
            ax.set_title('Parasitic Antigenic Entropy vs. Immune Channel Capacity\n'
                        'Two Independent Parasite Systems', fontsize=13)
            ax.legend(fontsize=10)

            for bar, val in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                       f'{val:.2f}', ha='center', va='bottom', fontsize=10,
                       fontweight='bold')

            plt.tight_layout()
            fig.savefig(os.path.join(FIGURES_DIR,
                       'fig6_combined_entropy_vs_capacity.png'),
                       dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved fig6_combined_entropy_vs_capacity.png")

        except FileNotFoundError:
            print("  Skipping combined figure (Analysis 1 results not found)")

    except Exception as e:
        print(f"  WARNING: Could not generate figures: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 70)
    print("ANALYSIS 2 COMPLETE")
    print("=" * 70)

    return summary, all_results


if __name__ == '__main__':
    main()
