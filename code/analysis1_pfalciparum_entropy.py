#!/usr/bin/env python3
"""
Analysis 1: P. falciparum var Gene Antigenic Entropy
=====================================================
Compute Shannon entropy H(var type) over DBLα domain architecture
distributions from Otto et al. (2019) varDB dataset.

Key question: Does antigenic entropy EXCEED immune channel capacity
(~1 bit for NF-κB, ~2.6 bits for NF-κB signaling codons)?

Data: Otto et al. (2019) Wellcome Open Research 4:193
       - 2,398 isolates across 15 countries
       - 168,738 var genes with domain architectures
"""

import numpy as np
import pandas as pd
from scipy.stats import entropy
import gzip
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
    'MAPK/ERK pulsatile (Nałęcz-Jawecki 2023)': 6.0,  # bits/hour
}

N_BOOTSTRAP = 1000
RANDOM_SEED = 42


def load_sample_metadata():
    """Load sample-to-country mapping from Overview_All_Samples.txt"""
    path = os.path.join(DATA_DIR, 'Overview_All_Samples.txt')
    df = pd.read_csv(path, sep='\t', skiprows=4)
    # Rename for convenience
    df = df.rename(columns={'Origin of sample': 'Country', 'ID': 'SampleID'})
    print(f"Loaded metadata: {len(df)} samples")
    print(f"Columns: {list(df.columns)}")
    print(f"Countries: {df['Country'].nunique()} unique")
    print(f"  {df['Country'].value_counts().to_dict()}")
    return df


def load_domain_architectures():
    """Load full domain architecture dataset (168K+ var genes)"""
    path = os.path.join(DATA_DIR, 'varDB.fulldataset.Domains_per_Gene.txt.gz')
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        print(f"Domain file columns: {header}")
        rows = []
        for line in f:
            parts = line.strip().split('\t')
            rows.append(parts)

    df = pd.DataFrame(rows, columns=header[:len(rows[0])] if len(header) >= len(rows[0]) else None)
    print(f"Loaded domain architectures: {len(df)} var genes")
    return df


def load_normalized_subdomains():
    """Load normalized subdomain classifications (controlled 60-sample subset)"""
    path = os.path.join(DATA_DIR, 'varDB.Normalised.Subdomains.txt.gz')
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        print(f"Subdomain file columns: {header}")
        rows = []
        for line in f:
            parts = line.strip().split('\t')
            rows.append(parts)

    df = pd.DataFrame(rows, columns=header[:len(rows[0])] if len(header) >= len(rows[0]) else None)
    print(f"Loaded normalized subdomains: {len(df)} entries")
    return df


def extract_sample_id(gene_id):
    """Extract sample ID from gene ID (e.g., 'PF0366-C.g495' -> 'PF0366-C')"""
    # Gene IDs have format SAMPLEID.gNNN
    parts = gene_id.rsplit('.g', 1)
    if len(parts) == 2:
        return parts[0]
    # Try alternative format
    parts = gene_id.rsplit('.', 1)
    return parts[0]


def compute_entropy_bits(counts):
    """Compute Shannon entropy in bits from count array"""
    counts = np.array(counts, dtype=float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    probs = counts / counts.sum()
    return entropy(probs, base=2)


def bootstrap_entropy(counts, n_bootstrap=N_BOOTSTRAP, seed=RANDOM_SEED):
    """Bootstrap Shannon entropy estimate with 95% CI"""
    rng = np.random.RandomState(seed)
    counts = np.array(counts, dtype=float)
    total = int(counts.sum())
    if total == 0:
        return 0.0, 0.0, 0.0

    # Create array of type indices for resampling
    types = np.repeat(np.arange(len(counts)), counts.astype(int))

    boot_entropies = []
    for _ in range(n_bootstrap):
        resample = rng.choice(types, size=total, replace=True)
        boot_counts = np.bincount(resample, minlength=len(counts))
        boot_entropies.append(compute_entropy_bits(boot_counts))

    boot_entropies = np.array(boot_entropies)
    return (
        np.mean(boot_entropies),
        np.percentile(boot_entropies, 2.5),
        np.percentile(boot_entropies, 97.5)
    )


def miller_madow_correction(n_types, n_samples):
    """Miller-Madow finite-sample bias correction for entropy"""
    return (n_types - 1) / (2 * n_samples * np.log(2))


# ============================================================
# Main Analysis
# ============================================================
def main():
    print("=" * 70)
    print("ANALYSIS 1: P. falciparum var Gene Antigenic Entropy")
    print("=" * 70)

    # ----------------------------------------------------------
    # Step 1: Load data
    # ----------------------------------------------------------
    print("\n--- Loading data ---")
    meta = load_sample_metadata()
    domains = load_domain_architectures()

    # Identify the gene ID and domain architecture columns
    print(f"\nFirst few rows of domain data:")
    print(domains.head())

    # Extract sample IDs from gene IDs
    gene_id_col = domains.columns[0]
    domains['sample_id'] = domains[gene_id_col].apply(extract_sample_id)

    # Identify domain architecture column (should be the one with NTS-DBL... strings)
    arch_col = None
    for col in domains.columns:
        sample_vals = domains[col].head(10).tolist()
        if any('DBL' in str(v) or 'NTS' in str(v) or 'CIDR' in str(v) for v in sample_vals):
            arch_col = col
            break

    if arch_col is None:
        # Use the last column as architecture
        arch_col = domains.columns[-1]
        print(f"WARNING: Could not auto-detect architecture column, using: {arch_col}")
    else:
        print(f"Domain architecture column: {arch_col}")

    print(f"Sample values: {domains[arch_col].head(5).tolist()}")

    # Map samples to countries
    sample_to_country = dict(zip(meta['SampleID'], meta['Country']))
    domains['country'] = domains['sample_id'].map(sample_to_country)

    n_mapped = domains['country'].notna().sum()
    n_total = len(domains)
    print(f"\nMapped {n_mapped}/{n_total} genes to countries ({100*n_mapped/n_total:.1f}%)")

    # Filter to mapped genes
    domains_mapped = domains[domains['country'].notna()].copy()

    # ----------------------------------------------------------
    # Step 2: Global entropy (all populations pooled)
    # ----------------------------------------------------------
    print("\n--- Global Entropy (all populations pooled) ---")

    global_arch_counts = Counter(domains_mapped[arch_col])
    n_unique_global = len(global_arch_counts)
    n_genes_global = sum(global_arch_counts.values())

    counts_array = np.array(list(global_arch_counts.values()))
    H_global = compute_entropy_bits(counts_array)
    H_max_global = np.log2(n_unique_global)
    efficiency_global = H_global / H_max_global if H_max_global > 0 else 0

    boot_mean, boot_lo, boot_hi = bootstrap_entropy(counts_array)
    mm_correction = miller_madow_correction(n_unique_global, n_genes_global)

    print(f"Total var genes: {n_genes_global:,}")
    print(f"Unique domain architectures: {n_unique_global:,}")
    print(f"H(global) = {H_global:.2f} bits")
    print(f"H_max = log2({n_unique_global}) = {H_max_global:.2f} bits")
    print(f"Efficiency H/H_max = {efficiency_global:.3f} ({100*efficiency_global:.1f}%)")
    print(f"Bootstrap 95% CI: [{boot_lo:.2f}, {boot_hi:.2f}] bits")
    print(f"Miller-Madow correction: {mm_correction:.4f} bits")
    print(f"H(corrected) = {H_global - mm_correction:.2f} bits")

    # ----------------------------------------------------------
    # Step 3: Per-country entropy
    # ----------------------------------------------------------
    print("\n--- Per-Country Entropy ---")

    country_results = []
    for country in sorted(domains_mapped['country'].unique()):
        country_data = domains_mapped[domains_mapped['country'] == country]
        arch_counts = Counter(country_data[arch_col])
        n_unique = len(arch_counts)
        n_genes = sum(arch_counts.values())
        n_samples = country_data['sample_id'].nunique()

        ca = np.array(list(arch_counts.values()))
        H = compute_entropy_bits(ca)
        H_max = np.log2(n_unique) if n_unique > 1 else 0
        eff = H / H_max if H_max > 0 else 0

        boot_mean, boot_lo, boot_hi = bootstrap_entropy(ca)
        mm = miller_madow_correction(n_unique, n_genes)

        country_results.append({
            'country': country,
            'n_isolates': n_samples,
            'n_genes': n_genes,
            'n_unique_types': n_unique,
            'H_bits': round(H, 3),
            'H_corrected': round(H - mm, 3),
            'H_max_bits': round(H_max, 3),
            'efficiency': round(eff, 4),
            'bootstrap_mean': round(boot_mean, 3),
            'bootstrap_CI_lo': round(boot_lo, 3),
            'bootstrap_CI_hi': round(boot_hi, 3),
            'MM_correction': round(mm, 4),
        })

        print(f"  {country:15s}: H = {H:.2f} bits  "
              f"[{boot_lo:.2f}, {boot_hi:.2f}]  "
              f"({n_unique:,} types, {n_genes:,} genes, {n_samples} isolates)  "
              f"H_max = {H_max:.2f}")

    country_df = pd.DataFrame(country_results)

    # ----------------------------------------------------------
    # Step 4: Per-genome entropy (theoretical)
    # ----------------------------------------------------------
    print("\n--- Per-Genome Entropy ---")

    genes_per_isolate = domains_mapped.groupby('sample_id')[arch_col].count()
    mean_genes = genes_per_isolate.mean()
    median_genes = genes_per_isolate.median()

    # Per-genome entropy assuming ~60 var genes with uniform expression
    H_per_genome_theoretical = np.log2(60)
    print(f"Genes per isolate: mean={mean_genes:.1f}, median={median_genes:.1f}")
    print(f"Theoretical per-genome entropy (uniform over 60): {H_per_genome_theoretical:.2f} bits")

    # Actual per-isolate entropy (how diverse are var genes within each genome?)
    isolate_entropies = []
    for sample_id, group in domains_mapped.groupby('sample_id'):
        arch_counts = Counter(group[arch_col])
        ca = np.array(list(arch_counts.values()))
        H = compute_entropy_bits(ca)
        n_genes = len(group)
        n_unique = len(arch_counts)
        isolate_entropies.append({
            'sample_id': sample_id,
            'country': group['country'].iloc[0],
            'n_genes': n_genes,
            'n_unique_architectures': n_unique,
            'H_bits': H,
            'H_max': np.log2(n_unique) if n_unique > 1 else 0,
        })

    isolate_df = pd.DataFrame(isolate_entropies)
    print(f"\nPer-isolate entropy stats:")
    print(f"  Mean H = {isolate_df['H_bits'].mean():.2f} bits")
    print(f"  Median H = {isolate_df['H_bits'].median():.2f} bits")
    print(f"  Std H = {isolate_df['H_bits'].std():.2f} bits")
    print(f"  Min H = {isolate_df['H_bits'].min():.2f}, Max H = {isolate_df['H_bits'].max():.2f}")
    print(f"  Mean unique architectures per genome: {isolate_df['n_unique_architectures'].mean():.1f}")

    # ----------------------------------------------------------
    # Step 5: Comparison to immune channel capacity
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print("KEY COMPARISON: Antigenic Entropy vs. Immune Channel Capacity")
    print("=" * 70)

    print(f"\n{'Quantity':<50s} {'Value (bits)':>12s}")
    print("-" * 65)
    print(f"{'P. falciparum global antigenic entropy':<50s} {H_global:>12.2f}")
    print(f"{'P. falciparum mean per-country entropy':<50s} {country_df['H_bits'].mean():>12.2f}")
    print(f"{'P. falciparum mean per-genome entropy':<50s} {isolate_df['H_bits'].mean():>12.2f}")
    print(f"{'Theoretical per-genome (60 var genes)':<50s} {H_per_genome_theoretical:>12.2f}")
    print("-" * 65)
    for name, cap in BENCHMARKS.items():
        print(f"{name:<50s} {cap:>12.2f}")
    print("-" * 65)

    # Ratio: antigenic entropy / channel capacity
    nfkb_cap = BENCHMARKS['NF-κB snapshot (Cheong 2011)']
    print(f"\nGlobal antigenic entropy / NF-κB capacity = "
          f"{H_global:.2f} / {nfkb_cap:.2f} = {H_global/nfkb_cap:.1f}x")
    print(f"Per-country mean entropy / NF-κB capacity = "
          f"{country_df['H_bits'].mean():.2f} / {nfkb_cap:.2f} = "
          f"{country_df['H_bits'].mean()/nfkb_cap:.1f}x")
    print(f"Per-genome mean entropy / NF-κB capacity = "
          f"{isolate_df['H_bits'].mean():.2f} / {nfkb_cap:.2f} = "
          f"{isolate_df['H_bits'].mean()/nfkb_cap:.1f}x")

    if H_global > nfkb_cap:
        print(f"\n>>> RESULT: Antigenic source entropy EXCEEDS immune channel capacity.")
        print(f"    H(antigen) = {H_global:.2f} bits > C(immune) = {nfkb_cap:.2f} bits")
        print(f"    By Shannon's converse theorem, reliable antigen discrimination is")
        print(f"    information-theoretically impossible. The parasite forces the immune")
        print(f"    system into the R > C regime where classification errors are guaranteed.")
        print(f"    Ratio H/C = {H_global/nfkb_cap:.1f}x at the population level.")
    else:
        print(f"\n>>> RESULT: Antigenic entropy does NOT exceed immune channel capacity.")
        print(f"    The framework prediction is NOT supported at the population level.")

    # ----------------------------------------------------------
    # Step 6: Sensitivity analysis - entropy at different resolutions
    # ----------------------------------------------------------
    print("\n--- Sensitivity Analysis: Entropy at Different Resolutions ---")

    # Resolution 1: Full domain architecture string (finest)
    # Already computed above

    # Resolution 2: First domain only (DBLα type)
    domains_mapped['first_domain'] = domains_mapped[arch_col].apply(
        lambda x: str(x).split('-')[1] if '-' in str(x) else str(x)
    )
    first_domain_counts = Counter(domains_mapped['first_domain'])
    H_first = compute_entropy_bits(np.array(list(first_domain_counts.values())))

    # Resolution 3: Number of domains
    domains_mapped['n_domains'] = domains_mapped[arch_col].apply(
        lambda x: len(str(x).split('-'))
    )
    n_domain_counts = Counter(domains_mapped['n_domains'])
    H_ndomain = compute_entropy_bits(np.array(list(n_domain_counts.values())))

    print(f"  Full architecture: H = {H_global:.2f} bits ({n_unique_global:,} types)")
    print(f"  First domain only: H = {H_first:.2f} bits ({len(first_domain_counts)} types)")
    print(f"  Domain count only: H = {H_ndomain:.2f} bits ({len(n_domain_counts)} types)")

    # ----------------------------------------------------------
    # Step 7: Save results
    # ----------------------------------------------------------
    print("\n--- Saving results ---")

    # Summary results
    summary = {
        'analysis': 'P. falciparum var gene antigenic entropy',
        'data_source': 'Otto et al. (2019) varDB dataset',
        'n_isolates': int(domains_mapped['sample_id'].nunique()),
        'n_var_genes': int(n_genes_global),
        'n_unique_architectures': int(n_unique_global),
        'global_entropy_bits': round(H_global, 3),
        'global_entropy_CI_lo': round(boot_lo, 3),
        'global_entropy_CI_hi': round(boot_hi, 3),
        'global_H_max_bits': round(H_max_global, 3),
        'global_efficiency': round(efficiency_global, 4),
        'mean_per_country_entropy': round(country_df['H_bits'].mean(), 3),
        'mean_per_genome_entropy': round(isolate_df['H_bits'].mean(), 3),
        'theoretical_per_genome_entropy': round(H_per_genome_theoretical, 3),
        'NF_kB_channel_capacity': nfkb_cap,
        'ratio_H_to_C_global': round(H_global / nfkb_cap, 2),
        'ratio_H_to_C_per_country': round(country_df['H_bits'].mean() / nfkb_cap, 2),
        'ratio_H_to_C_per_genome': round(isolate_df['H_bits'].mean() / nfkb_cap, 2),
        'interpretation': 'H(source) > C(channel) implies R > C regime: Shannon converse theorem guarantees unreliable discrimination',
        'sensitivity_full_arch_H': round(H_global, 3),
        'sensitivity_first_domain_H': round(H_first, 3),
        'sensitivity_domain_count_H': round(H_ndomain, 3),
    }

    with open(os.path.join(RESULTS_DIR, 'analysis1_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    country_df.to_csv(os.path.join(RESULTS_DIR, 'analysis1_per_country.csv'), index=False)
    isolate_df.to_csv(os.path.join(RESULTS_DIR, 'analysis1_per_isolate.csv'), index=False)

    print(f"  Saved analysis1_summary.json")
    print(f"  Saved analysis1_per_country.csv")
    print(f"  Saved analysis1_per_isolate.csv")

    # ----------------------------------------------------------
    # Step 8: Generate figures
    # ----------------------------------------------------------
    print("\n--- Generating figures ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        # Figure 1: Per-country entropy with immune capacity benchmarks
        fig, ax = plt.subplots(figsize=(12, 6))

        countries_sorted = country_df.sort_values('H_bits', ascending=True)
        y_pos = np.arange(len(countries_sorted))

        err_lo = np.maximum(countries_sorted['H_bits'].values - countries_sorted['bootstrap_CI_lo'].values, 0)
        err_hi = np.maximum(countries_sorted['bootstrap_CI_hi'].values - countries_sorted['H_bits'].values, 0)
        bars = ax.barh(y_pos, countries_sorted['H_bits'], color='#2196F3', alpha=0.8,
                       xerr=[err_lo, err_hi],
                       capsize=3, ecolor='gray')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(countries_sorted['country'], fontsize=10)
        ax.set_xlabel('Shannon Entropy H (bits)', fontsize=12)
        ax.set_title('P. falciparum var Gene Antigenic Entropy by Country\n'
                     '(Domain Architecture Distribution)', fontsize=13)

        # Add benchmark lines
        colors = ['red', 'orange', 'darkred', 'purple']
        for i, (name, cap) in enumerate(list(BENCHMARKS.items())[:4]):
            ax.axvline(x=cap, color=colors[i], linestyle='--', linewidth=1.5, alpha=0.7)
            short_name = name.split('(')[0].strip()
            ax.text(cap + 0.1, len(countries_sorted) - 1 - i * 0.8,
                    f'{short_name}\n({cap} bits)',
                    fontsize=7, color=colors[i], va='top')

        ax.set_xlim(0, max(countries_sorted['H_bits'].max() * 1.15,
                          max(list(BENCHMARKS.values())[:4]) * 1.1))
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig1_pfalciparum_entropy_by_country.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig1_pfalciparum_entropy_by_country.png")

        # Figure 2: Per-isolate entropy distribution
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(isolate_df['H_bits'], bins=50, color='#4CAF50', alpha=0.8,
                edgecolor='white', linewidth=0.5)
        ax.axvline(x=nfkb_cap, color='red', linestyle='--', linewidth=2,
                  label=f'NF-κB capacity ({nfkb_cap} bits)')
        ax.axvline(x=2.6, color='orange', linestyle='--', linewidth=2,
                  label='NF-κB codons (2.6 bits)')
        ax.axvline(x=isolate_df['H_bits'].mean(), color='blue', linestyle='-',
                  linewidth=2, label=f'Mean per-genome H ({isolate_df["H_bits"].mean():.2f} bits)')
        ax.set_xlabel('Shannon Entropy H (bits)', fontsize=12)
        ax.set_ylabel('Number of Isolates', fontsize=12)
        ax.set_title('Distribution of Per-Genome var Gene Entropy\n'
                     'P. falciparum (N={:,} isolates)'.format(len(isolate_df)),
                     fontsize=13)
        ax.legend(fontsize=10)
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig2_pfalciparum_per_genome_entropy.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig2_pfalciparum_per_genome_entropy.png")

        # Figure 3: Entropy comparison summary
        fig, ax = plt.subplots(figsize=(10, 5))

        categories = [
            'Per-genome\n(mean)',
            'Per-genome\n(theoretical)',
            'Per-country\n(mean)',
            'Global\n(all populations)',
        ]
        values = [
            isolate_df['H_bits'].mean(),
            H_per_genome_theoretical,
            country_df['H_bits'].mean(),
            H_global,
        ]

        x = np.arange(len(categories))
        bars = ax.bar(x, values, color=['#4CAF50', '#8BC34A', '#2196F3', '#1565C0'],
                     alpha=0.85, width=0.6)

        # Benchmark lines
        ax.axhline(y=nfkb_cap, color='red', linestyle='--', linewidth=2,
                  label=f'NF-κB capacity = {nfkb_cap} bits')
        ax.axhline(y=2.6, color='orange', linestyle='--', linewidth=1.5,
                  label='NF-κB codons = 2.6 bits')

        ax.set_xticks(x)
        ax.set_xticklabels(categories, fontsize=10)
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
        ax.set_title('P. falciparum Antigenic Entropy vs. Immune Channel Capacity',
                     fontsize=13)
        ax.legend(fontsize=10, loc='upper left')

        # Add value labels on bars
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                   f'{val:.1f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig3_entropy_vs_capacity.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig3_entropy_vs_capacity.png")

    except Exception as e:
        print(f"  WARNING: Could not generate figures: {e}")

    print("\n" + "=" * 70)
    print("ANALYSIS 1 COMPLETE")
    print("=" * 70)

    return summary, country_df, isolate_df


if __name__ == '__main__':
    main()
