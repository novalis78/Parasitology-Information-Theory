#!/usr/bin/env python3
"""
Sensitivity Analysis & Unified Comparison Table
================================================
Recompute entropy at multiple resolutions with bootstrap CIs
for all three parasite systems. Build the paper's master table.
"""

import numpy as np
import pandas as pd
from scipy.stats import entropy
from collections import Counter
import gzip
import os
import json

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')

N_BOOTSTRAP = 1000
RANDOM_SEED = 42

# Published immune channel capacity benchmarks
CAPACITY_BENCHMARKS = [
    {'pathway': 'TNF→NF-κB (snapshot)', 'C_bits': 0.92, 'method': 'Binned MI, single-cell',
     'reference': 'Cheong et al. 2011', 'n_cells': '~9,000'},
    {'pathway': 'TNF→NF-κB (dynamic)', 'C_bits': 1.5, 'method': 'Dynamic MI, temporal codons',
     'reference': 'Selimkhanov et al. 2014', 'n_cells': '~thousands'},
    {'pathway': 'NF-κB signaling codons', 'C_bits': 2.6, 'method': '6 temporal codons',
     'reference': 'Adelaja et al. 2021', 'n_cells': '~thousands'},
    {'pathway': 'TCR kinetic proofreading', 'C_bits': 1.0, 'method': 'Computational model',
     'reference': 'Ganti et al. 2020', 'n_cells': 'model'},
    {'pathway': 'ERK pathway', 'C_bits': 1.0, 'method': 'Amplitude MI',
     'reference': 'Uda et al. 2013', 'n_cells': '~thousands'},
]


def compute_entropy_bits(counts):
    counts = np.array(counts, dtype=float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    probs = counts / counts.sum()
    return entropy(probs, base=2)


def bootstrap_entropy(counts, n_bootstrap=N_BOOTSTRAP, seed=RANDOM_SEED):
    rng = np.random.RandomState(seed)
    counts = np.array(counts, dtype=float)
    total = int(counts.sum())
    if total == 0:
        return 0.0, 0.0, 0.0
    types = np.repeat(np.arange(len(counts)), counts.astype(int))
    boot_H = []
    for _ in range(n_bootstrap):
        resample = rng.choice(types, size=total, replace=True)
        bc = np.bincount(resample, minlength=len(counts))
        boot_H.append(compute_entropy_bits(bc))
    boot_H = np.array(boot_H)
    return np.mean(boot_H), np.percentile(boot_H, 2.5), np.percentile(boot_H, 97.5)


def bootstrap_entropy_from_probs(probs, n_eff=1000, n_bootstrap=N_BOOTSTRAP, seed=RANDOM_SEED):
    rng = np.random.RandomState(seed)
    probs = np.array(probs, dtype=float)
    probs = probs[probs > 0]
    if len(probs) == 0:
        return 0.0, 0.0, 0.0
    probs = probs / probs.sum()
    boot_H = []
    for _ in range(n_bootstrap):
        counts = rng.multinomial(n_eff, probs)
        boot_probs = counts / counts.sum()
        boot_H.append(compute_entropy_bits(boot_probs))
    boot_H = np.array(boot_H)
    return np.mean(boot_H), np.percentile(boot_H, 2.5), np.percentile(boot_H, 97.5)


def miller_madow(n_types, n_samples):
    return (n_types - 1) / (2 * n_samples * np.log(2))


# ============================================================
# Analysis 1 Sensitivity: P. falciparum at multiple resolutions
# ============================================================
def sensitivity_pfalciparum():
    print("=" * 70)
    print("SENSITIVITY: P. falciparum var Gene Entropy")
    print("=" * 70)

    # Load data
    meta = pd.read_csv(os.path.join(DATA_DIR, 'Overview_All_Samples.txt'),
                       sep='\t', skiprows=4)
    meta = meta.rename(columns={'Origin of sample': 'Country', 'ID': 'SampleID'})
    sample_to_country = dict(zip(meta['SampleID'], meta['Country']))

    # Full domain architectures
    path_full = os.path.join(DATA_DIR, 'varDB.fulldataset.Domains_per_Gene.txt.gz')
    with gzip.open(path_full, 'rt') as f:
        header = f.readline().strip().split('\t')
        rows = []
        for line in f:
            parts = line.strip().split('\t')
            rows.append(parts)
    df = pd.DataFrame(rows, columns=header[:len(rows[0])] if len(header) >= len(rows[0]) else None)
    gene_col = df.columns[0]
    arch_col = df.columns[1]
    df['sample_id'] = df[gene_col].apply(lambda x: x.rsplit('.g', 1)[0] if '.g' in x else x.rsplit('.', 1)[0])
    df['country'] = df['sample_id'].map(sample_to_country)
    df = df[df['country'].notna()].copy()

    # Normalized subdomains
    path_norm = os.path.join(DATA_DIR, 'varDB.Normalised.Subdomains.txt.gz')
    with gzip.open(path_norm, 'rt') as f:
        header_n = f.readline().strip().split('\t')
        rows_n = []
        for line in f:
            parts = line.strip().split('\t')
            rows_n.append(parts)
    df_norm = pd.DataFrame(rows_n, columns=header_n[:len(rows_n[0])] if len(header_n) >= len(rows_n[0]) else None)

    results = []

    # Resolution 1: Full domain architecture (finest, full dataset)
    arch_counts = Counter(df[arch_col])
    ca = np.array(list(arch_counts.values()))
    H = compute_entropy_bits(ca)
    _, ci_lo, ci_hi = bootstrap_entropy(ca)
    mm = miller_madow(len(arch_counts), sum(arch_counts.values()))
    results.append({
        'system': 'P. falciparum', 'level': 'Global (full architecture)',
        'H_bits': H, 'CI_lo': ci_lo, 'CI_hi': ci_hi,
        'H_corrected': H + mm, 'n_types': len(arch_counts),
        'n_observations': sum(arch_counts.values()),
        'data_source': 'Otto et al. 2019 (varDB, N=2,398 isolates)'
    })
    print(f"  Full architecture (global): H = {H:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]  "
          f"({len(arch_counts):,} types)")

    # Resolution 2: Per-country entropy (top 5 largest)
    top_countries = df['country'].value_counts().head(5).index.tolist()
    for country in top_countries:
        cdata = df[df['country'] == country]
        cc = Counter(cdata[arch_col])
        ca_c = np.array(list(cc.values()))
        H_c = compute_entropy_bits(ca_c)
        _, ci_lo_c, ci_hi_c = bootstrap_entropy(ca_c)
        results.append({
            'system': 'P. falciparum', 'level': f'Per-country ({country})',
            'H_bits': H_c, 'CI_lo': ci_lo_c, 'CI_hi': ci_hi_c,
            'H_corrected': H_c - miller_madow(len(cc), sum(cc.values())),
            'n_types': len(cc), 'n_observations': sum(cc.values()),
            'data_source': f'Otto et al. 2019 ({country})'
        })
        print(f"  {country}: H = {H_c:.3f} [{ci_lo_c:.3f}, {ci_hi_c:.3f}]  ({len(cc):,} types)")

    # Resolution 3: Per-genome entropy (mean ± std across isolates)
    isolate_H = []
    for sid, group in df.groupby('sample_id'):
        ic = Counter(group[arch_col])
        isolate_H.append(compute_entropy_bits(np.array(list(ic.values()))))
    mean_iso = np.mean(isolate_H)
    std_iso = np.std(isolate_H)
    # Bootstrap the mean
    rng = np.random.RandomState(RANDOM_SEED)
    boot_means = []
    for _ in range(N_BOOTSTRAP):
        sample = rng.choice(isolate_H, size=len(isolate_H), replace=True)
        boot_means.append(np.mean(sample))
    boot_means = np.array(boot_means)
    results.append({
        'system': 'P. falciparum', 'level': 'Per-genome (mean)',
        'H_bits': mean_iso, 'CI_lo': np.percentile(boot_means, 2.5),
        'CI_hi': np.percentile(boot_means, 97.5),
        'H_corrected': mean_iso, 'n_types': None,
        'n_observations': len(isolate_H),
        'data_source': f'Otto et al. 2019 (N={len(isolate_H)} isolates)'
    })
    print(f"  Per-genome mean: H = {mean_iso:.3f} [{np.percentile(boot_means, 2.5):.3f}, "
          f"{np.percentile(boot_means, 97.5):.3f}]  (N={len(isolate_H)} isolates)")

    # Resolution 4: First domain only (coarser)
    df['first_domain'] = df[arch_col].apply(lambda x: str(x).split('-')[1] if '-' in str(x) else str(x))
    fd_counts = Counter(df['first_domain'])
    ca_fd = np.array(list(fd_counts.values()))
    H_fd = compute_entropy_bits(ca_fd)
    _, ci_lo_fd, ci_hi_fd = bootstrap_entropy(ca_fd)
    results.append({
        'system': 'P. falciparum', 'level': 'Global (first domain only)',
        'H_bits': H_fd, 'CI_lo': ci_lo_fd, 'CI_hi': ci_hi_fd,
        'H_corrected': H_fd - miller_madow(len(fd_counts), sum(fd_counts.values())),
        'n_types': len(fd_counts), 'n_observations': sum(fd_counts.values()),
        'data_source': 'Otto et al. 2019 (DBLα type only)'
    })
    print(f"  First domain only: H = {H_fd:.3f} [{ci_lo_fd:.3f}, {ci_hi_fd:.3f}]  "
          f"({len(fd_counts)} types)")

    # Resolution 5: Normalized subdomain set (controlled 60-sample subset)
    sub_col = df_norm.columns[1] if len(df_norm.columns) > 1 else df_norm.columns[0]
    sub_counts = Counter(df_norm[sub_col])
    ca_sub = np.array(list(sub_counts.values()))
    H_sub = compute_entropy_bits(ca_sub)
    _, ci_lo_sub, ci_hi_sub = bootstrap_entropy(ca_sub)
    results.append({
        'system': 'P. falciparum', 'level': 'Normalized subset (subdomains)',
        'H_bits': H_sub, 'CI_lo': ci_lo_sub, 'CI_hi': ci_hi_sub,
        'H_corrected': H_sub - miller_madow(len(sub_counts), sum(sub_counts.values())),
        'n_types': len(sub_counts), 'n_observations': sum(sub_counts.values()),
        'data_source': 'Otto et al. 2019 (60-sample normalized set)'
    })
    print(f"  Normalized subdomains: H = {H_sub:.3f} [{ci_lo_sub:.3f}, {ci_hi_sub:.3f}]  "
          f"({len(sub_counts):,} types)")

    return results


# ============================================================
# Analysis 2 Sensitivity: T. brucei across mice and timepoints
# ============================================================
def sensitivity_tbrucei():
    print("\n" + "=" * 70)
    print("SENSITIVITY: T. brucei VSG Switching Entropy")
    print("=" * 70)

    results = []

    # Load all mice
    all_mouse_H = {}
    for i in range(1, 5):
        df = pd.read_csv(os.path.join(DATA_DIR, f'mugnier_database_s{i}.csv'))
        df = df[df['VSG'] != 'parasitemia'].copy()
        df['pct'] = pd.to_numeric(df['pct'], errors='coerce')
        df['day'] = pd.to_numeric(df['day'], errors='coerce')
        df = df.dropna(subset=['pct', 'day'])

        mouse_H_by_day = {}
        for day in sorted(df['day'].unique()):
            day_data = df[df['day'] == day]
            probs = day_data['pct'].values / 100.0
            probs = probs[probs > 0]
            probs = probs / probs.sum()
            H = compute_entropy_bits(probs)
            _, ci_lo, ci_hi = bootstrap_entropy_from_probs(probs)
            mouse_H_by_day[int(day)] = {'H': H, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
                                         'n_vsgs': len(probs)}

        all_mouse_H[f'Mouse{i}'] = mouse_H_by_day

        # Per-mouse overall
        all_H = [v['H'] for v in mouse_H_by_day.values()]
        mean_H = np.mean(all_H)
        max_H = max(all_H)

        # Bootstrap the per-mouse mean
        rng = np.random.RandomState(RANDOM_SEED + i)
        boot_means = [np.mean(rng.choice(all_H, size=len(all_H), replace=True))
                      for _ in range(N_BOOTSTRAP)]
        boot_means = np.array(boot_means)

        results.append({
            'system': 'T. brucei', 'level': f'Mouse {i} (mean across timepoints)',
            'H_bits': mean_H,
            'CI_lo': np.percentile(boot_means, 2.5),
            'CI_hi': np.percentile(boot_means, 97.5),
            'H_corrected': mean_H, 'n_types': df['VSG'].nunique(),
            'n_observations': len(df['day'].unique()),
            'data_source': f'Mugnier et al. 2015 (Mouse {i})'
        })
        print(f"  Mouse {i} mean: H = {mean_H:.3f} [{np.percentile(boot_means, 2.5):.3f}, "
              f"{np.percentile(boot_means, 97.5):.3f}]  ({df['VSG'].nunique()} VSGs)")

        results.append({
            'system': 'T. brucei', 'level': f'Mouse {i} (max timepoint)',
            'H_bits': max_H,
            'CI_lo': None, 'CI_hi': None,
            'H_corrected': max_H, 'n_types': None,
            'n_observations': None,
            'data_source': f'Mugnier et al. 2015 (Mouse {i})'
        })

    # Overall mean across all 4 mice
    all_all_H = []
    for mouse, days in all_mouse_H.items():
        all_all_H.extend([v['H'] for v in days.values()])
    overall_mean = np.mean(all_all_H)
    rng = np.random.RandomState(RANDOM_SEED)
    boot_overall = [np.mean(rng.choice(all_all_H, size=len(all_all_H), replace=True))
                    for _ in range(N_BOOTSTRAP)]
    boot_overall = np.array(boot_overall)

    results.append({
        'system': 'T. brucei', 'level': 'All mice (mean, all timepoints)',
        'H_bits': overall_mean,
        'CI_lo': np.percentile(boot_overall, 2.5),
        'CI_hi': np.percentile(boot_overall, 97.5),
        'H_corrected': overall_mean, 'n_types': None,
        'n_observations': len(all_all_H),
        'data_source': 'Mugnier et al. 2015 (4 mice, 8 timepoints each)'
    })
    print(f"  Overall mean: H = {overall_mean:.3f} [{np.percentile(boot_overall, 2.5):.3f}, "
          f"{np.percentile(boot_overall, 97.5):.3f}]  (N={len(all_all_H)} mouse-timepoints)")

    # Max across all mice
    overall_max = max(all_all_H)
    results.append({
        'system': 'T. brucei', 'level': 'All mice (max single timepoint)',
        'H_bits': overall_max,
        'CI_lo': None, 'CI_hi': None,
        'H_corrected': overall_max, 'n_types': None,
        'n_observations': None,
        'data_source': 'Mugnier et al. 2015'
    })
    print(f"  Overall max: H = {overall_max:.3f}")

    # Late infection (Mouse 3, days 96-105)
    df5 = pd.read_csv(os.path.join(DATA_DIR, 'mugnier_database_s5.csv'))
    df5 = df5[df5['VSG'] != 'parasitemia'].copy()
    df5['pct'] = pd.to_numeric(df5['pct'], errors='coerce')
    df5['day'] = pd.to_numeric(df5['day'], errors='coerce')
    df5 = df5.dropna(subset=['pct', 'day'])
    late_H = []
    for day in sorted(df5['day'].unique()):
        probs = df5[df5['day'] == day]['pct'].values / 100.0
        probs = probs[probs > 0]
        probs = probs / probs.sum()
        late_H.append(compute_entropy_bits(probs))
    mean_late = np.mean(late_H)
    results.append({
        'system': 'T. brucei', 'level': 'Late infection (days 96-105)',
        'H_bits': mean_late,
        'CI_lo': None, 'CI_hi': None,
        'H_corrected': mean_late, 'n_types': df5['VSG'].nunique(),
        'n_observations': len(late_H),
        'data_source': 'Mugnier et al. 2015 (Mouse 3 late)'
    })
    print(f"  Late infection mean: H = {mean_late:.3f}  ({df5['VSG'].nunique()} VSGs)")

    return results, all_mouse_H


# ============================================================
# Analysis 3 Sensitivity: HIV with subtype filtering attempt
# ============================================================
def sensitivity_hiv():
    print("\n" + "=" * 70)
    print("SENSITIVITY: HIV-1 env Entropy")
    print("=" * 70)

    AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
    H_MAX_AA = np.log2(20)

    # Parse FASTA
    sequences = []
    with open(os.path.join(DATA_DIR, 'HIV1_env_proteins.fasta')) as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sequences.append((header, ''.join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header:
            sequences.append((header, ''.join(seq_parts)))

    results = []

    # --- Attempt subtype B filtering ---
    subtype_b = [(h, s) for h, s in sequences if 'subtype B' in h or 'subtype b' in h]
    print(f"  Total sequences: {len(sequences)}")
    print(f"  Subtype B (header match): {len(subtype_b)}")

    # If not enough from header, try to deduplicate by keeping unique sequences
    # as a proxy for one-per-patient
    unique_seqs = {}
    for h, s in sequences:
        if s not in unique_seqs:
            unique_seqs[s] = (h, s)
    unique_list = list(unique_seqs.values())
    print(f"  Unique sequences (deduplicated): {len(unique_list)}")

    # Compute for all datasets: full, unique-only, subtype B (if available)
    datasets = [
        ('All sequences', sequences),
        ('Deduplicated (unique)', unique_list),
    ]
    if len(subtype_b) >= 50:
        datasets.append(('Subtype B only', subtype_b))

    for label, seqs in datasets:
        lengths = [len(s) for _, s in seqs]
        median_len = int(np.median(lengths))
        filtered = [(h, s) for h, s in seqs if abs(len(s) - median_len) < 0.1 * median_len]
        min_len = min(len(s) for _, s in filtered) if filtered else 0

        if len(filtered) < 20:
            print(f"  {label}: too few sequences ({len(filtered)}) after length filter")
            continue

        # Per-position entropy
        pos_H = []
        for pos in range(min_len):
            aa_counts = Counter()
            for _, seq in filtered:
                if pos < len(seq) and seq[pos] in AA_ALPHABET:
                    aa_counts[seq[pos]] += 1
            counts = np.array([aa_counts.get(aa, 0) for aa in AA_ALPHABET])
            if counts.sum() > 0:
                pos_H.append(compute_entropy_bits(counts))

        mean_pos = np.mean(pos_H) if pos_H else 0
        total_pos = sum(pos_H)

        # Bootstrap the mean positional entropy
        rng = np.random.RandomState(RANDOM_SEED)
        boot_means = [np.mean(rng.choice(pos_H, size=len(pos_H), replace=True))
                      for _ in range(N_BOOTSTRAP)]
        boot_means = np.array(boot_means)

        # Fraction of positions exceeding benchmarks
        n_above_092 = sum(1 for h in pos_H if h > 0.92)
        n_above_26 = sum(1 for h in pos_H if h > 2.6)

        results.append({
            'system': 'HIV-1 env', 'level': f'Mean per-position ({label})',
            'H_bits': mean_pos,
            'CI_lo': np.percentile(boot_means, 2.5),
            'CI_hi': np.percentile(boot_means, 97.5),
            'H_corrected': mean_pos, 'n_types': None,
            'n_observations': len(filtered),
            'data_source': f'NCBI ({label}, N={len(filtered)})'
        })
        print(f"  {label}: mean H(pos) = {mean_pos:.3f} [{np.percentile(boot_means, 2.5):.3f}, "
              f"{np.percentile(boot_means, 97.5):.3f}]  "
              f"({len(pos_H)} positions, {n_above_092}/{len(pos_H)} > 0.92 bits)")

        # 9-mer epitope entropy
        kmer_counts = Counter()
        for _, seq in seqs:
            for i in range(len(seq) - 8):
                kmer = seq[i:i+9]
                if all(aa in AA_ALPHABET for aa in kmer):
                    kmer_counts[kmer] += 1

        ca_k = np.array(list(kmer_counts.values()))
        H_9mer = compute_entropy_bits(ca_k)
        _, ci_lo_k, ci_hi_k = bootstrap_entropy(ca_k)

        results.append({
            'system': 'HIV-1 env', 'level': f'9-mer epitope entropy ({label})',
            'H_bits': H_9mer,
            'CI_lo': ci_lo_k, 'CI_hi': ci_hi_k,
            'H_corrected': H_9mer - miller_madow(len(kmer_counts), sum(kmer_counts.values())),
            'n_types': len(kmer_counts),
            'n_observations': sum(kmer_counts.values()),
            'data_source': f'NCBI ({label})'
        })
        print(f"  {label}: H(9-mer) = {H_9mer:.3f} [{ci_lo_k:.3f}, {ci_hi_k:.3f}]  "
              f"({len(kmer_counts):,} unique 9-mers)")

        # Max per-position entropy
        max_pos = max(pos_H) if pos_H else 0
        results.append({
            'system': 'HIV-1 env', 'level': f'Max per-position ({label})',
            'H_bits': max_pos,
            'CI_lo': None, 'CI_hi': None,
            'H_corrected': max_pos, 'n_types': None,
            'n_observations': None,
            'data_source': f'NCBI ({label})'
        })

    return results


# ============================================================
# Build unified comparison table
# ============================================================
def build_unified_table(pf_results, tb_results, hiv_results):
    print("\n" + "=" * 70)
    print("UNIFIED COMPARISON TABLE")
    print("=" * 70)

    # Select key rows for the paper table
    all_results = pf_results + tb_results + hiv_results

    # Add capacity benchmarks
    for bm in CAPACITY_BENCHMARKS:
        all_results.append({
            'system': 'Immune signaling',
            'level': bm['pathway'],
            'H_bits': None,
            'CI_lo': None, 'CI_hi': None,
            'H_corrected': None,
            'n_types': None,
            'n_observations': bm['n_cells'],
            'data_source': bm['reference'],
            'C_bits': bm['C_bits'],
        })

    # Build the table
    table_rows = []
    C_nfkb = 0.92

    for r in all_results:
        if r['system'] == 'Immune signaling':
            table_rows.append({
                'System': r['system'],
                'Quantity': r['level'],
                'H_bits': '—',
                'CI_95': '—',
                'C_bits': f"{r['C_bits']:.2f}",
                'Ratio_H_C': '—',
                'Source': r['data_source'],
            })
        else:
            H = r['H_bits']
            ci_str = f"[{r['CI_lo']:.2f}, {r['CI_hi']:.2f}]" if r['CI_lo'] is not None else '—'
            ratio = H / C_nfkb
            table_rows.append({
                'System': r['system'],
                'Quantity': r['level'],
                'H_bits': f"{H:.2f}",
                'CI_95': ci_str,
                'C_bits': f"{C_nfkb:.2f}",
                'Ratio_H_C': f"{ratio:.1f}×",
                'Source': r['data_source'],
            })

    table_df = pd.DataFrame(table_rows)

    # Print formatted
    print(f"\n{'System':<18s} {'Quantity':<40s} {'H (bits)':<10s} {'95% CI':<18s} "
          f"{'C (bits)':<10s} {'H/C':<8s} {'Source'}")
    print("-" * 140)
    for _, row in table_df.iterrows():
        print(f"{row['System']:<18s} {row['Quantity']:<40s} {row['H_bits']:<10s} "
              f"{row['CI_95']:<18s} {row['C_bits']:<10s} {row['Ratio_H_C']:<8s} "
              f"{row['Source']}")

    return table_df, all_results


# ============================================================
# Main
# ============================================================
def main():
    pf_results = sensitivity_pfalciparum()
    tb_results, tb_mouse_H = sensitivity_tbrucei()
    hiv_results = sensitivity_hiv()

    table_df, all_results = build_unified_table(pf_results, tb_results, hiv_results)

    # Save
    table_df.to_csv(os.path.join(RESULTS_DIR, 'unified_comparison_table.csv'), index=False)

    # Save full sensitivity results
    sensitivity_data = {
        'pfalciparum': [{k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in r.items()} for r in pf_results],
        'tbrucei': [{k: (round(v, 4) if isinstance(v, float) else v)
                    for k, v in r.items()} for r in tb_results],
        'hiv': [{k: (round(v, 4) if isinstance(v, float) else v)
                for k, v in r.items()} for r in hiv_results],
        'capacity_benchmarks': CAPACITY_BENCHMARKS,
    }
    with open(os.path.join(RESULTS_DIR, 'sensitivity_analysis.json'), 'w') as f:
        json.dump(sensitivity_data, f, indent=2, default=str)

    print(f"\n  Saved unified_comparison_table.csv")
    print(f"  Saved sensitivity_analysis.json")

    # Generate sensitivity figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Figure 10: Sensitivity - entropy across all resolutions with CIs
        fig, ax = plt.subplots(figsize=(14, 8))

        # Collect plottable entries
        plot_data = []
        for r in pf_results + tb_results + hiv_results:
            if r['CI_lo'] is not None:
                plot_data.append(r)

        # Sort by system then H
        plot_data.sort(key=lambda x: (x['system'], x['H_bits']))

        y_pos = np.arange(len(plot_data))
        labels = [f"{r['system']}\n{r['level']}" for r in plot_data]
        h_vals = [r['H_bits'] for r in plot_data]
        ci_los = [max(r['H_bits'] - r['CI_lo'], 0) for r in plot_data]
        ci_his = [max(r['CI_hi'] - r['H_bits'], 0) for r in plot_data]

        # Color by system
        colors = []
        for r in plot_data:
            if 'falciparum' in r['system']:
                colors.append('#2196F3')
            elif 'brucei' in r['system']:
                colors.append('#FF9800')
            else:
                colors.append('#9C27B0')

        ax.barh(y_pos, h_vals, color=colors, alpha=0.8,
                xerr=[ci_los, ci_his], capsize=3, ecolor='gray')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Shannon Entropy H (bits)', fontsize=12)

        # Benchmark lines
        ax.axvline(x=0.92, color='red', linestyle='--', linewidth=2, alpha=0.7,
                  label='NF-κB capacity (0.92 bits)')
        ax.axvline(x=2.6, color='orange', linestyle='--', linewidth=1.5, alpha=0.7,
                  label='NF-κB codons (2.6 bits)')

        ax.set_title('Sensitivity Analysis: Antigenic Source Entropy Across All Resolutions\n'
                     'with 95% Bootstrap Confidence Intervals', fontsize=13)
        ax.legend(fontsize=10, loc='lower right')

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig10_sensitivity_all_systems.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig10_sensitivity_all_systems.png")

    except Exception as e:
        print(f"  WARNING: Could not generate figure: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
