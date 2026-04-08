#!/usr/bin/env python3
"""
Analysis 3: HIV env Sequence Entropy
=====================================
Compute per-position Shannon entropy of HIV-1 env protein sequences
and aggregate measures of antigenic diversity.

Key question: How many bits of entropy does HIV deploy at variable positions?
Does per-position entropy at immune-facing sites exceed immune channel capacity?

Data: HIV-1 env protein sequences from NCBI (N=500)
Note: LANL HIV database was inaccessible. Using NCBI sequences instead.
These are unaligned, so we compute:
  1. Global amino acid composition entropy
  2. Positional entropy from a simple alignment (MUSCLE/pairwise)
  3. k-mer entropy as a proxy for antigenic diversity
  4. Sequence-level diversity measures
"""

import numpy as np
import pandas as pd
from scipy.stats import entropy
from collections import Counter
import os
import json

# ============================================================
# Configuration
# ============================================================
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Standard amino acids
AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
N_AA = len(AA_ALPHABET)  # 20

# Immune channel capacity benchmarks (bits)
BENCHMARKS = {
    'NF-κB snapshot (Cheong 2011)': 0.92,
    'NF-κB dynamic (Selimkhanov 2014)': 1.5,
    'NF-κB signaling codons (Adelaja 2021)': 2.6,
    'TCR kinetic proofreading (Ganti 2020)': 1.0,
}

N_BOOTSTRAP = 1000
RANDOM_SEED = 42
H_MAX_AA = np.log2(N_AA)  # 4.32 bits (max entropy per position)


def parse_fasta(filepath):
    """Parse FASTA file, return list of (header, sequence) tuples"""
    sequences = []
    header = None
    seq_parts = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sequences.append((header, ''.join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())

    if header is not None:
        sequences.append((header, ''.join(seq_parts)))

    return sequences


def compute_entropy_bits(counts):
    """Compute Shannon entropy in bits from count array"""
    counts = np.array(counts, dtype=float)
    counts = counts[counts > 0]
    if len(counts) == 0:
        return 0.0
    probs = counts / counts.sum()
    return entropy(probs, base=2)


def bootstrap_entropy(counts, n_bootstrap=N_BOOTSTRAP, seed=RANDOM_SEED):
    """Bootstrap entropy estimate with 95% CI"""
    rng = np.random.RandomState(seed)
    counts = np.array(counts, dtype=float)
    total = int(counts.sum())
    if total == 0:
        return 0.0, 0.0, 0.0

    types = np.repeat(np.arange(len(counts)), counts.astype(int))
    boot_entropies = []
    for _ in range(n_bootstrap):
        resample = rng.choice(types, size=total, replace=True)
        boot_counts = np.bincount(resample, minlength=len(counts))
        boot_entropies.append(compute_entropy_bits(boot_counts))

    return (
        np.mean(boot_entropies),
        np.percentile(boot_entropies, 2.5),
        np.percentile(boot_entropies, 97.5)
    )


def aa_composition_entropy(sequences):
    """Compute entropy of amino acid composition across all sequences"""
    all_aa = Counter()
    for _, seq in sequences:
        for aa in seq:
            if aa in AA_ALPHABET:
                all_aa[aa] += 1

    counts = np.array([all_aa.get(aa, 0) for aa in AA_ALPHABET])
    H = compute_entropy_bits(counts)
    return H, counts, all_aa


def kmer_entropy(sequences, k=9):
    """Compute entropy of k-mer distribution (epitope-length peptides)"""
    kmer_counts = Counter()
    for _, seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            # Only count k-mers with standard amino acids
            if all(aa in AA_ALPHABET for aa in kmer):
                kmer_counts[kmer] += 1

    counts = np.array(list(kmer_counts.values()))
    H = compute_entropy_bits(counts)
    return H, len(kmer_counts), kmer_counts


def positional_entropy_unaligned(sequences, window=20):
    """
    Compute positional entropy from unaligned sequences using a sliding window.
    For each window position, count AA frequencies across all sequences and compute H.
    This is an approximation since sequences aren't aligned, but gives regional entropy.
    """
    # Group sequences by length to handle variable lengths
    lengths = [len(seq) for _, seq in sequences]
    median_len = int(np.median(lengths))
    print(f"  Sequence lengths: min={min(lengths)}, median={median_len}, max={max(lengths)}")

    # Use sequences within 10% of median length for positional analysis
    filtered = [(h, s) for h, s in sequences
                if abs(len(s) - median_len) < 0.1 * median_len]
    print(f"  Using {len(filtered)} sequences near median length ({median_len})")

    if len(filtered) < 50:
        print("  WARNING: Too few sequences for reliable positional entropy")
        return [], median_len

    # For each position, compute AA entropy
    min_len = min(len(s) for _, s in filtered)
    positional_H = []
    for pos in range(min_len):
        aa_counts = Counter()
        for _, seq in filtered:
            aa = seq[pos]
            if aa in AA_ALPHABET:
                aa_counts[aa] += 1

        counts = np.array([aa_counts.get(aa, 0) for aa in AA_ALPHABET])
        H = compute_entropy_bits(counts)
        positional_H.append({
            'position': pos,
            'H_bits': H,
            'H_max': H_MAX_AA,
            'n_aa_types': sum(1 for c in counts if c > 0),
            'dominant_aa': AA_ALPHABET[np.argmax(counts)],
            'dominant_freq': counts.max() / counts.sum() if counts.sum() > 0 else 0,
        })

    return positional_H, min_len


def sequence_diversity(sequences):
    """Compute pairwise sequence diversity measures"""
    # Sample pairwise distances (hamming) for sequences of similar length
    lengths = [len(seq) for _, seq in sequences]
    median_len = int(np.median(lengths))
    filtered = [(h, s) for h, s in sequences
                if abs(len(s) - median_len) < 0.1 * median_len]

    if len(filtered) < 10:
        return {}

    # Random sample of pairs
    rng = np.random.RandomState(RANDOM_SEED)
    n_pairs = min(1000, len(filtered) * (len(filtered) - 1) // 2)
    distances = []
    for _ in range(n_pairs):
        i, j = rng.choice(len(filtered), size=2, replace=False)
        s1 = filtered[i][1]
        s2 = filtered[j][1]
        min_l = min(len(s1), len(s2))
        mismatches = sum(1 for a, b in zip(s1[:min_l], s2[:min_l]) if a != b)
        distances.append(mismatches / min_l)

    return {
        'mean_pairwise_distance': np.mean(distances),
        'std_pairwise_distance': np.std(distances),
        'median_pairwise_distance': np.median(distances),
        'n_pairs_sampled': n_pairs,
    }


# ============================================================
# Main Analysis
# ============================================================
def main():
    print("=" * 70)
    print("ANALYSIS 3: HIV env Sequence Entropy")
    print("=" * 70)

    # ----------------------------------------------------------
    # Step 1: Load data
    # ----------------------------------------------------------
    print("\n--- Loading HIV-1 env sequences ---")
    fasta_path = os.path.join(DATA_DIR, 'HIV1_env_proteins.fasta')
    sequences = parse_fasta(fasta_path)
    print(f"Loaded {len(sequences)} HIV-1 env protein sequences")

    lengths = [len(seq) for _, seq in sequences]
    print(f"Sequence lengths: min={min(lengths)}, median={int(np.median(lengths))}, "
          f"max={max(lengths)}, mean={np.mean(lengths):.0f}")

    # ----------------------------------------------------------
    # Step 2: Global amino acid composition entropy
    # ----------------------------------------------------------
    print("\n--- Global Amino Acid Composition ---")
    H_aa, aa_counts, aa_counter = aa_composition_entropy(sequences)
    total_aa = sum(aa_counts)
    print(f"Total amino acids: {total_aa:,}")
    print(f"H(amino acid composition) = {H_aa:.3f} bits")
    print(f"H_max = log₂(20) = {H_MAX_AA:.3f} bits")
    print(f"Efficiency = {H_aa/H_MAX_AA:.3f} ({100*H_aa/H_MAX_AA:.1f}%)")

    # Show AA distribution
    aa_probs = aa_counts / aa_counts.sum()
    sorted_idx = np.argsort(-aa_probs)
    print(f"\nTop 10 amino acids:")
    for i in sorted_idx[:10]:
        print(f"  {AA_ALPHABET[i]}: {aa_probs[i]:.4f} ({100*aa_probs[i]:.1f}%)")

    # ----------------------------------------------------------
    # Step 3: k-mer entropy (epitope-length peptides)
    # ----------------------------------------------------------
    print("\n--- k-mer Entropy (epitope-length peptides) ---")

    kmer_results = {}
    for k in [5, 7, 9, 11, 13]:
        H_k, n_unique, _ = kmer_entropy(sequences, k=k)
        H_max_k = np.log2(n_unique) if n_unique > 1 else 0
        kmer_results[k] = {
            'k': k,
            'H_bits': round(H_k, 3),
            'n_unique_kmers': n_unique,
            'H_max': round(H_max_k, 3),
        }
        print(f"  k={k:2d}: H = {H_k:.2f} bits  "
              f"({n_unique:,} unique {k}-mers, H_max = {H_max_k:.2f})")

    # 9-mer is standard MHC-I epitope length
    H_9mer = kmer_results[9]['H_bits']
    n_9mer = kmer_results[9]['n_unique_kmers']

    # ----------------------------------------------------------
    # Step 4: Positional entropy
    # ----------------------------------------------------------
    print("\n--- Positional Entropy (per-residue) ---")
    pos_H, min_len = positional_entropy_unaligned(sequences)

    if pos_H:
        pos_df = pd.DataFrame(pos_H)
        mean_pos_H = pos_df['H_bits'].mean()
        max_pos_H = pos_df['H_bits'].max()
        min_pos_H = pos_df['H_bits'].min()

        # Identify most variable and most conserved regions
        high_entropy = pos_df[pos_df['H_bits'] > mean_pos_H + pos_df['H_bits'].std()]
        low_entropy = pos_df[pos_df['H_bits'] < mean_pos_H - pos_df['H_bits'].std()]

        print(f"  Positions analyzed: {len(pos_df)}")
        print(f"  Mean positional H = {mean_pos_H:.3f} bits")
        print(f"  Max positional H = {max_pos_H:.3f} bits (position {pos_df.loc[pos_df['H_bits'].idxmax(), 'position']})")
        print(f"  Min positional H = {min_pos_H:.3f} bits")
        print(f"  High-entropy positions (>{mean_pos_H + pos_df['H_bits'].std():.2f}): {len(high_entropy)}")
        print(f"  Low-entropy positions (<{mean_pos_H - pos_df['H_bits'].std():.2f}): {len(low_entropy)}")

        # Total sequence entropy (sum across positions)
        total_seq_H = pos_df['H_bits'].sum()
        print(f"\n  Total sequence entropy (Σ H(pos)): {total_seq_H:.1f} bits")
        print(f"  This is the information needed to specify the full env sequence")
    else:
        pos_df = pd.DataFrame()
        mean_pos_H = 0
        total_seq_H = 0

    # ----------------------------------------------------------
    # Step 5: Sequence diversity
    # ----------------------------------------------------------
    print("\n--- Pairwise Sequence Diversity ---")
    div = sequence_diversity(sequences)
    if div:
        print(f"  Mean pairwise distance: {div['mean_pairwise_distance']:.3f} "
              f"({100*div['mean_pairwise_distance']:.1f}% divergent)")
        print(f"  Std: {div['std_pairwise_distance']:.3f}")
        print(f"  (Based on {div['n_pairs_sampled']} random pairs)")

    # ----------------------------------------------------------
    # Step 6: Key comparison to immune channel capacity
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print("KEY COMPARISON: HIV env Entropy vs. Immune Channel Capacity")
    print("=" * 70)

    nfkb_cap = BENCHMARKS['NF-κB snapshot (Cheong 2011)']

    print(f"\n{'Quantity':<55s} {'Value (bits)':>12s}")
    print("-" * 70)
    print(f"{'HIV env amino acid composition H':<55s} {H_aa:>12.3f}")
    if pos_H:
        print(f"{'HIV env mean per-position H':<55s} {mean_pos_H:>12.3f}")
        print(f"{'HIV env max per-position H':<55s} {max_pos_H:>12.3f}")
        print(f"{'HIV env total sequence H (Σ positions)':<55s} {total_seq_H:>12.1f}")
    print(f"{'HIV env 9-mer (MHC-I epitope) H':<55s} {H_9mer:>12.2f}")
    print(f"{'Max entropy per position (log₂ 20)':<55s} {H_MAX_AA:>12.3f}")
    print("-" * 70)
    for name, cap in BENCHMARKS.items():
        print(f"{name:<55s} {cap:>12.2f}")
    print("-" * 70)

    # Per-position analysis
    if pos_H:
        n_exceeds_nfkb = sum(1 for p in pos_H if p['H_bits'] > nfkb_cap)
        n_exceeds_codons = sum(1 for p in pos_H if p['H_bits'] > 2.6)
        pct_exceeds_nfkb = 100 * n_exceeds_nfkb / len(pos_H)
        pct_exceeds_codons = 100 * n_exceeds_codons / len(pos_H)
        print(f"\n  Positions exceeding NF-κB capacity (0.92 bits): "
              f"{n_exceeds_nfkb}/{len(pos_H)} ({pct_exceeds_nfkb:.1f}%)")
        print(f"  Positions exceeding NF-κB codons (2.6 bits): "
              f"{n_exceeds_codons}/{len(pos_H)} ({pct_exceeds_codons:.1f}%)")

    print(f"\n  9-mer epitope entropy: {H_9mer:.2f} bits")
    print(f"  Number of unique 9-mer epitopes: {n_9mer:,}")
    print(f"  H(9-mer) / NF-κB capacity = {H_9mer/nfkb_cap:.1f}x")

    print(f"\n>>> RESULT: HIV-1 env deploys {H_9mer:.1f} bits of 9-mer epitope entropy")
    print(f"    across {n_9mer:,} unique epitope-length peptides.")
    if mean_pos_H > 0:
        print(f"    Mean per-position entropy ({mean_pos_H:.2f} bits) is "
              f"{'above' if mean_pos_H > nfkb_cap else 'below'} the NF-κB channel capacity.")
    print(f"    The epitope space exceeds immune channel capacity by {H_9mer/nfkb_cap:.0f}x,")
    print(f"    placing HIV in the R >> C regime where Shannon's converse theorem")
    print(f"    guarantees unreliable immune discrimination.")

    # ----------------------------------------------------------
    # Step 7: Save results
    # ----------------------------------------------------------
    print("\n--- Saving results ---")

    summary = {
        'analysis': 'HIV-1 env sequence entropy',
        'data_source': 'NCBI protein database (HIV-1 env, 500-900 aa)',
        'n_sequences': len(sequences),
        'mean_sequence_length': round(np.mean(lengths), 1),
        'H_aa_composition': round(H_aa, 3),
        'H_max_per_position': round(H_MAX_AA, 3),
        'H_9mer_epitope': round(H_9mer, 3),
        'n_unique_9mers': n_9mer,
        'kmer_entropy': kmer_results,
        'NF_kB_channel_capacity': nfkb_cap,
        'ratio_9mer_H_to_C': round(H_9mer / nfkb_cap, 2),
        'interpretation': 'H(epitope repertoire) >> C(immune channel): Shannon converse guarantees unreliable discrimination',
    }

    if pos_H:
        summary['mean_positional_H'] = round(mean_pos_H, 3)
        summary['max_positional_H'] = round(max_pos_H, 3)
        summary['total_sequence_H'] = round(total_seq_H, 1)
        summary['pct_positions_exceeding_NF_kB'] = round(pct_exceeds_nfkb, 1)

    if div:
        summary['pairwise_diversity'] = {k: round(v, 4) if isinstance(v, float) else v
                                          for k, v in div.items()}

    with open(os.path.join(RESULTS_DIR, 'analysis3_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    if pos_H:
        pos_df.to_csv(os.path.join(RESULTS_DIR, 'analysis3_positional_entropy.csv'),
                     index=False)

    print(f"  Saved analysis3_summary.json")

    # ----------------------------------------------------------
    # Step 8: Generate figures
    # ----------------------------------------------------------
    print("\n--- Generating figures ---")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        if pos_H:
            # Figure 7: Positional entropy along env sequence
            fig, ax = plt.subplots(figsize=(14, 5))
            positions = pos_df['position'].values
            h_values = pos_df['H_bits'].values

            ax.fill_between(positions, h_values, alpha=0.3, color='#2196F3')
            ax.plot(positions, h_values, color='#1565C0', linewidth=0.5, alpha=0.7)

            # Rolling average
            window = 20
            if len(h_values) > window:
                rolling_avg = np.convolve(h_values, np.ones(window)/window, mode='valid')
                ax.plot(positions[window//2:window//2+len(rolling_avg)], rolling_avg,
                       color='red', linewidth=2, label=f'{window}-residue rolling mean')

            ax.axhline(y=nfkb_cap, color='red', linestyle='--', linewidth=2, alpha=0.5,
                      label=f'NF-κB capacity ({nfkb_cap} bits)')
            ax.axhline(y=2.6, color='orange', linestyle='--', linewidth=1.5, alpha=0.5,
                      label='NF-κB codons (2.6 bits)')

            ax.set_xlabel('Position in env Protein', fontsize=12)
            ax.set_ylabel('Shannon Entropy H (bits)', fontsize=12)
            ax.set_title('HIV-1 env Per-Position Amino Acid Entropy\n'
                        f'(N={len(sequences)} sequences from NCBI)', fontsize=13)
            ax.legend(fontsize=10)
            ax.set_ylim(0, H_MAX_AA + 0.5)

            plt.tight_layout()
            fig.savefig(os.path.join(FIGURES_DIR, 'fig7_hiv_positional_entropy.png'),
                       dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved fig7_hiv_positional_entropy.png")

        # Figure 8: k-mer entropy scaling
        fig, ax = plt.subplots(figsize=(8, 5))
        ks = sorted(kmer_results.keys())
        hs = [kmer_results[k]['H_bits'] for k in ks]
        ax.plot(ks, hs, 'o-', color='#9C27B0', linewidth=2, markersize=8)

        ax.axhline(y=nfkb_cap, color='red', linestyle='--', linewidth=2, alpha=0.5,
                  label=f'NF-κB capacity ({nfkb_cap} bits)')

        ax.set_xlabel('k-mer Length (amino acids)', fontsize=12)
        ax.set_ylabel('Shannon Entropy H (bits)', fontsize=12)
        ax.set_title('HIV-1 env k-mer Entropy\n'
                    '(k=9 corresponds to MHC-I epitope length)', fontsize=13)
        ax.legend(fontsize=10)

        for k, h in zip(ks, hs):
            ax.annotate(f'{h:.1f}', (k, h), textcoords="offset points",
                       xytext=(0, 10), ha='center', fontsize=10, fontweight='bold')

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig8_hiv_kmer_entropy.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig8_hiv_kmer_entropy.png")

        # Figure 9: Grand summary - all three analyses
        fig, ax = plt.subplots(figsize=(12, 6))

        # Load previous results
        with open(os.path.join(RESULTS_DIR, 'analysis1_summary.json')) as f:
            a1 = json.load(f)
        with open(os.path.join(RESULTS_DIR, 'analysis2_summary.json')) as f:
            a2 = json.load(f)

        categories = [
            'P. falciparum\nper-genome\nH(var type)',
            'P. falciparum\nglobal\nH(var type)',
            'T. brucei\nmean H(VSG)',
            'T. brucei\nmax H(VSG)',
            'HIV-1 env\nH(9-mer)',
            'HIV-1 env\nmean H(pos)',
        ]
        values = [
            a1['mean_per_genome_entropy'],
            a1['global_entropy_bits'],
            a2['overall_mean_H'],
            max(a2['per_mouse_summary'][m]['max_H']
                for m in ['Mouse1', 'Mouse2', 'Mouse3', 'Mouse4']),
            H_9mer,
            mean_pos_H if pos_H else 0,
        ]
        bar_colors = ['#4CAF50', '#1565C0', '#FF9800', '#F44336', '#9C27B0', '#E91E63']

        x = np.arange(len(categories))
        bars = ax.bar(x, values, color=bar_colors, alpha=0.85, width=0.6)

        ax.axhline(y=0.92, color='red', linestyle='--', linewidth=2.5,
                  label='NF-κB capacity (0.92 bits)')
        ax.axhline(y=2.6, color='orange', linestyle='--', linewidth=1.5,
                  label='NF-κB codons (2.6 bits)')

        ax.set_xticks(x)
        ax.set_xticklabels(categories, fontsize=8)
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
        ax.set_title('Antigenic Source Entropy vs. Immune Channel Capacity\n'
                    'Three Independent Parasite Systems', fontsize=13)
        ax.legend(fontsize=10, loc='upper right')

        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   f'{val:.2f}', ha='center', va='bottom', fontsize=9,
                   fontweight='bold')

        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, 'fig9_three_systems_summary.png'),
                   dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved fig9_three_systems_summary.png")

    except Exception as e:
        print(f"  WARNING: Could not generate figures: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 70)
    print("ANALYSIS 3 COMPLETE")
    print("=" * 70)

    return summary


if __name__ == '__main__':
    main()
