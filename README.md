# The Parasitic Channel

**An Information-Theoretic Measurement of Immune Evasion Across Three Parasite Systems**

Lennart Lopin, Protogen Bio Corporation (April 2026)

## Summary

We compute the Shannon entropy of antigenic repertoires from three independent parasite systems and compare each against published measurements of immune signaling channel capacity. All three systems exceed the immune channel's capacity to discriminate -- placing the immune system in Shannon's R > C regime where classification errors are guaranteed.

| System | Antigenic Entropy H (bits) | 95% CI | H/C* |
|--------|---------------------------|--------|------|
| *P. falciparum* per-genome | 4.49 | [4.47, 4.52] | 4.9x |
| *P. falciparum* global | 6.79 | [6.76, 6.79] | 7.4x |
| *T. brucei* mean VSG (4 mice) | 1.70 | [1.34, 2.07] | 1.8x |
| *T. brucei* max timepoint | 3.99 | -- | 4.3x |
| HIV-1 per-position env | 2.98 | [2.94, 3.03] | 3.2x |
| HIV-1 9-mer epitope | 13.94 | [13.75, 13.77] | 15.1x |

*H/C computed against C = 0.92 bits (NF-kB snapshot capacity, Cheong et al. 2011)

## Repository Structure

```
code/                           # Analysis scripts
  analysis1_pfalciparum_entropy.py   # P. falciparum var gene entropy
  analysis2_tbrucei_vsg_entropy.py   # T. brucei VSG switching entropy
  analysis3_hiv_env_entropy.py       # HIV-1 env sequence entropy
  sensitivity_analysis.py            # Cross-system sensitivity & unified table

data/                           # Input datasets (see Data Sources below)
results/                        # Computed results (JSON, CSV)
figures/                        # Generated figures (PNG, 300 dpi)

parasitic_channel.tex           # LaTeX source
references.bib                  # BibTeX bibliography
The_Parasitic_Channel.pdf       # Compiled paper
```

## Data Sources

| Dataset | Source | Access |
|---------|--------|--------|
| *P. falciparum* var gene repertoires | Otto et al. (2019), Wellcome Open Research 4:193 | Zenodo: doi:10.5281/zenodo.3549732 |
| *T. brucei* VSG expression dynamics | Mugnier et al. (2015), Science 349:1470 | Supplementary databases S1-S5 |
| *T. brucei* VSG reference genome | Mugnier lab | github.com/mugnierlab/VSGSeqPipeline |
| HIV-1 env protein sequences | NCBI Protein Database | Entrez: txid11676[Organism] AND env[Gene Name] |
| Immune channel capacity benchmarks | Cheong et al. 2011, Selimkhanov et al. 2014, Adelaja et al. 2021 | Published measurements |

## Reproducing the Analysis

```bash
# Install dependencies
pip install numpy scipy pandas matplotlib seaborn scikit-learn

# Run analyses (in order)
python code/analysis1_pfalciparum_entropy.py
python code/analysis2_tbrucei_vsg_entropy.py
python code/analysis3_hiv_env_entropy.py
python code/sensitivity_analysis.py

# Compile paper
pdflatex parasitic_channel && bibtex parasitic_channel && pdflatex parasitic_channel && pdflatex parasitic_channel
```

## Companion Paper

This work is a companion to [The Democratic Channel](https://github.com/novalis78/The-Democratic-Channel) (Lopin, 2026), which measured the channel efficiency of democratic preference-policy transmission at 2.7% of Shannon capacity using the same information-theoretic framework applied to political data.

## License

The analysis code is released under the MIT License. The paper text and figures are copyright the author.
