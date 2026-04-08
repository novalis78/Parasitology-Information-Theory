# Data Sources and Provenance

All datasets used in this analysis are derived from published, peer-reviewed research
and are included here to facilitate replication. Please cite the original authors.

## Dataset 1: P. falciparum var gene repertoires

**Source:** Otto et al. (2019) "Evolutionary analysis of the most polymorphic gene family
in falciparum malaria." *Wellcome Open Research*, 4:193.

**License:** CC-BY 4.0 (Wellcome Open Research open data policy; Zenodo deposit)

**Original deposit:** https://doi.org/10.5281/zenodo.3549732

**Files:**
- `varDB.fulldataset.Domains_per_Gene.txt.gz` — 168,738 var genes with domain architectures
- `varDB.Normalised.Subdomains.txt.gz` — Subdomain classifications (60-sample normalized set)
- `varDB.Normalised.Domains.txt.gz` — Domain classifications (normalized set)
- `Overview_All_Samples.txt` — Sample metadata (2,398 isolates, 18 countries)
- `Additional_file_1.xlsx` — Supplementary tables S1-S7

## Dataset 2: T. brucei VSG expression dynamics

**Source:** Mugnier et al. (2015) "The in vivo dynamics of antigenic variation in
*Trypanosoma brucei*." *Science*, 349(6255):1470-1473.

**DOI:** 10.1126/science.aaa4502

**Files:**
- `mugnier_database_s1.csv` — Mouse 1 VSG expression (days 6-30)
- `mugnier_database_s2.csv` — Mouse 2 VSG expression (days 6-30)
- `mugnier_database_s3.csv` — Mouse 3 VSG expression (days 7-30)
- `mugnier_database_s4.csv` — Mouse 4 VSG expression (days 7-30)
- `mugnier_database_s5.csv` — Mouse 3 late infection (days 96-105)
- `Mugnier2015_SRA_metadata.tsv` — SRA sample metadata (SRP051697)

**Note:** Supplementary databases S1-S5 from the Science paper, provided for replication.
Reference VSG sequences available from https://github.com/mugnierlab/VSGSeqPipeline

## Dataset 3: HIV-1 envelope protein sequences

**Source:** NCBI Protein Database

**License:** Public domain (no use restrictions on NCBI/GenBank data)

**Query:** `txid11676[Organism] AND env[Gene Name] AND 500:900[Sequence Length]`

**Files:**
- `HIV1_env_proteins.fasta` — 500 HIV-1 env protein sequences (500-900 aa)

**Note:** These are cross-subtype sequences. For subtype-specific analyses,
filter by header annotations or obtain curated alignments from the
LANL HIV Sequence Database (https://www.hiv.lanl.gov/).
