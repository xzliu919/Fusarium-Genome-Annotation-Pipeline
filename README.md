# Fusarium Genome Annotation Pipeline

A comprehensive pipeline for automated gene prediction and annotation of Fusarium genomes using BRAKER3 with multi-evidence integration (RNA-seq evidence and protein homology).

## Overview

This pipeline provides an end-to-end solution for annotating novel Fusarium genomes by integrating three lines of evidence:
- **RNA-seq data**: Direct transcription evidence from mapped reads
- **Protein homology**: Cross-species protein evidence from related taxa
- **Ab initio prediction**: AUGUSTUS gene model predictions

The pipeline implements an adaptive workflow that automatically selects the best available evidence and handles cases where certain data types are unavailable.

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                           FUSARIUM GENOME ANNOTATION                                │
└─────────────────────────────────────────────────────────────────────────────────────┘
                                      │
        ┌─────────────────────────────┼─────────────────────────────┐
        ▼                             ▼                             ▼
┌───────────────┐            ┌────────────────┐            ┌─────────────────┐
│   Step 1      │            │   Step 2       │            │   Step 3        │
│   Genome      │            │   RNA-seq      │            │   BRAKER3       │
│   Preprocessing│            │   Mapping      │            │   Annotation    │
└───────────────┘            └────────────────┘            └─────────────────┘
        │                             │                             │
        ▼                             ▼                             ▼
┌───────────────┐            ┌────────────────┐            ┌─────────────────┐
│  Repeat       │            │   STAR        │            │   Strategy A   │
│  Masking      │            │   Alignment   │            │   Full Evidence│
│               │            │               │            │   (BAM+Hints+  │
│  -BuildDatabase              │  -MapQ filter│            │    Protein)    │
│  -RepeatModeler│            │  -Merge BAMs │            │                │
│  -RepeatMasker│            │               │            │   Strategy B   │
│               │            │   StringTie   │            │   RNA+BAM only │
│  Soft-mask    │            │   Assembly    │            │                │
│  .fna.masked  │            │               │            │   Strategy C   │
└───────────────┘            │  -GTF hints   │            │   Protein only │
                              └────────────────┘            │                │
                                      │                      │   TSEBRA       │
                                      ▼                      │   Merging      │
                             ┌────────────────┐               │                │
                             │   Evidence    │               └─────────────────┘
                             │   Generation  │
                             └────────────────┘
                                      │
                                      ▼
                             ┌────────────────┐
                             │   Final Output │
                             ├────────────────┤
                             │ - final.gff3   │
                             │ - _cds.fa      │
                             │ - _pep.fa      │
                             └────────────────┘
```

## Directory Structure

```
.
├── 00.genomes/                  # Raw genome sequences (FASTA)
│   ├── 00.nomal_fnas/          # Original genome files
│   └── 01.Fusarium_302_masked_fnas/  # Soft-masked genomes
│
├── 01.RNAseqs/                  # RNA-seq processing
│   ├── 00.RNAseq_path.txt      # Sample ID list
│   ├── 02.run.sh               # Batch processing script
│   ├── 02.map_and_stringtie.sh # STAR mapping + StringTie assembly
│   └── work_*/                 # Working directories per species
│       ├── star_index/
│       ├── bams/
│       └── *_evidence_rnaseq.gtf
│
├── 02.Uniprot_AAs/              # Protein evidence
│   └── Hypocreales_filtered_*.fasta  # Filtered UniProt proteins
│
├── 03.Run_braker3/             # BRAKER3 annotation
│   ├── 02.run_CODE3_ROBUST.sh  # Main annotation script (ADAPTIVE mode)
│   └── braker_out_*/           # Output directories
│       ├── run_rna/            # BRAKER1 (RNA mode)
│       ├── run_prot/           # BRAKER2 (Protein mode)
│       └── final_annotation.gff3
│
└── 04.Reannotation_323genome_files/  # Final file conversion
    ├── 00.run.sh               # Batch conversion runner
    ├── 05.gtf2gff__cds_faa.sh  # Format conversion (GFF3 + CDS + Proteins)
    └── 00.fnas/ 01.gffs/      # Input/output for 323 genomes
```

## Detailed Steps

### Step 1: Genome Preprocessing (Repeat Masking)

```bash
# Usage: repeat_masking.sh <genome.fasta> <species_name>
bash 00.genomes/01.repeat_masking.sh input.fasta species_name
```

**Process:**
1. BuildRepeat database from genome
2. Predict repeat families with RepeatModeler (LTRStruct)
3. Soft-mask repeats with RepeatMasker (-xsmall flag)

**Output:** `*.fna.masked` file (repeats in lowercase)

### Step 2: RNA-seq Mapping and Assembly

```bash
# Usage: map_and_stringtie.sh <species_id> <genome.fa> <sample_list.txt>
bash 01.RNAseqs/02.map_and_stringtie-v2.sh species_id genome.fa sample_list.txt
```

**Process:**
1. Build STAR index
2. Align RNA-seq reads (MapQ >= 10 filter)
3. Merge all sample BAMs
4. Subsample to 20% forStringTie efficiency
5. Assemble with StringTie (min 200bp, coverage 3.0)

**Output:** `${species_id}_evidence_rnaseq.gtf` (hints file)

### Step 3: BRAKER3 Annotation

Three strategies are available (CODE3_ADAPTIVE is recommended):

#### Strategy A: Full Evidence Mode
```bash
# Usage: run_code1_full_evidence.sh <species> <genome> <proteins> <gtf> <bam_dir>
bash 99.test/01.run_code1_full_evidence.sh species_id genome.masked proteins.fa gtf bam_dir
```

#### Strategy B: BAM-based Mode
```bash
# See 99.test/02.run_code2_BAM-v1.sh
```

#### Strategy C: Adaptive Mode (Recommended)
```bash
# Usage: run_CODE3_ROBUST.sh <species> <masked_genome> <bam_dir> <protein_fa>
bash 03.Run_braker3/02.run_CODE3_ROBUST.sh species_id genome.masked bam_dir proteins.fa
```

**CODE3_ADAPTIVE Features:**
- Auto-detects genome fragmentation (scaffolds > 8000 → single thread mode)
- Smart BAM selection: tries exact match first, then best nearby relative (MapQ >= 30)
- Graceful fallback: RNA fails → protein-only; both fail → error
- TSEBRA merging when both RNA and protein evidence succeed

> **Note on fungal branch point model:** The `--fungus` flag (passed to GeneMark-EX for fungal-specific branch point model) has been enabled in both BRAKER1 (RNA) and BRAKER2 (Protein) calls. This improves splice site recognition and gene prediction accuracy in Fusarium genomes.

### Step 4: Output Conversion

```bash
# Convert GTF+GFF to standard format with sequences
# Usage: gtf2gff__cds_faa.sh <genome.fa> <input.gtf> <prefix>
bash 04.Reannotation_323genome_files/05.gtf2gff__cds_faa.sh genome.fa annotation.gtf species_name
```

**Outputs:**
- `${prefix}_standard.gff3` - Standard GFF3 format
- `${prefix}_cds.fa` - CDS nucleotide sequences
- `${prefix}_pep.fa` - Protein amino acid sequences

## Prerequisites

### Software Requirements (see environment.yml)
- **BRAKER3** (via Singularity container)
- **STAR** - RNA-seq aligner
- **StringTie** - Transcript assembler
- **SAMTools** - BAM processing
- **SeqKit** - Sequence manipulation
- **GFFRead** - Format conversion
- **RepeatModeler & RepeatMasker** - Repeat annotation

### System Requirements
- Linux (CentOS/Ubuntu)
- Singularity runtime
- 32+ CPU cores recommended
- 64GB+ RAM
- 500GB+ disk space

### Singularity Image
Download from: https://github.com/Gaius-Augustus/BRAKER

## Quick Start

### 1. Prepare Your Data
```
/path/to/
├── genomes/
│   └── species_name.fasta
├── rnaseq/
│   ├── sample1_1.fq.gz
│   └── sample1_2.fq.gz
└── proteins/
    └── reference_proteins.fasta
```

### 2. Run Full Pipeline
```bash
# Repeat masking
bash 00.genomes/01.repeat_masking.sh genome.fasta species_id

# RNA-seq mapping
bash 01.RNAseqs/02.map_and_stringtie-v2.sh species_id genome.fna data_list.txt

# BRAKER annotation
bash 03.Run_braker3/02.run_CODE3_ROBUST.sh species_id genome.masked bam_dir proteins.fa

# Generate final files
bash 04.Reannotation_323genome_files/05.gtf2gff__cds_faa.sh genome.fna final.gff3 output_prefix
```

### 3. Batch Processing
Edit `00.run.sh` in each step directory to add your species:

```bash
# In 04.Reannotation_323genome_files/00.run.sh
bash 05.gtf2gff__cds_faa.sh 00.fnas/YOUR_SPECIES.fna 01.gffs/YOUR_SPECIES.gtf YOUR_SPECIES
```

## Key Parameters

| Step | Parameter | Default | Description |
|------|------------|---------|-------------|
| STAR | `--outFilterMismatchNmax` | 15 | Max mismatches |
| STAR | `--alignIntronMax` | 50000 | Max intron length |
| StringTie | `-m` | 200 | Min transcript length |
| StringTie | `-c` | 3.0 | Min coverage |
| BRAKER | `--threads` | 32 | CPU threads |
| BRAKER | `--fungus` | enabled | Fungal branch point model (GeneMark-EX) |
| Filter | Min BAM reads | 100,000 | Minimum RNA evidence |

## Output Files

| File | Description |
|------|-------------|
| `final_annotation.gff3` | Complete gene annotation |
| `${species}_cds.fa` | Nucleotide sequences |
| `${species}_pep.fa` | Protein sequences |
| `braker.gff3` | Individual BRAKER outputs |
| `*.log` | Execution logs |

## Performance Notes

- **Runtime**: ~2-4 hours per genome (depending on genome size and data availability)
- **Memory**: Peak ~60GB during TSEBRA merging
- **Parallelization**: Each species processed independently

## Troubleshooting

### Common Issues

1. **GeneMark failure**: Check genome masking and scaffold count
2. **TSEBRA memory error**: Reduce thread count or subsample
3. **Missing hints**: Verify RNA-seq coverage and StringTie output
4. **Augustus species conflict**: Clear old species config before re-run

### Logs Location
- BRAKER: `braker_out_*/run_*/braker.log`
- STAR: `work_*/bams/*.star.log`
- Main script: Check stdout/stderr

## Citation

If you use this pipeline, please cite:
- BRAKER3: Bruna, T. et al. (2021)
- AUGUSTUS: Stanke, M. et al. (2006)
- TSEBRA: Gabriel, L. et al. (2021)

## License

MIT License

## Contact

For issues, please open an issue on the project repository.
