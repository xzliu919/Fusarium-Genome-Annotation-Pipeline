# Fusarium Genome Annotation Pipeline

This repository contains scripts for annotating 323 Fusarium genomes using BRAKER3.

## Pipeline Overview

1. **Repeat Masking** - `00.genomes/01.repeat_masking.sh`
2. **RNA-seq Processing** - `01.RNAseqs/` (STAR + StringTie)
3. **BRAKER3 Gene Prediction** - `03.Run_braker3/02.run_CODE3_ROBUST.sh`
4. **Post-processing** - `04.Reannotation_323genome_files/`

## BRAKER3 Parameters

The BRAKER3 pipeline runs both RNA-mode (BRAKER1, GeneMark-ET + Augustus) and Protein-mode (BRAKER2), then merges results via TSEBRA.

### ⚠️ Fungal Branch Point Model

**Not used in current implementation.** The `--fungus` flag (which instructs GeneMark-EX to use the fungal branch point model) was **not** passed to `braker.pl`. For fungal genome annotation, it is recommended to add this parameter:

```
--fungus
```

Add it to both `braker.pl` calls in `02.run_CODE3_ROBUST.sh`:

**BRAKER1 (RNA) call (line 171):**
```
/opt/BRAKER/scripts/braker.pl \
    --genome=/genome.fa \
    --bam=/rnaseq.bam \
    --workingdir=/output \
    --threads=${THREADS} \
    --species=${TARGET_ID}_rna \
    --gff3 \
    --softmasking \
    --fungus \
    --AUGUSTUS_CONFIG_PATH=$TARGET_CONFIG
```

**BRAKER2 (Protein) call (line 233):**
```
/opt/BRAKER/scripts/braker.pl \
    --genome=/genome.fa \
    --prot_seq=/proteins.fa \
    --workingdir=/output \
    --threads=${THREADS} \
    --species=${TARGET_ID}_prot \
    --gff3 \
    --softmasking \
    --fungus \
    --AUGUSTUS_CONFIG_PATH=$TARGET_CONFIG
```

The `--fungus` flag enables the fungal-specific branch point model in GeneMark-EX, which better recognizes splice sites in fungi and can improve gene prediction accuracy, particularly for intron boundary detection.
