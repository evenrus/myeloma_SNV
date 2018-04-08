# myeloma_snv

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]

Program for processing SNV or indel calls from myTYPE with myeloma-specific annotations and filters.

Development on-going.

## Options

myeloma_SNV:

  --mode [snv|indel]   Set input variant type: snv or indel  [required]
  --outdir PATH        Path to output directory.  [required]
  --infile PATH        Path to input file with merged SNV calls in tsv.gz or
                       csv format  [required]
  --skiplines INTEGER  Number of lines to skip at the top of input file when
                       importing.  [default: 0]
  --genes PATH         Excel file with column 'GENES'. Used to filter out
                       variants in other genes.  [required]
  --genes_bed PATH     Bed file containing panel regions, to filter out
                       outside calls.  [required]
  --igh PATH           BED file with IGH locus to filter out overlapping
                       calls.  [required]
  --mmrf PATH          Path to MMRF reference file, tab separated text
                       [required]
  --bolli PATH         Path to Bolli reference file, tab separated text
                       [required]
  --lohr PATH          Path to Lohr reference file - tab separated hg19
                       format.  [required]
  --normals PATH       Path to good normal calls in tsv.gz format  [required]
  --mytype PATH        Path to manually annotated myTYPE data in csv format
  --version            Show the version and exit.
  --help               Show this message and exit.

## Description
Takes as input a file containing indel or SNV calls in .csv or .tsv.gz format.

###Mandatory columns in input file
ID_VARIANT = Unique identifier of each variant
CHR = chromosome (int)
START = variant start position (int)
STOP = variant stop position (int)
GENE = Gene name (str)
TARGET_VAF = Variant allele frequency of variant (num)
EFFECT = predicted effect of variant, e.g. 'non_synonymous_codon'
BIDIR = 1 if reads supporting variant is found on both strands, 0 otherwise
FILTER = 'PASS' if variant has passed all previously applied variant caller filters.
COSMIC = Cosmic annotation
*_MAF = Frequency of mutation in EXAC populations and 1000 genomes database.

###Annotations:
Data from the following datasets are incorporated as annotations for each variant:
-16 internal normals sequenced by myTYPE
-MMRF CoMMpass IA9 (889 WES, matched normal, not manually curated)
-Bolli Leukemia 2017 (418 targeted sequencing, manually curated)
-Lohr Cell 2014 203 WES/WGS
-myTYPE: Database of previously manually annotated variants 

###Filtering by myeloma (M) FLAGs:
-MFLAG_PANEL: Gene not in panel
-MFLAG_IGH: In IGH locus
-MFLAG_MAF: MAF > 3 % in exac/1000genomes
-MFLAG_MAFCOS: MAF > 0.1 % and not in COSMIC
                For SNVs: Only EXACT/POS in COSMIC counts as match.
-MFLAG_NONPASS: NON-PASS IF not in COSMIC and not previously known in MM.
                For SNVs: Only EXACT/POS in COSMIC counts as match.
                          Only missense mutations can be removed by this filter (i.e. 'non_synonymous_codon')
-MFLAG_NORM: Variant in 1 or more good normal control run by myTYPE
-MFLAG_VAF: Remove variants with target VAF < 1 %
-MFLAG_BIDIR: Remove variants BIDIR = 0 (only reads on one strand)

###Output files
Variants with no MFLAGs are output into a file for 'good calls'. Other variants are output as 'bad calls'. A run summary report including flag statistics is created into the same folder. 

## Credits

This package was created using [Cookiecutter] and the
[leukgen/cookiecutter-toil] project template.

<!-- References -->
[singularity]: http://singularity.lbl.gov/
[docker2singularity]: https://github.com/singularityware/docker2singularity
[cookiecutter]: https://github.com/audreyr/cookiecutter
[leukgen/cookiecutter-toil]: https://github.com/leukgen/cookiecutter-toil
[`--batchSystem`]: http://toil.readthedocs.io/en/latest/developingWorkflows/batchSystem.html?highlight=BatchSystem

<!-- Badges -->
[docker_base]: https://hub.docker.com/r/evenrus/myeloma_snv
[docker_badge]: https://img.shields.io/docker/build/evenrus/myeloma_snv.svg
[automated_badge]: https://img.shields.io/docker/automated/leukgen/myeloma_snv.svg
[codecov_badge]: https://codecov.io/gh/evenrus/myeloma_snv/branch/master/graph/badge.svg
[codecov_base]: https://codecov.io/gh/evenrus/myeloma_snv
[pypi_badge]: https://img.shields.io/pypi/v/myeloma_snv.svg
[pypi_base]: https://pypi.python.org/pypi/myeloma_snv
[travis_badge]: https://img.shields.io/travis/evenrus/myeloma_snv.svg
[travis_base]: https://travis-ci.org/evenrus/myeloma_snv
