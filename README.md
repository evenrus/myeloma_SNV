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
  --version            Show the version and exit.
  --help               Show this message and exit.

## Description
Takes as input a file containing indel or SNV calls in .csv or .tsv.gz format.

Mandatory input columns:
COSMIC, EXAC, 1000 genomes, VAF and other variant metadata

Annotation columns:
16 internal normals sequenced by myTYPE
MMRF (889 WES, matched normal, not manually curated)
Bolli (418 targeted seq, manually curated)
Lohr ~200 WES/WGS

Filtering criteria:
 -calls in IGH locus
 -calls in genes not in myTYPE panel
 -calls with > 3 % MAF in Exac or 1000 genomes
 -calls with >0.1 % MAF in Exac or 1000 genomes if not present in COSMIC
 -non-PASS calls not present in COSMIC or previous cohorts (MMRF, Bolli, etc.)
 -calls present in at least 1 internal normals
 -calls with VAF < 1 %

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
