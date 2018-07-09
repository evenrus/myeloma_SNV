# myeloma_snv

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]


## Description

Post-processing of targeted snv or indel data. Adapted to output-format
from leukgen/click_annotvcf. Custom panel-agnostic filtering with options
to include/exclude  genes and/or regions.  Optional annotation with
myeloma databases.

Development on-going.

### Features

Applies the following filters:  
    1. Pass by 1 or more callers  
    2. Target VAF > threshold (defalut 2 %)  
    3. Target depth > 100  
    4. Mutation supported by 1 or more reads on each strand  
    5. No SE or SR flag  
    6. Maximal population frequency < 0.005. If frequency 0.0025-0.005, remove missense mutations except in BRCA1, BRCA2 and TP53.  

Optional:  
    1. Filtering based on panel-specific good normals  
    2. Filtering based on gene lists and BED regions.  
    3. Annotation with myeloma datasets.  
      - 16 internal normals sequenced by myTYPE  
      - MMRF CoMMpass IA9 (889 WES, matched normal, not manually curated)  
      - Bolli Leukemia 2017 (418 targeted sequencing, manually curated)  
      - Lohr Cell 2014 203 WES/WGS  
      - myTYPE: Database of previously manually annotated variants  

### Output files

Variants passing all filters are annotated with optional datasets and written to csv ("goodcalls"). Failed variants are not annotated, but but also provided as output ("badcalls").

## Options

Required arguments:  

  --outdir:        Path to output directory.   
  --infile:        Path to input file with merged snv or indel calls in tsv.gz or
                   csv format.  

Optional arguments:  

  --normals:       Path to good normal calls in tsv.gz format  
  --vaf:           VAF threshold for filtering [default 0.02]  
  --genes:         Mutation frequencies by gene in excel-file. Colnames:
                   GENE and *freq*  
  --genes_bed:     BED-file of regions to keep.  
  --igh:           BED-file of regions to remove (e.g. IGH-region).  
  --mmrf:          MMRF reference file, tab separated text.  
  --bolli:         Bolli reference file, tab separated text.   
  --lohr:          Lohr reference file, tab separated text, hg19.  
  --mytype:        Manually annotated myTYPE data, csv file.    
  --version:       Show the version and exit.  
  --help:          Show this message and exit.  

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
