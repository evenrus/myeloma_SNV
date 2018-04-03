"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will
cause problems, the code will get executed twice:

    - When you run `python -m myeloma_snv` python will execute
      `__main__.py` as a script. That means there won't be any
      `myeloma_snv.__main__` in `sys.modules`.

    - When you import __main__ it will get executed again (as a module) because
      there's no `myeloma_snv.__main__` in `sys.modules`.

Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

# Program for post-processing of CNV data from myTYPE with myeloma-specific annotations and filters
# v 29.3.18
# Current version takes input file already subjected to some processing:
## - Annotated with COSMIC, 1000g and EXAC.
## - Removed synonymous (AND IGH locus?)

import click

from myeloma_snv import __version__
from myeloma_snv import commands

@click.command()
@click.option(
    "--outdir",
    required=True,
    type=click.Path(exists=True),
    help="Path to output directory.")
@click.option(
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="Path to input file with merged SNV calls in tsv.gz or csv format")
@click.option(
    "--skiplines",
    default=0,
    show_default=True,
    help="Number of lines to skip at the top of input file when importing.")
@click.option(
    "--genes",
    default='references/myType_genes.xlsx',
    show_default=True,
    type=click.Path(exists=True),
    help="Excel file with column 'GENES'. Used to filter out variants in other genes.")
@click.option(
    "--mmrf",
    default='references/MMRF_CoMMpass_IA9_All_Canonical_Variants.txt',
    show_default=True,
    type=click.Path(exists=True),
    help="Path to MMRF reference file in .tsv format")
@click.option(
    "--lohr",
    default='references/lohr_variants.xlsx',
    show_default=True,
    type=click.Path(exists=True),
    help="Excel file with raw data from lohr et al 2014 - changed from hg18 to hg19 format.")
@click.version_option(__version__)

def main(outdir, infile, skiplines, genes, mmrf, lohr):
    r"""
    rogram for post-processing of CNV data from myTYPE with myeloma-specific annotations and filters
    """
    commands.process(
        infile=infile,
        skiplines=skiplines,
        outdir=outdir,
        genes=genes,
        mmrf=mmrf,
        lohr=lohr
        )

#### ANNOTATIONS OF INPUT FILE ####
## Another general purpose post-processing script annotates with:
## COSMIC, EXAC, 1000 genomes, VAF and other variant metadata


#### ANNOTATIONS IN THIS SCRIPT ####
#16 internal normals sequenced by myTYPE
#MMRF (889 WES, matched normal, not manually curated)
#Bolli (418 targeted seq, manually curated)

#Lohr ~200 WES/WGS


#### FILTERS ####
#Remove:
# -calls in IGH locus
# -synonymous variants
##NEW## -calls in genes not in myTYPE panel
##DONE## -calls with > 3 % MAF in Exac or 1000 genomes
##DONE## -calls with >0.1 % MAF in Exac or 1000 genomes if not present in COSMIC
# -non-PASS calls not present in COSMIC or previous cohorts (MMRF, Bolli, etc.)
# -calls present in at least 4 internal normals
# -non-pass calls present in at least 1 normal AND not in COSMIC, MMRF, etc. ##Redundant?
