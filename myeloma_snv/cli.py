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

"""
Post-processing script for SNV and indel data from myTYPE
with myeloma-specific annotations and filters
Version 1.0
"""

import click

from myeloma_snv import __version__
from myeloma_snv import commands

@click.command()
@click.option(
    "--outdir",
    required=True,
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    help="Path to output directory.")
@click.option(
    "--infile",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to input file with merged SNV calls in tsv.gz or csv format")
@click.option(
    "--normals",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to good normal calls csv file")
@click.option(
    "--vaf",
    required=False,
    default=0.02,
    type=float,
    help="VAF threshold for filtering (default 0.02)")
@click.option(
    "--genes",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Mutation frequencies by gene in excel-file. Colnames: GENE and *freq*")
@click.option(
    "--genes_drop",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Name of genes to filter out in excel-file. Colname: GENE")
@click.option(
    "--genes_bed",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="BED-file of regions to keep.")
@click.option(
    "--igh",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="BED-file of regions to remove (e.g. IGH-region).")
@click.option(
    "--mmrf",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help=" MMRF reference file, tab separated text")
@click.option(
    "--bolli",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Bolli reference file, tab separated text")
@click.option(
    "--lohr",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Lohr reference file, tab separated text, hg19.")
@click.option(
    "--mytype",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Manually annotated myTYPE data, csv file")
@click.version_option(__version__)

def main(outdir, infile, normals, vaf, genes, genes_drop, genes_bed,
         igh, mmrf, bolli, lohr, mytype):
    r"""
    Post-processing of targeted snv or indel data.
    Adapted to output format from leukgen/click_annotvcf.
    Custom panel-agnostic filtering with options to include/exclude 
    genes and/or regions. 
    Optional annotation with myeloma databases. 

    """
    commands.process(
        infile=infile,
        outdir=outdir,
        normals=normals,
        vaf=vaf,
        genes=genes,
        genes_drop=genes_drop,
        genes_bed=genes_bed,
        igh=igh,
        mmrf=mmrf,
        bolli=bolli,
        lohr=lohr,
        mytype=mytype
        )
