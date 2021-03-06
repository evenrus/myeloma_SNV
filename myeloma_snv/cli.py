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
    "--mode",
    required=True,
    type=click.Choice(["snv", "indel"]),
    help="Set input variant type: snv or indel")
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
    "--genes",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Mutation frequencies by gene in excel-file. Cnames: GENE and *freq*")
@click.option(
    "--genes_drop",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="List of genes to filter out. Cname: GENE")
@click.option(
    "--genes_bed",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Bed file genes to keep.")
@click.option(
    "--igh",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="BED file with IGH locus to filter out overlapping calls.")
@click.option(
    "--mmrf",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to MMRF reference file, tab separated text")
@click.option(
    "--bolli",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to Bolli reference file, tab separated text")
@click.option(
    "--lohr",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to Lohr reference file - tab separated hg19 format.")
@click.option(
    "--normals",
    required=True,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to good normal calls in tsv.gz format")
@click.option(
    "--mytype",
    required=False,
    type=click.Path(file_okay=True, readable=True, resolve_path=True),
    help="Path to manually annotated myTYPE data in csv format")
@click.version_option(__version__)

def main(mode, outdir, infile, genes, genes_drop, genes_bed,
         igh, mmrf, bolli, lohr, normals, mytype):
    r"""
    Post-processing script for SNV and indel data from myTYPE
    with myeloma-specific annotations and filters
    """
    commands.process(
        mode=mode,
        infile=infile,
        outdir=outdir,
        genes=genes,
        genes_drop=genes_drop,
        genes_bed=genes_bed,
        igh=igh,
        mmrf=mmrf,
        bolli=bolli,
        lohr=lohr,
        normals=normals,
        mytype=mytype
        )
