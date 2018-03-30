"""snv_process main command."""

import click

import pandas as pd
import numpy as np
import re

def import_snv(path, skiplines):
    #determine filetype and import, returns pandas dataFrame
    if re.search('.csv$', path):
        snv = pd.read_csv(
            filepath_or_buffer=path,
            skiprows=skiplines)
        try:
            snv
        except NameError:
            raise Exception(f'Error when importing file {path}')
        else:
            print(f'Loaded csv-file containing {snv.shape[0]} SNV calls. Annotating...\n')
            return(snv)
    elif re.search('.tsv.gz$', path):
        snv = pd.read_csv(
            filepath_or_buffer=path,
            compression='gzip',
            sep='\t',
            skiprows=skiplines)
        try:
            snv
        except NameError:
            raise Exception(f'Error when importing file {path}')
        else:
            print(f'Loaded csv-file containing {snv.shape[0]} SNV calls. Annotating...\n')
            return(snv)
    else:
        raise Exception(f'Input file {path} has unsupported extension: try .csv or .tsv.gz')

def annotate_maf(snv):
    #adds column with maximal MAF of variant in any normal database
    snv['MAX_MAF']=snv.filter(regex='MAF').max(axis=1)
    return(snv)

def annotate_lohr(snv, lohr):
    #args snv = snv file; lohr = raw data from lohr 2014 in hg19 format.
    return(snv) 

def filter_panel(snv, genes):
    #args: snv = snv file; genes = path to file with genes to include
    panel = pd.read_excel(io=genes)
    keep = panel['GENE']
    snv_out = snv[snv.GENE.isin(keep)]
    print('Removed %(removed)d calls in genes outside of myTYPE panel. Remaining: %(remaining)d\n' % 
            {'removed' : snv.shape[0]-snv_out.shape[0], 'remaining' : snv_out.shape[0]})
    return(snv_out)

def filter_MAF(snv):
    snv_out=snv.loc[snv['MAX_MAF'] <= 0.03]
    print('Removed %(removed)d calls with MAF > 3 %% in EXAC or 1000 genomes. Remaining: %(remaining)d\n' %
            {'removed' : snv.shape[0]-snv_out.shape[0], 'remaining' : snv_out.shape[0]})
    return(snv_out)

    #print(f'This is a {snv} blabla {snv_out}') -- new way of doing the %% thing. 

def filter_MAF_COSMIC(snv):
    snv_out=snv.drop(snv[(snv['MAX_MAF'] > 0.001) & (snv['COSMIC'].isnull())].index)
    print('Removed %(removed)d calls with MAF > 0.1 %% in EXAC or 1000 genomes that are not present in COSMIC. Remaining %(remaining)d\n' %
            {'removed' : snv.shape[0]-snv_out.shape[0], 'remaining' : snv_out.shape[0]})
    return(snv_out)

def namecore(infile):
    name = infile.split('/')[-1]
    if re.search('.csv$', name):
        return(re.sub('.csv$', '', name))
    else:
        return(re.sub('.tsv.gz$', '', name))

def filter_export(snv, outdir, name):
    good=snv
    bad=snv
    if not re.search('/$', outdir):
        outdir =''.join([outdir,'/'])
    goodname=''.join([outdir, name, '_goodcalls.csv'])
    badname=''.join([outdir, name, '_badcalls.csv'])
    good.to_csv(
        path_or_buf=goodname,
        index=False)
    bad.to_csv(
        path_or_buf=badname,
        index=False)
    print('####\nSNV processing complete')
    return()

# Main Function
def process(
        infile,
        skiplines,
        outdir,
        genes,
        lohr):
    """Main function to process myTYPE SNV output"""
    ##IMPORTING DATA
    snv = import_snv(infile, skiplines)

    ##ANNOTATIONS
    #snv = annotate_COSMIC(snv) - make column HEME_EXACT
    snv = annotate_maf(snv)
    #snv = annotate_lohr(snv, lohr)
    print('####\nAnnotation complete, running filters...\n')

    ##FILTERS
    snv = filter_panel(snv, genes) #Remove calls outside of myTYPE panel
    snv = filter_MAF(snv) #Remove if MAF > 3 %.
    snv = filter_MAF_COSMIC(snv) #Remove if >0.1 % MAF and not in COSMIC
    #snv <- filter_nonpass(snv) awaiting function
    #snv <- filter_normals(snv) awaiting function
    #snv = filter_IGH(snv) awaiting function
    #snv = filter_synonymous(snv) awaiting function

    ##OUTPUT
    name = namecore(infile)
    filter_export(snv, outdir, name)
