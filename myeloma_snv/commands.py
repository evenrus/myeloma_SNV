"""snv_process main command."""

import click

import pandas as pd
import numpy as np

def annotate_maf(snv):
    #adds column with maximal MAF of variant in any normal database
    snv['MAX_MAF']=snv.filter(regex='MAF').max(axis=1)
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

def filter_MAF_COSMIC(snv):
    snv_out=snv.drop(snv[(snv['MAX_MAF'] > 0.001) & (snv['COSMIC'].isnull())].index)
    print('Removed %(removed)d calls with MAF > 0.1 %% in EXAC or 1000 genomes that are not present in COSMIC. Remaining %(remaining)d\n' %
            {'removed' : snv.shape[0]-snv_out.shape[0], 'remaining' : snv_out.shape[0]})
    return(snv_out)

# Main Function
def process(
        infile,
        skiplines,
        outfile,
        genes):
    """Main function to process myTYPE SNV output"""
    ##IMPORTING DATA
    snv = pd.read_csv(
        filepath_or_buffer=infile,
        skiprows=skiplines)
    print('Loaded file containing %d SNV calls. Annotating...\n' % snv.shape[0])

    ##FILTERS
    #snv = filter_IGH(snv) awaiting function
    #snv = filter_synonymous(snv) awaiting function

    ##ANNOTATIONS
    #snv = annotate_COSMIC(snv) - make column HEME_EXACT
    snv = annotate_maf(snv)
    print('####\nAnnotation complete, running filters...\n')

    ##FILTERS
    snv = filter_panel(snv, genes) #Remove calls outside of myTYPE panel
    snv = filter_MAF(snv) #Remove if MAF > 3 %.
    snv = filter_MAF_COSMIC(snv) #Remove if >0.1 % MAF and not in COSMIC
    #snv <- filter_nonpass(snv) awaiting function
    #snv <- filter_normals(snv) awaiting function

    snv.to_csv(
        path_or_buf=outfile,
        index=False)
    print('####\nsnv post-processing complete. Output stored as %s' % outfile)
