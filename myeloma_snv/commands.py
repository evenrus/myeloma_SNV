"""snv_process main command."""

import click

import pandas as pd
import numpy as np
import re
import os

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
            print(f'Loaded csv-file containing {snv.shape[0]} SNV calls. Processing...')
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
            print(f'Loaded csv-file containing {snv.shape[0]} SNV calls. Processing...')
            return(snv)
    else:
        raise Exception(f'Input file {path} has unsupported extension: try .csv or .tsv.gz')

def annotate_COSMIC(snv):
    #make column HEME_EXACT
    heme_exact = []
    cosmic = snv['COSMIC'].tolist()
    for entry in cosmic:
        if pd.isnull(entry):
            heme_exact.append(None)
        else:
            first=entry.split('|')[0]
            if re.search('^GENOMIC_EXACT', first):
                if re.search('HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', first):
                    count=re.search('(?<=HAEMATOPOIETIC_AND_LYMPHOID_TISSUE=)\w+', first)[0]
                    heme_exact.append(count)
                else:
                    heme_exact.append(None)
            else:
                heme_exact.append(None)
    snv['HEME_EXACT']=heme_exact
    return(snv)

def annotate_genefreq(snv, genes):
    #adds column with maximal mutation frequency in gene as previously published in large MM studies.
    freqlist = pd.read_excel(io=genes)
    freqlist['MAX_MUTFREQ']=round(freqlist.filter(regex='freq').max(axis=1),1)
    freqlist=freqlist[['GENE', 'MAX_MUTFREQ']]
    snv=pd.merge(snv, freqlist)
    return(snv)

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
    snv['MFLAG_PANEL'] = np.where(snv.GENE.isin(keep), 0, 1)
    return(snv)

def filter_MAF(snv):
    snv['MFLAG_MAF'] = np.where(snv['MAX_MAF'] <= 0.03, 0, 1)
    return(snv)

def filter_MAF_COSMIC(snv):
    snv['MFLAG_MAFCOS'] = np.where((snv['MAX_MAF'] > 0.001) & (snv['COSMIC'].isnull()), 1, 0)
    return(snv)

def namecore(infile):
    name = infile.split('/')[-1]
    if re.search('.csv$', name):
        return(re.sub('.csv$', '', name))
    else:
        return(re.sub('.tsv.gz$', '', name))

def filter_export(snv, outdir, name):
    good=snv[snv.filter(regex='MFLAG').sum(axis=1) == 0]
    bad=snv[snv.filter(regex='MFLAG').sum(axis=1) > 0]
    if not re.search('/$', outdir):
        outdir =''.join([outdir,'/'])
    goodname=''.join([outdir, name, '_goodcalls.csv'])
    badname=''.join([outdir, name, '_badcalls.csv'])
    textname=''.join([outdir, name, '_report.txt'])
    good.to_csv(
        path_or_buf=goodname,
        index=False)
    bad.to_csv(
        path_or_buf=badname,
        index=False)
    with open(textname, 'w') as f:
        f.write(f'Imported SNV calls: {snv.shape[0]}\n')
        f.write('Flagging variants for filtering:\n')
        f.write(f'MFLAG: Gene not in panel: {snv["MFLAG_PANEL"].sum()}\n')
        f.write(f'MFLAG: MAF > 3 % in exax/1000genomes: {snv["MFLAG_MAF"].sum()}\n')
        f.write(f'MFLAG: MAF > 0.1 % and not in COSMIC: {snv["MFLAG_MAFCOS"].sum()}\n')
        f.write(f'Calls filtered out: {bad.shape[0]}\n')
        f.write(f'Calls passed filters: {good.shape[0]}\n')
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
    snv = annotate_COSMIC(snv) 
    snv = annotate_genefreq(snv, genes)
    snv = annotate_maf(snv)
    #snv = annotate_lohr(snv, lohr)

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
    print('SNV processing complete')
