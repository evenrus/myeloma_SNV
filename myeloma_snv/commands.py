"""snv_process main command."""

import click

import pandas as pd
import numpy as np
import re
import os
import pybedtools as pyb 

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

def annotate_normals(snv, path_normals):
    normals = pd.read_csv(
        filepath_or_buffer=path_normals,
        compression='gzip',
        sep='\t',
        skiprows=1,
        low_memory=False)
    normals_counts=normals.groupby(["CHR","START","REF","ALT"]).size().reset_index(name="count")
    normals_counts=normals_counts[["CHR", "START","count"]].set_index(['CHR','START'])
    normals_median=normals.groupby(['CHR','START'])['TARGET_VAF'].median()
    normalC = []
    normalVAF = [] 
    for record in snv.itertuples(index=False, name=None):
        try:
            chrom=str(record[3])
            pos=int(record[4])
            tempC = normals_counts.loc[(chrom,pos),"count"]
            tempC = tempC.ix[0]
            normalC.append(int(tempC))
            normalVAF.append(str(normals_median.loc[(chrom,pos)]))
        except:
            normalC.append(0)
            normalVAF.append("0")
    snv["Normals_Frequency"] = normalC
    snv["Normals_median_VAF"] = normalVAF
    return(snv)
    
def annotate_bolli(snv, path_bolli):
    bolli = pd.read_csv(filepath_or_buffer=path_bolli, sep='\t')
    bolli = bolli[["CHR", "START", "WT", "MT", "Variant_class"]]
    bolli_counts = bolli.groupby(['CHR','START'])['MT'].count()
    bolli_var = bolli.drop(['WT','MT'], axis = 1) 
    bolli_var = bolli_var.set_index(['CHR', 'START'])
    cl = [] 
    freq = [] 
    positions = []
    annot = []
    for record in snv.itertuples(index=False, name=None):
        flag = 0
        try:
            chrom = str(record[3])
            pos = int(record[4])
            start = int(record[4]) - 9
            end = int(record[4]) + 9
            if (chrom, pos) in  bolli_counts.index:
                cl.append("genomic_exact")
                freq.append(str(bolli_counts.loc[(chrom,pos)]))
                positions.append(str(pos))
                annot.append(str(bolli_var.loc[chrom, pos]['Variant_class'].values[0]))
                flag = 1
            if flag == 0: 
                bolli_counts_sub=bolli_counts.loc[chrom]
                if not bolli_counts_sub[(bolli_counts_sub.index >= start) & (bolli_counts_sub.index <= end)].empty:
                    fr = []
                    posit = []
                    ann = []
                    for i in bolli_counts_sub[(bolli_counts_sub.index >= start) & (bolli_counts_sub.index <= end)].index.values:
                        fr.append(str(bolli_counts.loc[(chrom,i)]))
                        posit.append(str(i))
                        ann.append(str(bolli_var.loc[(chrom,i)]['Variant_class'].values[0]))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    positions.append((":".join(posit)))
                    annot.append((":".join(ann)))
                else:
                    cl.append(None)
                    freq.append(None)
                    positions.append(None)
                    annot.append(None)
        except:
            cl.append(None)
            freq.append(None)
            positions.append(None)
            annot.append(None)
    snv["Bolli_Class"] = cl
    snv["Bolli_Frequency"] = freq
    snv["Bolli_Positions"] = positions
    snv["Bolli_Annotation"] = annot
    return(snv)

def annotate_mmrf(snv, path_mmrf):
    mmrf = pd.read_csv(filepath_or_buffer=path_mmrf, sep='\t')
    mmrf=mmrf[["Sample", "CHROM", "POS", "REF", "ALT", "GEN[0].AR", "GEN[1].AR"]]
    mmrf=mmrf.drop_duplicates()
    mmrfM=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].median()
    mmrfC=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].count()
    mmrfQ25=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].quantile(q=0.25)
    mmrfQ75=mmrf.groupby(['CHROM','POS'])['GEN[1].AR'].quantile(q=0.75)
    cl = [] 
    freq = [] 
    medVAF = [] 
    Q25 = [] 
    Q75 = [] 
    positions = [] 
    for record in snv.itertuples(index=False, name=None):
        flag = 0
        try:    #what does try/except pairs do?
            # Define position and chrom variable positions outside of loop based on varnames.
            chrom = str(record[3])
            pos = int(record[4])
            start = int(record[4]) - 9
            end = int(record[4]) + 9
            # 
            if (chrom, pos) in mmrfC.index:
                cl.append("genomic_exact")
                freq.append(str(mmrfC.loc[(chrom,pos)]))
                medVAF.append(str(mmrfM.loc[(chrom,pos)]))
                Q25.append(str(mmrfQ25.loc[(chrom,pos)]))
                Q75.append(str(mmrfQ75.loc[(chrom,pos)]))
                positions.append(str(pos))
                flag = 1
            if flag == 0:
                mmrfCsub=mmrfC.loc[chrom]
                if not mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].empty:
                    fr = []
                    mv = []
                    Q2 = []
                    Q7 = []
                    posit = []
                    for i in mmrfCsub[(mmrfCsub.index >= start) & (mmrfCsub.index <= end)].index.values:
                        fr.append(str(mmrfC.loc[(chrom,i)]))
                        mv.append(str(mmrfM.loc[(chrom,i)]))
                        Q2.append(str(mmrfQ25.loc[(chrom,i)]))
                        Q7.append(str(mmrfQ75.loc[(chrom,i)]))
                        posit.append(str(i))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    medVAF.append((":".join(mv)))
                    Q25.append((":".join(Q2)))
                    Q75.append((":".join(Q7)))
                    positions.append((":".join(posit)))
                else:
                    cl.append(None)
                    freq.append(None)
                    medVAF.append(None)
                    Q25.append(None)
                    Q75.append(None)
                    positions.append(None)
        except:
            cl.append(None)
            freq.append(None)
            medVAF.append(None)
            Q25.append(None)
            Q75.append(None)
            positions.append(None)
    snv["MMRF_Class"] = cl
    snv["MMRF_Frequency"] = freq
    snv["MMRF_VAF"] = medVAF
    snv["MMRF_Q25"] = Q25
    snv["MMRF_Q75"] = Q75
    snv["MMRF_Positions"] = positions
    return(snv)

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
    snv=pd.merge(snv, freqlist, how='left')
    return(snv)

def annotate_maf(snv):
    #adds column with maximal MAF of variant in any normal database
    snv['MAX_MAF']=0 #sets variable to 0 if frequency is not reported
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

def filter_IGH(snv, IGH_path):
    #Remove calls in IGH locus
    snv_bed = snv[["CHR", "START", "END", "ID_VARIANT"]]
    snv_bed = pyb.BedTool.from_dataframe(snv_bed)
    igh = pyb.BedTool(IGH_path)
    snv_inter = snv_bed.intersect(igh, u=True)
    snv_inter = pyb.BedTool.to_dataframe(snv_inter)['name']
    snv['MFLAG_IGH'] = np.where(snv.ID_VARIANT.isin(snv_inter), 1, 0)
    return(snv)

def filter_MAF(snv):
    snv['MFLAG_MAF'] = np.where(snv['MAX_MAF'] <= 0.03, 0, 1)
    return(snv)

def filter_MAF_COSMIC(snv):
    #Remove if >0.1 % MAF and not in COSMIC
    snv['MFLAG_MAFCOS'] = np.where((snv['MAX_MAF'] > 0.001) & (snv['COSMIC'].isnull()), 1, 0)
    return(snv)

def filter_nonpass(snv):
    #Remove non-PASS calls not present in COSMIC or previous cohorts (MMRF, Bolli, etc.)
    snv['MFLAG_NONPASS'] = np.where((snv['FILTER'] != "PASS") & (snv['COSMIC'].isnull()) & (snv['MMRF_Class'].isnull()) & (snv['Bolli_Class'].isnull()), 1, 0)
    return(snv)

def filter_normals(snv):
    #Remove calls present in at least 4 internal normals
    snv['MFLAG_NORM'] = np.where((snv['Normals_Frequency'] >= 4), 1, 0)
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
        f.write(f'MFLAG_PANEL: Gene not in panel: {snv["MFLAG_PANEL"].sum()}\n')
        f.write(f'MFLAG_IGH: In IGH locus: {snv["MFLAG_IGH"].sum()}\n')
        f.write(f'MFLAG_MAF: MAF > 3 % in exax/1000genomes: {snv["MFLAG_MAF"].sum()}\n')
        f.write(f'MFLAG_MAFCOS: MAF > 0.1 % and not in COSMIC: {snv["MFLAG_MAFCOS"].sum()}\n')
        f.write(f'MFLAG_NONPASS: NON-PASS and not in COSMIC, MMRF or Bolli: {snv["MFLAG_NONPASS"].sum()}\n')
        f.write(f'MFLAG_NORM: Variant in >= 4 good normals: {snv["MFLAG_NORM"].sum()}\n')
        f.write(f'Removing calls with >= 1 MFLAG: {bad.shape[0]}\n')
        f.write(f'Calls passed filters: {good.shape[0]}\n')
    return()

# Main Function
def process(
        infile,
        skiplines,
        outdir,
        genes,
        igh,
        mmrf,
        bolli,
        lohr,
        normals):
    """Main function to process myTYPE SNV output"""
    ##IMPORTING DATA
    snv = import_snv(infile, skiplines)

    ##ANNOTATIONS
    snv = annotate_COSMIC(snv) 
    snv = annotate_genefreq(snv, genes) #Replace this with mutation frequency from MMRF? (and other raw data?)
    snv = annotate_maf(snv)
    snv = annotate_mmrf(snv, mmrf)
    snv = annotate_bolli(snv, bolli)
    snv = annotate_normals(snv, normals)
    #snv = annotate_lohr(snv, lohr) #Waiting for wile from Teja.

    ##FILTERS
    snv = filter_panel(snv, genes) #Remove calls outside of myTYPE panel
    snv = filter_IGH(snv, igh) #Remove calls in IGH locus
    snv = filter_MAF(snv) #Remove if MAF > 3 %.
    snv = filter_MAF_COSMIC(snv) #Remove if >0.1 % MAF and not in COSMIC
    snv = filter_nonpass(snv) #Remove non-PASS calls not present in COSMIC or previous cohorts (MMRF, Bolli, etc.)
    snv = filter_normals(snv) #Remove calls present in at least 4 internal normals
    #snv = filter_synonymous(snv) ## Already done with input - Not necessary?

    ##OUTPUT
    name = namecore(infile)
    filter_export(snv, outdir, name)
    print('SNV processing complete')