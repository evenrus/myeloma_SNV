"""variants_process main command."""

import click

import pandas as pd
import numpy as np
import re
import os
import pybedtools as pyb
from datetime import datetime

run_time = str(datetime.today())

## IMPORT VARIANTS FILE
def import_variants(path, skiplines):
    #determine filetype and import, returns pandas dataFrame
    if re.search('.csv$', path):
        variants = pd.read_csv(
            filepath_or_buffer=path,
            skiprows=skiplines,
            low_memory=False)
        try:
            variants
        except NameError:
            raise Exception(f'Error when importing file {path}')
        else:
            print(f'Loaded csv-file containing {variants.shape[0]} variant calls. Processing...')
            return(variants)
    elif re.search('.tsv.gz$', path):
        variants = pd.read_csv(
            filepath_or_buffer=path,
            compression='gzip',
            sep='\t',
            skiprows=skiplines,
            low_memory=False)
        try:
            variants
        except NameError:
            raise Exception(f'Error when importing file {path}')
        else:
            print(f'Loaded tsv.gz-file containing {variants.shape[0]} variant calls. Processing...')
            return(variants)
    else:
        raise Exception(f'Input file {path} has unsupported extension: try .csv or .tsv.gz')

## ANNOTATION FUNCTIONS
def annotate_COSMIC(variants):
    #make column HEME_EXACT and ANY_EXACT_POS
    heme_exact = []
    cosmic = variants['COSMIC'].tolist()
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
    variants['HEME_EXACT']=heme_exact
    any_exact_pos = []
    for entry in cosmic:
        if pd.isnull(entry):
            any_exact_pos.append(0)
        elif re.search('GENOMIC_EXACT', entry) or re.search('GENOMIC_POS', entry):
            any_exact_pos.append(1)
        else:
            any_exact_pos.append(0)
    variants['ANY_EXACT_POS']=any_exact_pos
    return(variants)

def annotate_genefreq(variants, genes):
    #adds column with maximal mutation frequency in gene as previously published in large MM studies.
    freqlist = pd.read_excel(io=genes)
    freqlist['MAX_MUTFREQ']=round(freqlist.filter(regex='freq').max(axis=1),1)
    freqlist=freqlist[['GENE', 'MAX_MUTFREQ']]
    variants=pd.merge(variants, freqlist, how='left')
    return(variants)

def annotate_maf(variants):
    #adds column with maximal MAF of variant in any normal database
    variants['MAX_MAF']=0 #sets variable to 0 if frequency is not reported
    variants['MAX_MAF']=variants.filter(regex='MAF').max(axis=1)
    return(variants)

def annotate_normals(variants, path_normals):
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
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START") 
    for record in variants.itertuples(index=False, name=None):
        try:
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            tempC = normals_counts.loc[(chrom,pos),"count"]
            tempC = tempC.ix[0]
            normalC.append(int(tempC))
            normalVAF.append(str(normals_median.loc[(chrom,pos)]))
        except:
            normalC.append(0)
            normalVAF.append("0")
    variants["Normals_Frequency"] = normalC
    variants["Normals_median_VAF"] = normalVAF
    return(variants)

def annotate_mmrf(variants, path_mmrf):
    mmrf = pd.read_csv(filepath_or_buffer=path_mmrf, sep='\t')
    mmrf=mmrf[["Sample", "CHROM", "POS", "REF", "ALT", "GEN[0].AR", "GEN[1].AR"]]
    mmrf=mmrf.drop_duplicates() ## What are these duplicates?
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
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START")
    for record in variants.itertuples(index=False, name=None):
        flag = 0
        try:   
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            start = int(record[poscol]) - 9
            end = int(record[poscol]) + 9
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
    variants["MMRF_Class"] = cl
    variants["MMRF_Frequency"] = freq
    variants["MMRF_VAF"] = medVAF
    variants["MMRF_Q25"] = Q25
    variants["MMRF_Q75"] = Q75
    variants["MMRF_Positions"] = positions
    return(variants)

def annotate_bolli(variants, path_bolli):
    bolli = pd.read_csv(filepath_or_buffer=path_bolli, sep='\t')
    bolli = bolli[["CHR", "START", "WT", "MT", "Variant_class"]]
    bolli_counts = bolli.groupby(['CHR','START'])['MT'].count()
    bolli_var = bolli.drop(['WT','MT'], axis = 1) 
    bolli_var = bolli_var.set_index(['CHR', 'START'])
    cl = [] 
    freq = [] 
    positions = []
    annot = []
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START")
    for record in variants.itertuples(index=False, name=None):
        flag = 0
        try:
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            start = int(record[poscol]) - 9
            end = int(record[poscol]) + 9
            if (chrom, pos) in  bolli_counts.index:
                cl.append("genomic_exact")
                freq.append(str(bolli_counts.loc[(chrom,pos)]))
                positions.append(str(pos))
                #Annotating each position with all unique Bolli classes
                annot.append(', '.join(set(bolli_var.loc[chrom, pos]['Variant_class'].values)))
                #annot.append(str(bolli_var.loc[chrom, pos]['Variant_class'].values[0]))
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
                        #Annotating each position with all unique Bolli classes
                        ann.append(', '.join(set(bolli_var.loc[chrom, i]['Variant_class'].values)))
                        #ann.append(str(bolli_var.loc[(chrom,i)]['Variant_class'].values[0]))
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
    variants["Bolli_Class"] = cl
    variants["Bolli_Frequency"] = freq
    variants["Bolli_Positions"] = positions
    variants["Bolli_Annotation"] = annot
    return(variants)

def annotate_lohr(variants, lohr_path):
    lohr = pd.read_csv(filepath_or_buffer=lohr_path, sep='\t')
    lohr=lohr[["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification"]]
    lohrC=lohr.groupby(['Chromosome','Start_Position'])['Tumor_Seq_Allele2'].count()
    cl = [] 
    freq = [] 
    positions = [] 
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START")
    for record in variants.itertuples(index=False, name=None):
        flag = 0
        try:    
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            start = int(record[poscol]) - 9
            end = int(record[poscol]) + 9
            if (chrom, pos) in lohrC.index:
                cl.append("genomic_exact")
                freq.append(str(lohrC.loc[(chrom,pos)]))
                positions.append(str(pos))
                flag = 1
            if flag == 0:
                lohrCsub=lohrC.loc[chrom]
                if not lohrCsub[(lohrCsub.index >= start) & (lohrCsub.index <= end)].empty:
                    fr = []
                    posit = []
                    for i in lohrCsub[(lohrCsub.index >= start) & (lohrCsub.index <= end)].index.values:
                        fr.append(str(lohrC.loc[(chrom,i)]))
                        posit.append(str(i))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    positions.append((":".join(posit)))
                else:
                    cl.append(None)
                    freq.append(None)
                    positions.append(None)
        except:
            cl.append(None)
            freq.append(None)
            positions.append(None)
    variants["Lohr_Class"] = cl
    variants["Lohr_Frequency"] = freq
    variants["Lohr_Positions"] = positions
    return(variants)

def annotate_mytype(variants, path_mytype):
    mytype = pd.read_csv(filepath_or_buffer=path_mytype, sep=',')
    mytype = mytype[["CHR", "START", "REF", "ALT", "CONSENSUS_ANNOTATION", "TARGET_VAF"]]
    mytype_counts = mytype.groupby(['CHR','START'])['ALT'].count()
    mytype_var = mytype.drop(['REF','ALT', "TARGET_VAF"], axis = 1) 
    mytype_var = mytype_var.set_index(['CHR', 'START'])
    mytype_alt = mytype.drop(['REF','CONSENSUS_ANNOTATION', "TARGET_VAF"], axis = 1) 
    mytype_alt = mytype_alt.set_index(['CHR', 'START'])
    mytype_med=mytype.groupby(['CHR','START'])['TARGET_VAF'].median()
    mytype_Q25=mytype.groupby(['CHR','START'])['TARGET_VAF'].quantile(q=0.25)
    mytype_Q75=mytype.groupby(['CHR','START'])['TARGET_VAF'].quantile(q=0.75)
    cl = [] 
    freq = [] 
    medVAF = [] 
    Q25 = [] 
    Q75 = [] 
    positions = []
    alt = []
    annot = []
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START")
    for record in variants.itertuples(index=False, name=None):
        flag = 0
        try:
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            start = int(record[poscol]) - 9
            end = int(record[poscol]) + 9
            if (chrom, pos) in  mytype_counts.index:
                cl.append("genomic_exact")
                freq.append(str(mytype_counts.loc[(chrom,pos)]))
                medVAF.append(str(mytype_med.loc[(chrom,pos)]))
                Q25.append(str(mytype_Q25.loc[(chrom,pos)]))
                Q75.append(str(mytype_Q75.loc[(chrom,pos)]))
                positions.append(str(pos))
                #Annotating each position with all unique alternative alleles
                alt.append(', '.join(set(mytype_alt.loc[chrom, pos]['ALT'].values)))
                #Annotating each position with all unique myTYPE consensus annotation
                annot.append(', '.join(set(mytype_var.loc[chrom, pos]['CONSENSUS_ANNOTATION'].values)))
                flag = 1
            if flag == 0: 
                mytype_counts_sub=mytype_counts.loc[chrom]
                if not mytype_counts_sub[(mytype_counts_sub.index >= start) & (mytype_counts_sub.index <= end)].empty:
                    fr = []
                    mv = []
                    Q2 = []
                    Q7 = []
                    posit = []
                    al = []
                    ann = []
                    for i in mytype_counts_sub[(mytype_counts_sub.index >= start) & (mytype_counts_sub.index <= end)].index.values:
                        fr.append(str(mytype_counts.loc[(chrom,i)]))
                        mv.append(str(mytype_med.loc[(chrom,i)]))
                        Q2.append(str(mytype_Q25.loc[(chrom,i)]))
                        Q7.append(str(mytype_Q75.loc[(chrom,i)]))
                        posit.append(str(i))
                        #Annotating each position with all unique alternative alleles
                        al.append(', '.join(set(mytype_alt.loc[chrom, i]['ALT'].values)))
                        #Annotating with all unique myTYPE consensus annotations for position
                        ann.append(', '.join(set(mytype_var.loc[chrom, i]['CONSENSUS_ANNOTATION'].values)))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    medVAF.append((":".join(mv)))
                    Q25.append((":".join(Q2)))
                    Q75.append((":".join(Q7)))
                    positions.append((":".join(posit)))
                    alt.append((":".join(al)))
                    annot.append((":".join(ann)))
                else:
                    cl.append(None)
                    freq.append(None)
                    medVAF.append(None)
                    Q25.append(None)
                    Q75.append(None)
                    positions.append(None)
                    alt.append(None)
                    annot.append(None)
        except:
            cl.append(None)
            freq.append(None)
            medVAF.append(None)
            Q25.append(None)
            Q75.append(None)
            positions.append(None)
            alt.append(None)
            annot.append(None)
    variants["myTYPE_Class"] = cl
    variants["myTYPE_Frequency"] = freq
    variants["myTYPE_VAF"] = medVAF
    variants["myTYPE_Q25"] = Q25
    variants["myTYPE_Q75"] = Q75
    variants["myTYPE_Positions"] = positions
    variants["myTYPE_Alt"] = alt
    variants["myTYPE_Annotation"] = annot
    return(variants)

def annotate_known(variants, mytype):
    #KNOWN_MM = 1 if previously found in MM. Includes any match in MMRF, Bolli and Lohr, and UNKNOWN/LIKELY/ONCOGENIC by mytype

    # Define column with data on whether somatic variant has been identified by myTYPE for that position
    if (mytype): #Only run function if data is passed to the optional variable "mytype"
        mytype_annot = variants['myTYPE_Annotation'].tolist()
        myTYPE_somatic = []
        for entry in mytype_annot:
            if pd.isnull(entry):
                myTYPE_somatic.append(0)
            elif re.search('ONCOGENIC', entry) or re.search('LIKELY', entry) or re.search('UNKNOWN', entry):
                myTYPE_somatic.append(1)
            else:
                myTYPE_somatic.append(0)
        variants['myTYPE_somatic'] = myTYPE_somatic
    else:
        variants['myTYPE_somatic'] = 0

    # Define column KNOWN_MM based on annotation data         
    variants['KNOWN_MM'] = np.where((variants['myTYPE_somatic'] == 1) | (variants['MMRF_Class'].notnull()) | (variants['Bolli_Class'].notnull()) | (variants['Lohr_Class'].notnull()), 1, 0)
    variants=variants.drop('myTYPE_somatic', axis = 1)
    return(variants)

## APPLY FLAGS FOR FILTERING
def filter_gene(variants, genes_bed):
    #args: variants = variants file; genes_bed = path to file with beds to include
    variants_bed = variants[["CHR", "START", "END", "ID_VARIANT"]]
    variants_bed = pyb.BedTool.from_dataframe(variants_bed) #Turning variants file into bed format
    genes = pyb.BedTool(genes_bed) #import list of genes in panel as bed format
    variants_inter = variants_bed.intersect(genes, u=True) #bed file with intersection of panel and input file
    flaglist = [] #empty list for names of variants in intersection bed file

    if not variants_inter.head(n=1, as_string = True) == '': #if bed file is not empty
        flaglist = pyb.BedTool.to_dataframe(variants_inter)['name'] #convert intersect bed file to data frame and subset col with variant ID
    variants['MFLAG_GENE'] = np.where(variants.ID_VARIANT.isin(flaglist), 0, 1) #flag variant if ID is not in overlap list
    return(variants)

def filter_IGH(variants, IGH_path):
    #Remove calls in IGH locus
    variants_bed = variants[["CHR", "START", "END", "ID_VARIANT"]]
    variants_bed = pyb.BedTool.from_dataframe(variants_bed)
    igh = pyb.BedTool(IGH_path)
    variants_inter = variants_bed.intersect(igh, u=True)
    flaglist = []
    if not variants_inter.head(n=1, as_string = True) == '':
        flaglist = pyb.BedTool.to_dataframe(variants_inter)['name']
    variants['MFLAG_IGH'] = np.where(variants.ID_VARIANT.isin(flaglist), 1, 0)
    return(variants)

def filter_MAF(variants):
    variants['MFLAG_MAF'] = np.where(variants['MAX_MAF'] > 0.03, 1, 0)
    return(variants)

def filter_MAF_COSMIC(variants, mode):
    #Remove if >0.1 % MAF and not in COSMIC
    if mode == 'snv':
        #Counts as COSMIC only exact and pos
        variants['MFLAG_MAFCOS'] = np.where((variants['MAX_MAF'] > 0.001) & (variants['ANY_EXACT_POS'] == 0), 1, 0)
    if mode == 'indel':
        #Counts all COSMIC
        variants['MFLAG_MAFCOS'] = np.where((variants['MAX_MAF'] > 0.001) & (variants['COSMIC'].isnull()), 1, 0)
    
    return(variants)

def filter_nonpass(variants, mode):
    #Remove non-PASS calls not present in COSMIC or previous cohorts (MMRF, Bolli, etc.).
    if mode == 'snv':
        #If SNVs, remove only non-synonymous mutations. Counts as COSMIC only exact and pos. 
        drop = ['non_synonymous_codon']
        variants['MFLAG_NONPASS'] = np.where((variants['FILTER'] != "PASS") & (variants['EFFECT'].isin(drop)) & (variants['ANY_EXACT_POS'] == 0) & (variants['KNOWN_MM'] == 0), 1, 0)
        return(variants)
    elif mode == 'indel': 
        #If indels: do not take into account mutation effect, also allow all cosmic mentions as "in cosmic"
        variants['MFLAG_NONPASS'] = np.where((variants['FILTER'] != "PASS") & (variants['COSMIC'].isnull()) & (variants['KNOWN_MM'] == 0), 1, 0)
        return(variants)

def filter_normals(variants):
    #Remove calls present in at least 1 internal normals
    variants['MFLAG_NORM'] = np.where((variants['Normals_Frequency'] > 0), 1, 0)
    return(variants)

def filter_VAF(variants):
    variants['MFLAG_VAF'] = np.where(variants['TARGET_VAF'] < 0.01, 1, 0)
    return(variants)

def filter_BIDIR(variants):
    variants['MFLAG_BIDIR'] = np.where(variants['BIDIR'] == 0, 1, 0)
    return(variants)

## FILTER AND EXPORT
def namecore(infile):
    name = infile.split('/')[-1]
    if re.search('.csv$', name):
        return(re.sub('.csv$', '', name))
    else:
        return(re.sub('.tsv.gz$', '', name))

def filter_export(variants, outdir, name, mode):
    good=variants[variants.filter(regex='MFLAG').sum(axis=1) == 0]
    bad=variants[variants.filter(regex='MFLAG').sum(axis=1) > 0]
    if not re.search('/$', outdir):
        outdir =''.join([outdir,'/'])
    date=str(datetime.today()).split()[0].split("-")
    name='_'.join([name, '_'.join(date)])
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
        #Call the "Version" file for version info?
        f.write(f'Somatic variant processing for myTYPE\nv.1.0\nTime of run start: {run_time.split(".")[0]}\n')
        f.write(f'####\nMode: {mode}\n')
        f.write(f'Imported calls: {variants.shape[0]}\n')
        f.write('Flagging variants for filtering:\n')
        f.write(f'MFLAG_GENE: Gene not in list of genes to keep: {variants["MFLAG_GENE"].sum()}\n')
        f.write(f'MFLAG_IGH: In IGH locus: {variants["MFLAG_IGH"].sum()}\n')
        f.write(f'MFLAG_MAF: MAF > 3 % in exac/1000genomes: {variants["MFLAG_MAF"].sum()}\n')
        f.write(f'MFLAG_MAFCOS: MAF > 0.1 % and not in COSMIC (exact/pos): {variants["MFLAG_MAFCOS"].sum()}\n')
        f.write(f'MFLAG_NONPASS: NON-PASS IF not in cosmic, previously known in MM, not stopgain, splicesite..: {variants["MFLAG_NONPASS"].sum()}\n')
        f.write(f'MFLAG_NORM: Variant in 1 or more good normal: {variants["MFLAG_NORM"].sum()}\n')
        f.write(f'MFLAG_VAF: Remove variants with target VAF < 1 %: {variants["MFLAG_VAF"].sum()}\n')
        f.write(f'MFLAG_BIDIR: Remove variants BIDIR = 0 (only reads on one strand): {variants["MFLAG_BIDIR"].sum(0)}\n')
        f.write(f'Removing calls with >= 1 MFLAG: {bad.shape[0]}\n')
        f.write(f'Calls passed filters: {good.shape[0]}\n')
    return()

# Main Function
def process(
        mode,
        infile,
        skiplines,
        outdir,
        genes,
        genes_bed,
        igh,
        mmrf,
        bolli,
        lohr,
        normals,
        mytype):
    """Main function to process myTYPE SNV and indel output"""
    ##IMPORTING DATA
    variants = import_variants(infile, skiplines)

    ##ANNOTATIONS
    variants = annotate_COSMIC(variants) 
    if (genes): #Only runs function if a path was passed to optional argument "gene"
        variants = annotate_genefreq(variants, genes) #Replace this with mutation frequency from MMRF? (and other raw data?)
    variants = annotate_maf(variants)
    variants = annotate_normals(variants, normals)
    variants = annotate_mmrf(variants, mmrf)
    variants = annotate_bolli(variants, bolli)
    variants = annotate_lohr(variants, lohr)
    if (mytype): #Only runs function if a path was passed to optional argument "mytype"
        variants = annotate_mytype(variants, mytype)
    variants = annotate_known(variants, mytype)

    ##FILTERS
    variants = filter_gene(variants, genes_bed) #Remove calls outside of genes in 'genes_bed'
    variants = filter_IGH(variants, igh) #Remove calls in IGH locus
    variants = filter_MAF(variants) #Remove if MAF > 3 %.
    variants = filter_MAF_COSMIC(variants, mode) #Remove if >0.1 % MAF and not in COSMIC (exact or pos)
    variants = filter_nonpass(variants, mode) #Remove non-PASS calls not present in COSMIC (exact or pos) or previous cohorts (MMRF, Bolli, etc.). EXCEPT nonsense, etc.
    variants = filter_normals(variants) #Remove calls present in at least 1 internal normal
    variants = filter_VAF(variants) #Remove calls with VAF < 1 %
    variants = filter_BIDIR(variants) #Remove calls with BIDIR = 0

    ##OUTPUT
    name = namecore(infile)
    filter_export(variants, outdir, name, mode)
    print('Variant processing complete')

