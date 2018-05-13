"""variants_process main command."""

import re
from datetime import datetime
from os.path import join
import pandas as pd
import numpy as np
import pybedtools as pyb

RUN_TIME = str(datetime.today())

## IMPORT VARIANTS FILE
def import_variants(path, skiplines):
    """
    Determine filetype and import, returns pandas dataFrame
    """
    if re.search('.csv$', path):
        try:
            variants = pd.read_csv(
                filepath_or_buffer=path,
                skiprows=skiplines,
                low_memory=False)
        except NameError:
            raise Exception(f'Error when importing file {path}')

        print(f'Loaded file containing {variants.shape[0]} '
              f'variant calls. Processing...')
        return(variants)
    elif re.search('.tsv.gz$', path):
        try:
            variants = pd.read_csv(
                filepath_or_buffer=path,
                compression='gzip',
                sep='\t',
                skiprows=skiplines,
                low_memory=False)
        except NameError:
            raise Exception(f'Error when importing file {path}')
        print(f'Loaded file containing {variants.shape[0]} '
              f'variant calls. Processing...')
        return(variants)
    else:
        raise Exception(f'Input file {path} has unsupported '
                        f'extension: try .csv or .tsv.gz')

## ANNOTATION FUNCTIONS
def annotate_cosmic(variants):
    """
    Generate columns:
    HEME_EXACT: Number of exact matches for hematopoietic and
    lymphoid tissue in cosmic.
    ANY_EXACT_POS: YES/NO for any EXACT or POS match in cosmic.
    """
    heme_exact = []
    cosmic = variants['COSMIC'].tolist()
    search_1 = 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'
    search_2 = r'(?<=HAEMATOPOIETIC_AND_LYMPHOID_TISSUE=)\w+'
    for entry in cosmic:
        if pd.isnull(entry):
            heme_exact.append(None)
        else:
            first = entry.split('|')[0]
            if re.search('^GENOMIC_EXACT', first):
                if re.search(search_1, first):
                    count = re.search(search_2, first)[0]
                    heme_exact.append(count)
                else:
                    heme_exact.append(None)
            else:
                heme_exact.append(None)
    variants['HEME_EXACT'] = heme_exact
    any_exact_pos = []
    for entry in cosmic:
        if pd.isnull(entry):
            any_exact_pos.append(0)
        elif re.search(
                'GENOMIC_EXACT', entry) or re.search(
                    'GENOMIC_POS', entry):
            any_exact_pos.append(1)
        else:
            any_exact_pos.append(0)
    variants['ANY_EXACT_POS'] = any_exact_pos
    return(variants)

def annotate_genefreq(variants, genes):
    """
    Generate column:
    MAX_MUTFREQ: Maximal mutation frequency in gene
    as previously published in large MM studies.
    """
    freqlist = pd.read_excel(io=genes)
    freqlist['MAX_MUTFREQ'] = round(
        freqlist.filter(regex='freq').max(axis=1), 1)
    freqlist = freqlist[['GENE', 'MAX_MUTFREQ']]
    variants = pd.merge(variants, freqlist, how='left')
    return(variants)

def annotate_maf(variants):
    """
    Generate column:
    MAX_MAF: Maximal MAF of variant in any normal database
    """
    variants['MAX_MAF'] = 0 # Sets variable to 0 if frequency is not reported
    variants['MAX_MAF'] = variants.filter(regex='MAF').max(axis=1)
    return(variants)

def annotate_normals(variants, path_normals):
    """
    Generate columns:
    Normals_Frequency: Number of good normals sequenced by myTYPE who have a
    variant with the same CHR and START.
    Normals_median_VAF: Median VAF of matching variants in normals.
    """
    normals = pd.read_csv(
        filepath_or_buffer=path_normals,
        compression='gzip',
        sep='\t',
        skiprows=1,
        low_memory=False)
    normals_counts = normals.groupby(
        ["CHR", "START", "REF", "ALT"]).size().reset_index(name="count")
    normals_counts = normals_counts[["CHR", "START", "count"]].set_index(
        ['CHR', 'START'])
    normals_median = normals.groupby(
        ['CHR', 'START'])['TARGET_VAF'].median()
    normalC = []
    normalVAF = []
    chrcol = variants.columns.get_loc("CHR")
    poscol = variants.columns.get_loc("START")
    for record in variants.itertuples(index=False, name=None):
        try:
            chrom = str(record[chrcol])
            pos = int(record[poscol])
            tempC = normals_counts.loc[(chrom, pos), "count"]
            tempC = tempC.ix[0]
            normalC.append(int(tempC))
            normalVAF.append(str(normals_median.loc[(chrom, pos)]))
        except KeyError:
            normalC.append(0)
            normalVAF.append("0")
    variants["Normals_Frequency"] = normalC
    variants["Normals_median_VAF"] = normalVAF
    return(variants)

def annotate_mmrf(variants, path_mmrf):
    """
    Generate columns:
    MMRF_Class: Exact or close (within 10 bp) match for CHR and START
    MMRF_Frequency: Number of matches
    MMRF_VAF: Median VAF of matches
    MMRF_Q25: 25th VAFF-quartile of matches
    MMRF_Q75: 75th VAF-quartile of matches
    MMRF_Positions: START position of matches
    """
    mmrf = pd.read_csv(filepath_or_buffer=path_mmrf, sep='\t')
    mmrf = mmrf[["Sample", "#CHROM", "POS", "REF", "ALT",
                 "GEN[0].AR", "GEN[1].AR"]]
    mmrf = mmrf.drop_duplicates() ## What are these duplicates?
    mmrfM = mmrf.groupby(['#CHROM', 'POS'])['GEN[1].AR'].median()
    mmrfC = mmrf.groupby(['#CHROM', 'POS'])['GEN[1].AR'].count()
    mmrfQ25 = mmrf.groupby(['#CHROM', 'POS'])['GEN[1].AR'].quantile(q=0.25)
    mmrfQ75 = mmrf.groupby(['#CHROM', 'POS'])['GEN[1].AR'].quantile(q=0.75)
    cl = []
    freq = []
    medVAF = []
    q25 = []
    q75 = []
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
                freq.append(str(mmrfC.loc[(chrom, pos)]))
                medVAF.append(str(mmrfM.loc[(chrom, pos)]))
                q25.append(str(mmrfQ25.loc[(chrom, pos)]))
                q75.append(str(mmrfQ75.loc[(chrom, pos)]))
                positions.append(str(pos))
                flag = 1
            if flag == 0:
                mmrfCsub = mmrfC.loc[chrom]
                if not mmrfCsub[(mmrfCsub.index >= start) &
                                (mmrfCsub.index <= end)].empty:
                    fr = []
                    mv = []
                    q2 = []
                    q7 = []
                    posit = []
                    for i in mmrfCsub[(mmrfCsub.index >= start) &
                                      (mmrfCsub.index <= end)].index.values:
                        fr.append(str(mmrfC.loc[(chrom, i)]))
                        mv.append(str(mmrfM.loc[(chrom, i)]))
                        q2.append(str(mmrfQ25.loc[(chrom, i)]))
                        q7.append(str(mmrfQ75.loc[(chrom, i)]))
                        posit.append(str(i))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    medVAF.append((":".join(mv)))
                    q25.append((":".join(q2)))
                    q75.append((":".join(q7)))
                    positions.append((":".join(posit)))
                else:
                    cl.append(None)
                    freq.append(None)
                    medVAF.append(None)
                    q25.append(None)
                    q75.append(None)
                    positions.append(None)
        except KeyError:
            cl.append(None)
            freq.append(None)
            medVAF.append(None)
            q25.append(None)
            q75.append(None)
            positions.append(None)
    variants["MMRF_Class"] = cl
    variants["MMRF_Frequency"] = freq
    variants["MMRF_VAF"] = medVAF
    variants["MMRF_Q25"] = q25
    variants["MMRF_Q75"] = q75
    variants["MMRF_Positions"] = positions
    return(variants)

def annotate_bolli(variants, path_bolli):
    """
    Generate columns:
    Bolli_Class: Exact or close (within 10 bp) match for CHR and START
    Bolli_Frequency: Number of matches
    Bolli_Positions: START position of matches
    Bolli_Annotation: Manual annotation category, may include more than one if
    variant has been non-uniformly annotated
    """
    bolli = pd.read_csv(filepath_or_buffer=path_bolli, sep='\t')
    bolli = bolli[["CHR", "START", "WT", "MT", "Variant_class"]]
    bolli_counts = bolli.groupby(['CHR', 'START'])['MT'].count()
    bolli_var = bolli.drop(['WT', 'MT'], axis=1)
    bolli_var = bolli_var.set_index(['CHR', 'START'])
    bolli_var = bolli_var.sort_index()
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
                freq.append(str(bolli_counts.loc[(chrom, pos)]))
                positions.append(str(pos))
                # Annotating each position with all unique Bolli classes
                annot.append(', '.join(
                    set(bolli_var.loc[chrom, pos].values.flat)))
                flag = 1
            if flag == 0:
                bolli_sub = bolli_counts.loc[chrom]
                bolli_sub = bolli_sub[(bolli_sub.index >= start) &
                                      (bolli_sub.index <= end)]
                if not bolli_sub.empty:
                    fr = []
                    posit = []
                    ann = []
                    for i in bolli_sub.index.values:
                        fr.append(str(bolli_counts.loc[(chrom, i)]))
                        posit.append(str(i))
                        ann.append(', '.join(set(
                            bolli_var.loc[chrom, i].values.flat)))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    positions.append((":".join(posit)))
                    annot.append((":".join(ann)))
                else:
                    cl.append(None)
                    freq.append(None)
                    positions.append(None)
                    annot.append(None)
        except KeyError:
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
    """
    Generate columns:
    Lohr_Class: Exact or close (within 10 bp) match for CHR and START
    Lohr_Frequency: Number of matches
    Lohr_Positions: START position of matches
    """
    lohr = pd.read_csv(filepath_or_buffer=lohr_path, sep='\t')
    lohr = lohr[["Tumor_Sample_Barcode", "Chromosome",
                 "Start_Position", "Reference_Allele",
                 "Tumor_Seq_Allele2", "Variant_Classification"]]
    lohrC = lohr.groupby(
        ['Chromosome', 'Start_Position'])['Tumor_Seq_Allele2'].count()
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
                freq.append(str(lohrC.loc[(chrom, pos)]))
                positions.append(str(pos))
                flag = 1
            if flag == 0:
                lohrCsub = lohrC.loc[chrom]
                lohrCsub = lohrCsub[(lohrCsub.index >= start) &
                                    (lohrCsub.index <= end)]
                if not lohrCsub.empty:
                    fr = []
                    posit = []
                    for i in lohrCsub.index.values:
                        fr.append(str(lohrC.loc[(chrom, i)]))
                        posit.append(str(i))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    positions.append((":".join(posit)))
                else:
                    cl.append(None)
                    freq.append(None)
                    positions.append(None)
        except KeyError:
            cl.append(None)
            freq.append(None)
            positions.append(None)
    variants["Lohr_Class"] = cl
    variants["Lohr_Frequency"] = freq
    variants["Lohr_Positions"] = positions
    return(variants)

def annotate_mytype(variants, path_mytype):
    """
    Generate columns:
    myTYPE_Class: Exact or close (within 10 bp) match for CHR and START
    myTYPE_Frequency: Number of matches
    myTYPE_VAF: Median VAF of matches
    myTYPE_Q25: 25th VAF-quartile of matches
    myTYPE_Q75: 75th VAF-quartile of matches
    myTYPE_Positions: START position of matches
    myTYPE_Alt: All unique alternative alleles with START == position
    myTYPE_Annotation: Manual annotation category, may include more than one if
    variant has been non-uniformly annotated
    """
    mytype = pd.read_csv(filepath_or_buffer=path_mytype, sep=',')
    mytype = mytype[["CHR", "START", "REF", "ALT",
                     "MANUAL_ANNOTATION", "TARGET_VAF"]]
    mytype_counts = mytype.groupby(['CHR', 'START'])['ALT'].count()
    mytype_var = mytype.drop(['REF', 'ALT', "TARGET_VAF"], axis=1)
    mytype_var = mytype_var.set_index(['CHR', 'START'])
    mytype_alt = mytype.drop(
        ['REF', 'MANUAL_ANNOTATION', "TARGET_VAF"], axis=1)
    mytype_alt = mytype_alt.set_index(['CHR', 'START'])
    mytype_med = mytype.groupby(['CHR', 'START'])['TARGET_VAF'].median()
    mytype_Q25 = mytype.groupby(['CHR', 'START'])['TARGET_VAF'].quantile(q=0.25)
    mytype_Q75 = mytype.groupby(['CHR', 'START'])['TARGET_VAF'].quantile(q=0.75)
    mytype_alt = mytype_alt.sort_index()
    mytype_var = mytype_var.sort_index()
    cl = []
    freq = []
    medVAF = []
    q25 = []
    q75 = []
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
                freq.append(str(mytype_counts.loc[(chrom, pos)]))
                medVAF.append(str(mytype_med.loc[(chrom, pos)]))
                q25.append(str(mytype_Q25.loc[(chrom, pos)]))
                q75.append(str(mytype_Q75.loc[(chrom, pos)]))
                positions.append(str(pos))
                # Annotating with all unique alternative alleles
                alt.append(', '.join(
                    set(mytype_alt.loc[chrom, pos].values.flat)))
                # Annotating with all unique myTYPE consensus annotations
                annot.append(', '.join(
                    set(mytype_var.loc[chrom, pos].values.flat)))
                flag = 1
            if flag == 0:
                mytype_sub = mytype_counts.loc[chrom]
                mytype_sub = mytype_sub[(mytype_sub.index >= start) &
                                        (mytype_sub.index <= end)]
                if not mytype_sub.empty:
                    fr = []
                    mv = []
                    q2 = []
                    q7 = []
                    posit = []
                    al = []
                    ann = []
                    for i in mytype_sub.index.values:
                        fr.append(str(mytype_counts.loc[(chrom, i)]))
                        mv.append(str(mytype_med.loc[(chrom, i)]))
                        q2.append(str(mytype_Q25.loc[(chrom, i)]))
                        q7.append(str(mytype_Q75.loc[(chrom, i)]))
                        posit.append(str(i))
                        al.append(', '.join(
                            set(mytype_alt.loc[chrom, i].values.flat)))
                        ann.append(', '.join(
                            set(mytype_var.loc[chrom, i].values.flat)))
                    cl.append("genomic_close")
                    freq.append((":".join(fr)))
                    medVAF.append((":".join(mv)))
                    q25.append((":".join(q2)))
                    q75.append((":".join(q7)))
                    positions.append((":".join(posit)))
                    alt.append((":".join(al)))
                    annot.append((":".join(ann)))
                else:
                    cl.append(None)
                    freq.append(None)
                    medVAF.append(None)
                    q25.append(None)
                    q75.append(None)
                    positions.append(None)
                    alt.append(None)
                    annot.append(None)
        except KeyError:
            cl.append(None)
            freq.append(None)
            medVAF.append(None)
            q25.append(None)
            q75.append(None)
            positions.append(None)
            alt.append(None)
            annot.append(None)
    variants["myTYPE_Class"] = cl
    variants["myTYPE_Frequency"] = freq
    variants["myTYPE_VAF"] = medVAF
    variants["myTYPE_Q25"] = q25
    variants["myTYPE_Q75"] = q75
    variants["myTYPE_Positions"] = positions
    variants["myTYPE_Alt"] = alt
    variants["myTYPE_Annotation"] = annot
    return(variants)

def annotate_known(variants, mytype):
    """
    Generate columns:
    KNOWN_MM = 1 if previously found in MM. Includes any match in MMRF,
    Bolli and Lohr, and UNKNOWN/LIKELY/ONCOGENIC by mytype
    """

    # Only run function if data is passed to the optional variable "mytype"
    if mytype:
        mytype_annot = variants['myTYPE_Annotation'].tolist()
        myTYPE_somatic = []
        for entry in mytype_annot:
            if pd.isnull(entry):
                myTYPE_somatic.append(0)
            else:
                search_1 = re.search('ONCOGENIC', entry)
                search_2 = re.search('LIKELY', entry)
                search_3 = re.search('UNKNOWN', entry)
                if search_1 or search_2 or search_3:
                    myTYPE_somatic.append(1)
                else:
                    myTYPE_somatic.append(0)
        variants['myTYPE_somatic'] = myTYPE_somatic
    else:
        variants['myTYPE_somatic'] = 0

    # Define column KNOWN_MM based on annotation data
    variants['KNOWN_MM'] = np.where((variants['myTYPE_somatic'] == 1) |
                                    (variants['MMRF_Class'].notnull()) |
                                    (variants['Bolli_Class'].notnull()) |
                                    (variants['Lohr_Class'].notnull()), 1, 0)

    variants = variants.drop('myTYPE_somatic', axis=1)
    return(variants)

## APPLY FLAGS FOR FILTERING
def filter_panel(variants, genes_bed):
    """
    Filter MFLAG_PANEL: 1 if variant is not in BED file of regions to keep
    """
    variants_bed = variants[["CHR", "START", "END", "ID_VARIANT"]]
    # Turning variants file into bed format
    variants_bed = pyb.BedTool.from_dataframe(variants_bed)
    # Import list of genes in panel as bed format
    genes = pyb.BedTool(genes_bed)
    # Bed file with intersection of panel and input file
    variants_inter = variants_bed.intersect(genes, u=True)
    # Empty list for names of variants in intersection bed file
    flaglist = []

    # If bed file is not empty
    if not variants_inter.head(n=1, as_string=True) == '':
        # Convert intersect bed file to data frame; subset col with variant ID
        flaglist = pyb.BedTool.to_dataframe(variants_inter)['name']
    # Flag variant if ID is not in overlap list
    variants['MFLAG_PANEL'] = np.where(variants.ID_VARIANT.isin(flaglist), 0, 1)
    return(variants)

def filter_drop(variants, genes_drop):
    """
    Filter MFLAG_DROP: 1 if variant is in list of genes to drop.
    """
    drop = pd.read_excel(io=genes_drop)['GENE']
    variants['MFLAG_DROP'] = np.where(variants.GENE.isin(drop), 1, 0)
    return(variants)

def filter_igh(variants, igh_path):
    """
    Filter MFLAG_IGH: 1 if variant in IGH locus
    """
    variants_bed = variants[["CHR", "START", "END", "ID_VARIANT"]]
    variants_bed = pyb.BedTool.from_dataframe(variants_bed)
    igh = pyb.BedTool(igh_path)
    variants_inter = variants_bed.intersect(igh, u=True)
    flaglist = []
    if not variants_inter.head(n=1, as_string=True) == '':
        flaglist = pyb.BedTool.to_dataframe(variants_inter)['name']
    variants['MFLAG_IGH'] = np.where(variants.ID_VARIANT.isin(flaglist), 1, 0)
    return(variants)

def filter_maf(variants):
    """
    Filter MFLAG_MAF: 1 if variant MAF > 3 % in exac/1000genomes
    """
    variants['MFLAG_MAF'] = np.where(variants['MAX_MAF'] > 0.03, 1, 0)
    return(variants)

def filter_maf_cosmic(variants, mode):
    """
    Filter MFLAG_MAFCOS: 1 if variant has >0.1 % MAF and not in COSMIC
    For SNVs: Only counts exact and pos as in cosmic
    For Indels: Counts all COSMIC.
    """
    if mode == 'snv':
        variants['MFLAG_MAFCOS'] = np.where(
            (variants['MAX_MAF'] > 0.001) &
            (variants['ANY_EXACT_POS'] == 0), 1, 0)
    if mode == 'indel':
        variants['MFLAG_MAFCOS'] = np.where(
            (variants['MAX_MAF'] > 0.001) &
            (variants['COSMIC'].isnull()), 1, 0)
    return(variants)

def filter_nonpass(variants, mode):
    """
    Filter MFLAG_MAF: 1 if NON-PASS AND not in cosmic or previously known in MM
    Counts SNVs and Indels as "in cosmic" like for MAFCOS flag.
    For SNV: Only removes missense mutations with this flag
    """
    if mode == 'snv':
        drop = ['non_synonymous_codon']
        variants['MFLAG_NONPASS'] = np.where(
            (variants['FILTER'] != "PASS") &
            (variants['EFFECT'].isin(drop)) &
            (variants['ANY_EXACT_POS'] == 0) &
            (variants['KNOWN_MM'] == 0), 1, 0)
        return(variants)
    variants['MFLAG_NONPASS'] = np.where(
        (variants['FILTER'] != "PASS") &
        (variants['COSMIC'].isnull()) &
        (variants['KNOWN_MM'] == 0), 1, 0)
    return(variants)

def filter_normals(variants):
    """
    Filter MFLAG_NORM: 1 if present in 1 or more good normal
    """
    variants['MFLAG_NORM'] = np.where((variants['Normals_Frequency'] > 0), 1, 0)
    return(variants)

def filter_vaf(variants):
    """
    Filter MFLAG_VAF: 1 if target VAF < 1 %
    """
    variants['MFLAG_VAF'] = np.where(variants['TARGET_VAF'] < 0.01, 1, 0)
    return(variants)

def filter_bidir(variants):
    """
    Filter MFLAG_BIDIR: 1 if BIDIR = 0
    """
    variants['MFLAG_BIDIR'] = np.where(variants['BIDIR'] == 0, 1, 0)
    return(variants)

## FILTER AND EXPORT
def namecore(infile):
    """
    Returns the "core" of the input file name, for use in output files.
    """
    name = infile.split('/')[-1]
    if re.search('.csv$', name):
        return(re.sub('.csv$', '', name))
    return(re.sub('.tsv.gz$', '', name))

def filter_export(variants, outdir, name, mode):
    """
    Function properties:
    1. Filters variants into "good" or "bad" based on flags.
    2. Writes files with good and bad variants.
    3. Creates processing summary report.
    """
    # Filtering
    good = variants[variants.filter(regex='MFLAG').sum(axis=1) == 0]
    bad = variants[variants.filter(regex='MFLAG').sum(axis=1) > 0]

    # Define output names
    date = str(datetime.today()).split()[0].split("-")
    name = '_'.join([name, '_'.join(date)])
    goodname = join(outdir, name + '_goodcalls.csv')
    badname = join(outdir, name + '_badcalls.csv')
    textname = join(outdir, name + '_report.txt')

    # Export files
    good.to_csv(
        path_or_buf=goodname,
        index=False)
    bad.to_csv(
        path_or_buf=badname,
        index=False)

    # Summary report
    with open(textname, 'w') as f:
        # Call the "Version" file for version info?
        f.write(
            f'Somatic variant processing for myTYPE\nv.1.0\nTime '
            f'of run start: {RUN_TIME.split(".")[0]}\n'
        )
        f.write(f'####\nMode: {mode}\n')
        f.write(f'Imported calls: {variants.shape[0]}\n')
        f.write('Flagging variants for filtering:\n')
        f.write(f'MFLAG_PANEL: Variant not in BED file of '
                f'regions to keep: {variants["MFLAG_PANEL"].sum()}\n')
        f.write(f'MFLAG_DROP: Variant in excluded gene: '
                f'{variants["MFLAG_DROP"].sum()}\n')
        f.write(f'MFLAG_IGH: In IGH locus: {variants["MFLAG_IGH"].sum()}\n')
        f.write(f'MFLAG_MAF: MAF > 3 % in exac/1000genomes: '
                f'{variants["MFLAG_MAF"].sum()}\n')
        f.write(f'MFLAG_MAFCOS: MAF > 0.1 % and not in COSMIC '
                f'(exact/pos): {variants["MFLAG_MAFCOS"].sum()}\n')
        f.write(f'MFLAG_NONPASS: NON-PASS IF not in cosmic, previously '
                f'known in MM, not stopgain, splicesite..: '
                f'{variants["MFLAG_NONPASS"].sum()}\n')
        f.write(f'MFLAG_NORM: Variant in 1 or more good normal: '
                f'{variants["MFLAG_NORM"].sum()}\n')
        f.write(f'MFLAG_VAF: Remove variants with target '
                f'VAF < 1 %: {variants["MFLAG_VAF"].sum()}\n')
        f.write(f'MFLAG_BIDIR: Remove variants BIDIR = 0 (only reads '
                f'on one strand): {variants["MFLAG_BIDIR"].sum(0)}\n')
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
        genes_drop,
        genes_bed,
        igh,
        mmrf,
        bolli,
        lohr,
        normals,
        mytype):
    """Main function to process myTYPE SNV and indel output"""
    ## IMPORTING DATA
    variants = import_variants(infile, skiplines)

    ## ANNOTATIONS
    variants = annotate_cosmic(variants)
    if genes:
        # Only runs if a path was passed to optional argument "gene"
        variants = annotate_genefreq(variants, genes)
        # Replace this with mutation frequency from MMRF? (and other raw data?)
    variants = annotate_maf(variants)
    variants = annotate_normals(variants, normals)
    variants = annotate_mmrf(variants, mmrf)
    variants = annotate_bolli(variants, bolli)
    variants = annotate_lohr(variants, lohr)
    if mytype:
        # Only runs if a path was passed to optional argument "mytype"
        variants = annotate_mytype(variants, mytype)
    variants = annotate_known(variants, mytype)

    ## FILTERS
    variants = filter_panel(variants, genes_bed)
    variants = filter_drop(variants, genes_drop)
    variants = filter_igh(variants, igh)
    variants = filter_maf(variants)
    variants = filter_maf_cosmic(variants, mode)
    variants = filter_nonpass(variants, mode)
    variants = filter_normals(variants)
    variants = filter_vaf(variants)
    variants = filter_bidir(variants)

    ## OUTPUT
    name = namecore(infile)
    filter_export(variants, outdir, name, mode)
    print('Variant processing complete')
    return(variants) # Added this here - may be necessary for test?
