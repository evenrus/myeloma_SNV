"""variants_process main command."""

import timeit
import re
from datetime import datetime
from os.path import join
import pandas as pd
import numpy as np
import pybedtools as pyb

START = timeit.default_timer()

## IMPORT VARIANTS FILE
def import_variants(path):
    """
    Determine filetype and import, returns pandas dataFrame
    """
    if re.search('.csv$', path):
        try:
            variants = pd.read_csv(
                filepath_or_buffer=path,
                comment='#',
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
                comment='#',
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
    # Rename GENE variable to VAG_GENE
    variants = pd.merge(variants, freqlist, how='left')
    return(variants)

def annotate_maf(variants):
    """
    Generate column:
    MAX_MAF: Maximal MAF of variant in any normal database
    """
    variants['MAX_MAF'] = 0 # Sets variable to 0 if frequency is not reported
    variants['MAX_MAF'] = variants.filter(regex='_AF$|_AF_').max(axis=1)
    return(variants)

def annotate_vaf_dep(variants):
    """
    Generate columns:
    TARGET_VAF: Mean TARGET_VAF from different callers
    REFERENCE_VAF: Mean REFERENCE_VAF from different callers
    TARGET_DEPTH: Mean TARGET_DEPTH from different callers
    REFERENCE_DEPTH: Mean REFERENCE_DEPTH from different callers
    """
    variants['TARGET_VAF'] = 0
    variants['TARGET_VAF'] = variants.filter(regex='TARGET_VAF').mean(axis=1)
    variants['REFERENCE_VAF'] = 0
    variants['REFERENCE_VAF'] = variants.filter(regex='REFERENCE_VAF').mean(axis=1)
    variants['TARGET_DEPTH'] = 0
    variants['TARGET_DEPTH'] = variants.filter(regex='TARGET_DEPTH').mean(axis=1)
    variants['REFERENCE_DEPTH'] = 0
    variants['REFERENCE_DEPTH'] = variants.filter(regex='REFERENCE_DEPTH').mean(axis=1)
    variants['DIRPROP'] = variants.filter(regex='DIRPROP').mean(axis=1)
    return(variants)

def annotate_normals(variants, path_normals):
    """
    Annotates variants with internal normal controls:
    Class:      Close   (chr, start within 10 bp),
                Pos     (chr, start),
                Exact   (chr, start, ref, alt)
    Positions:  START position
    Change:     REF > ALT
    """
    normals = pd.read_csv(
        filepath_or_buffer=path_normals)

    normals = normals.set_index(['CHROM','POS']).sort_index()
    
    def annot_row(row, data):
        thres = 10
        chrom = str(row['CHR'])
        start = row['START']
        po = (chrom, start) in data.index
        close = data.ix[(chrom, start-thres):(chrom, start+thres)]
        if po:
            pos = data.loc[(chrom, start)]
            exact = pos[(pos['REF'] == row['REF']) & (pos['ALT'] == row['ALT'])]
            if len(exact) > 0:
                ex_out = ['genomic_exact',
                          start,
                          exact['REF'].iloc[0] + '>' + exact['ALT'].iloc[0]
                         ]
                return pd.Series(ex_out)
            else:
                pos_out = ['genomic_pos',
                           ', '.join([str(i) for i in pos.index.\
                                      get_level_values('POS').tolist()]),
                           ', '.join([str(a) + '>' + str(b) for a, b in \
                           zip(pos['REF'], pos['ALT'])])
                          ]
                return pd.Series(pos_out)
        elif close.shape[0] > 0:
            cl_out = ['genomic_close',
                      ', '.join([str(i) for i in close.index.\
                                 get_level_values('POS').tolist()]),
                      ', '.join(set([str(a) + '>' + str(b) for a, b in \
                      zip(close['REF'].tolist(), close['ALT'].tolist())]))
                     ]
            return pd.Series(cl_out)
        else:
            return pd.Series([None]*3)

    out_names = ["_Class", "_Position", "_Change"]
    out_names = ['Normals' + s for s in out_names]

    variants[out_names] = variants.apply(lambda row: annot_row(row, normals),
                                         axis=1)
    return(variants)

def annotate_mmrf(variants, path_mmrf):
    """
    Annotates variants with MMRF data:
    Class:      Close   (chr, start within 10 bp),
                Pos     (chr, start),
                Exact   (chr, start, ref, alt)
    Frequency:  Number of matches
    VAF:        Median VAF
    Q25:        25th VAF-quartile
    Q75:        75th VAF-quartile
    Positions:  START position
    Change:     REF > ALT
    """
    mmrf = pd.read_csv(filepath_or_buffer=path_mmrf, sep='\t')
    mmrf = mmrf[["#CHROM", "POS", "REF", "ALT", "GEN[1].AR"]]
    mmrf = mmrf.drop_duplicates() ## What are these duplicates?
    mmrf.columns = ["CHR", "START", "REF", "ALT", "TARGET_VAF"]

    def annot_row(row, data):
        thres = 10
        subdat = data[data['CHR'].astype(str) == str(row['CHR'])]
        po = row['START'] in subdat['START'].as_matrix().astype(int)
        close = (abs(subdat['START'].as_matrix() \
                 .astype(int) - row['START']) < thres)
        if po:
            pos = subdat[subdat['START'] == row['START']]
            exact = pos[(pos['REF'] == row['REF']) & (pos['ALT'] == row['ALT'])]
            if len(exact) > 0:
                ex_out = ['genomic_exact',
                          exact['REF'].count(),
                          exact['TARGET_VAF'].median(),
                          exact['TARGET_VAF'].quantile(q=0.25),
                          exact['TARGET_VAF'].quantile(q=0.75),
                          ', '.join(set(exact['START'].astype(str))),
                          ', '.join(set([str(a) + '>' + str(b) for a, b in \
                          zip(exact['REF'].tolist(), exact['ALT'].tolist())]))
                         ]
                return pd.Series(ex_out)
            else:
                pos_out = ['genomic_pos',
                           pos['REF'].count(),
                           pos['TARGET_VAF'].median(),
                           pos['TARGET_VAF'].quantile(q=0.25),
                           pos['TARGET_VAF'].quantile(q=0.75),
                           ', '.join(set(pos['START'].astype(str))),
                           ', '.join(set([str(a) + '>' + str(b) for a, b in \
                           zip(pos['REF'].tolist(), pos['ALT'].tolist())]))
                          ]
                return pd.Series(pos_out)
        elif close.any():
            close = subdat[close]
            cl_out = ['genomic_close',
                      ', '.join(close.groupby(['ALT', 'REF']).size() \
                      .astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .median().astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .quantile(q=0.25).astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .quantile(q=0.75).astype(str).tolist()),
                      ', '.join(set(close['START'].astype(str))),
                      ', '.join(set([str(a) + '>' + str(b) for a, b in \
                      zip(close['REF'].tolist(), close['ALT'].tolist())]))
                     ]
            return pd.Series(cl_out)
        else:
            return pd.Series([None]*7)

    out_names = ["_Class", "_Frequency", "_VAF", "_Q25", "_Q75",
                 "_Position", "_Change"]
    out_names = ['MMRF' + s for s in out_names]

    variants[out_names] = variants.apply(lambda row: annot_row(row, mmrf),
                                         axis=1)
    return(variants)

def annotate_bolli(variants, path_bolli):
    """
    Annotates variants with Bolli data:
    Class:      Close   (chr, start within 10 bp),
                Pos     (chr, start),
                Exact   (chr, start, ref, alt)
    Frequency:  Number of matches
    Positions:  START position
    Change:     REF > ALT
    Annotation: Manual annotation category.
    """
    bolli = pd.read_csv(filepath_or_buffer=path_bolli, sep='\t')
    bolli = bolli[["CHR", "START", "WT", "MT", "Variant_class"]]
    bolli.columns = ["CHR", "START", "REF", "ALT", "ANNOTATION"]
    def annot_row(row, data):
        thres = 10
        subdat = data[data['CHR'].astype(str) == str(row['CHR'])]
        po = row['START'] in subdat['START'].as_matrix().astype(int)
        close = (abs(subdat['START'].as_matrix() \
                 .astype(int) - row['START']) < thres)
        if po:
            pos = subdat[subdat['START'] == row['START']]
            exact = pos[(pos['REF'] == row['REF']) & (pos['ALT'] == row['ALT'])]
            if len(exact) > 0:
                ex_out = ['genomic_exact',
                          exact['REF'].count(),
                          ', '.join(set(exact['START'].astype(str))),
                          ', '.join(set([str(a) + '>' + str(b) for a, b in \
                          zip(exact['REF'].tolist(), exact['ALT'].tolist())])),
                          ', '.join(set(exact['ANNOTATION']))
                         ]
                return pd.Series(ex_out)
            else:
                pos_out = ['genomic_pos',
                           pos['REF'].count(),
                           ', '.join(set(pos['START'].astype(str))),
                           ', '.join(set([str(a) + '>' + str(b) for a, b in \
                           zip(pos['REF'].tolist(), pos['ALT'].tolist())])),
                           ', '.join(set(pos['ANNOTATION']))
                          ]
                return pd.Series(pos_out)
        elif close.any():
            close = subdat[close]
            cl_out = ['genomic_close',
                      ', '.join(close.groupby(['ALT', 'REF']).size() \
                      .astype(str).tolist()),
                      ', '.join(set(close['START'].astype(str))),
                      ', '.join(set([str(a) + '>' + str(b) for a, b in \
                      zip(close['REF'].tolist(), close['ALT'].tolist())])),
                      ', '.join(set(close['ANNOTATION']))
                     ]
            return pd.Series(cl_out)
        else:
            return pd.Series([None]*5)

    out_names = ["_Class", "_Frequency",
                 "_Position", "_Change", "_Annotation"]
    out_names = ['Bolli' + s for s in out_names]

    variants[out_names] = variants.apply(lambda row: annot_row(row, bolli),
                                         axis=1)
    return(variants)

def annotate_lohr(variants, lohr_path):
    """
    Annotates variants with lohr data:
    Class:      Close   (chr, start within 10 bp),
                Pos     (chr, start),
                Exact   (chr, start, ref, alt)
    Frequency:  Number of matches
    Positions:  START position
    Change:     REF > ALT
    """
    lohr = pd.read_csv(filepath_or_buffer=lohr_path, sep='\t')
    lohr = lohr[["Chromosome", "Start_Position", "Reference_Allele",
                 "Tumor_Seq_Allele2"]]
    lohr.columns = ["CHR", "START", "REF", "ALT"]

    def annot_row(row, data):
        thres = 10
        subdat = data[data['CHR'].astype(str) == str(row['CHR'])]
        po = row['START'] in subdat['START'].as_matrix().astype(int)
        close = (abs(subdat['START'].as_matrix() \
                 .astype(int) - row['START']) < thres)
        if po:
            pos = subdat[subdat['START'] == row['START']]
            exact = pos[(pos['REF'] == row['REF']) & (pos['ALT'] == row['ALT'])]
            if len(exact) > 0:
                ex_out = ['genomic_exact',
                          exact['REF'].count(),
                          ', '.join(set(exact['START'].astype(str))),
                          ', '.join(set([str(a) + '>' + str(b) for a, b in \
                          zip(exact['REF'].tolist(), exact['ALT'].tolist())]))
                         ]
                return pd.Series(ex_out)
            else:
                pos_out = ['genomic_pos',
                           pos['REF'].count(),
                           ', '.join(set(pos['START'].astype(str))),
                           ', '.join(set([str(a) + '>' + str(b) for a, b in \
                           zip(pos['REF'].tolist(), pos['ALT'].tolist())]))
                          ]
                return pd.Series(pos_out)
        elif close.any():
            close = subdat[close]
            cl_out = ['genomic_close',
                      ', '.join(close.groupby(['ALT', 'REF']).size() \
                      .astype(str).tolist()),
                      ', '.join(set(close['START'].astype(str))),
                      ', '.join(set([str(a) + '>' + str(b) for a, b in \
                      zip(close['REF'].tolist(), close['ALT'].tolist())]))
                     ]
            return pd.Series(cl_out)
        else:
            return pd.Series([None]*4)

    out_names = ["_Class", "_Frequency", 
                 "_Position", "_Change"]
    out_names = ['Lohr' + s for s in out_names]

    variants[out_names] = variants.apply(lambda row: annot_row(row, lohr),
                                         axis=1)
    return(variants)

def annotate_mytype(variants, path_mytype):
    """
    Annotates variants with previous myTYPE data:
    Class:      Close   (chr, start within 10 bp),
                Pos     (chr, start),
                Exact   (chr, start, ref, alt)
    Frequency:  Number of matches
    VAF:        Median VAF
    Q25:        25th VAF-quartile
    Q75:        75th VAF-quartile
    Positions:  START position
    Change:     REF > ALT
    Annotation: Manual annotation category.
    """
    mytype = pd.read_csv(filepath_or_buffer=path_mytype, sep=',')
    mytype = mytype[["CHR", "START", "REF", "ALT",
                     "MANUAL_ANNOTATION", "TARGET_VAF"]]
    mytype.columns = ["CHR", "START", "REF", "ALT",
                      "ANNOTATION", "TARGET_VAF"]

    def annot_row(row, data):
        thres = 10
        subdat = data[data['CHR'].astype(str) == str(row['CHR'])]
        po = row['START'] in subdat['START'].as_matrix().astype(int)
        close = (abs(subdat['START'].as_matrix() \
                 .astype(int) - row['START']) < thres)
        if po:
            pos = subdat[subdat['START'] == row['START']]
            exact = pos[(pos['REF'] == row['REF']) & (pos['ALT'] == row['ALT'])]
            if len(exact) > 0:
                ex_out = ['genomic_exact',
                          exact['REF'].count(),
                          exact['TARGET_VAF'].median(),
                          exact['TARGET_VAF'].quantile(q=0.25),
                          exact['TARGET_VAF'].quantile(q=0.75),
                          ', '.join(set(exact['START'].astype(str))),
                          ', '.join(set([str(a) + '>' + str(b) for a, b in \
                          zip(exact['REF'].tolist(), exact['ALT'].tolist())])),
                          ', '.join(set(exact['ANNOTATION']))
                         ]
                return pd.Series(ex_out)
            else:
                pos_out = ['genomic_pos',
                           pos['REF'].count(),
                           pos['TARGET_VAF'].median(),
                           pos['TARGET_VAF'].quantile(q=0.25),
                           pos['TARGET_VAF'].quantile(q=0.75),
                           ', '.join(set(pos['START'].astype(str))),
                           ', '.join(set([str(a) + '>' + str(b) for a, b in \
                           zip(pos['REF'].tolist(), pos['ALT'].tolist())])),
                           ', '.join(set(pos['ANNOTATION']))
                          ]
                return pd.Series(pos_out)
        elif close.any():
            close = subdat[close]
            cl_out = ['genomic_close',
                      ', '.join(close.groupby(['ALT', 'REF']).size() \
                      .astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .median().astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .quantile(q=0.25).astype(str).tolist()),
                      ', '.join(close.groupby(['ALT', 'REF'])['TARGET_VAF'] \
                      .quantile(q=0.75).astype(str).tolist()),
                      ', '.join(set(close['START'].astype(str))),
                      ', '.join(set([str(a) + '>' + str(b) for a, b in \
                      zip(close['REF'].tolist(), close['ALT'].tolist())])),
                      ', '.join(set(close['ANNOTATION']))
                     ]
            return pd.Series(cl_out)
        else:
            return pd.Series([None]*8)

    out_names = ["_Class", "_Frequency", "_VAF", "_Q25", "_Q75",
                 "_Position", "_Change", "_Annotation"]
    out_names = ['myTYPE' + s for s in out_names]

    variants[out_names] = variants.apply(lambda row: annot_row(row, mytype),
                                         axis=1)
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
def filter_call(variants):
    """
    Filter MFLAG_CALL: 1 if variant is not called by either caller
                          (i.e. caveman, strelka or mutect)
    """
    variants['MFLAG_CALL'] = np.where(variants['NUMBER_OF_CALLERS'] == 0, 1, 0)
    return(variants)

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
    variants['MFLAG_DROP'] = np.where(variants.VAG_GENE.isin(drop), 1, 0)
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
    Indels: Counts as "in cosmic" like for MAFCOS flag.
    SNV:    Only removes missense mutations with this flag
            Counts as "in cosmic" only if HEME_EXACT
    """
    if mode == 'snv':
        drop = ['non_synonymous_codon']
        variants['MFLAG_NONPASS'] = np.where(
            (variants['FLAGS_ALL'] != 'PASS') &
            (variants['VAG_EFFECT'].isin(drop)) &
            (variants['ANY_EXACT_POS'] == 0) &
            (variants['KNOWN_MM'] == 0), 1, 0)
        return(variants)
    variants['MFLAG_NONPASS'] = np.where(
        (variants['FLAGS_ALL'] != 'PASS') &
        (variants['COSMIC'].isnull()) &
        (variants['KNOWN_MM'] == 0), 1, 0)
    return(variants)

def filter_normals(variants):
    """
    Filter MFLAG_NORM: 1 if variant has genomic exact or pos in normals
    """
    match = ['genomic_exact', 'genomic_pos']
    variants['MFLAG_NORM'] = np.where(variants['Normals_Class'] \
                                      .isin(match), 1, 0)
    return(variants)

def filter_vaf(variants):
    """
    Filter MFLAG_VAF: 1 if non-pass and 
                    target VAF < 1 %, or 
                    reference VAF is > 1% or
                    VAF not estimated (i.e. 0).
    """
    variants['MFLAG_VAF'] = np.where(
        (variants['NUMBER_OF_CALLERS'] == 0) &
        (variants['TARGET_VAF'] < 0.01) & 
        (variants['REFERENCE_VAF'] > 0.01), 1, 0)
    return(variants)

def filter_dir(variants):
    """
    Filter MFLAG_DIR: 1 if DIRPROP = 0
    """
    variants['MFLAG_DIR'] = np.where(variants['DIRPROP'] == 0, 1, 0)
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
    stop = timeit.default_timer()

    with open(textname, 'w') as f:
        # Call the "Version" file for version info?
        f.write(
            f'Somatic variant processing for myTYPE\nv.1.0\n '
            f'Completed time: {str(datetime.today()).split(".")[0]}\n')
        f.write(f'Run time: {round(stop-START, 3)}\n')
        f.write(f'####\nMode: {mode}\n')
        f.write(f'Imported calls: {variants.shape[0]}\n')
        f.write('Flagging variants for filtering:\n')
        f.write(f'MFLAG_CALL: Not called by any callers i.e. caveman, '
                f'strelka or mutect: {variants["MFLAG_CALL"].sum()}\n')
        f.write(f'MFLAG_PANEL: Variant not in BED file of '
                f'regions to keep: {variants["MFLAG_PANEL"].sum()}\n')
        f.write(f'MFLAG_DROP: Variant in excluded gene: '
                f'{variants["MFLAG_DROP"].sum()}\n')
        f.write(f'MFLAG_IGH: In IGH locus: {variants["MFLAG_IGH"].sum()}\n')
        f.write(f'MFLAG_MAF: MAF > 3 % in normal databases: '
                f'{variants["MFLAG_MAF"].sum()}\n')
        f.write(f'MFLAG_MAFCOS: MAF > 0.1 % and not in COSMIC '
                f'(exact/pos): {variants["MFLAG_MAFCOS"].sum()}\n')
        f.write(f'MFLAG_NONPASS: NON-PASS non-synonymous mutations'
                f' IF not exact/pos cosmic or previously known in MM: '
                f'{variants["MFLAG_NONPASS"].sum()}\n')
        f.write(f'MFLAG_NORM: Variant exact or pos in >0 good normals: '
                f'{variants["MFLAG_NORM"].sum()}\n')
        f.write(f'MFLAG_VAF: Remove NON-PASS calls with target '
                f'VAF < 1 %, reference VAF > 1 % or '
                f'VAF not reported: {variants["MFLAG_VAF"].sum()}\n')
        f.write(f'MFLAG_DIR: Remove variants DIRPROP = 0 (only reads '
                f'on one strand): {variants["MFLAG_DIR"].sum(0)}\n')
        f.write(f'Removing calls with >= 1 MFLAG: {bad.shape[0]}\n')
        f.write(f'Calls passed filters: {good.shape[0]}\n')
    return()

# Main Function
def process(
        mode,
        infile,
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
    variants = import_variants(infile)

    ## ANNOTATIONS
    variants = annotate_cosmic(variants)
    #if genes:
        # Only runs if a path was passed to optional argument "gene"
    #    variants = annotate_genefreq(variants, genes)
        # Replace this with mutation frequency from MMRF? (and other raw data?)
    variants = annotate_maf(variants)
    variants = annotate_vaf_dep(variants)
    variants = annotate_normals(variants, normals)
    variants = annotate_mmrf(variants, mmrf)
    variants = annotate_bolli(variants, bolli)
    variants = annotate_lohr(variants, lohr)
    if mytype:
        # Only runs if a path was passed to optional argument "mytype"
        variants = annotate_mytype(variants, mytype)
    variants = annotate_known(variants, mytype)

    ## FILTERS
    variants = filter_call(variants)
    variants = filter_panel(variants, genes_bed)
    if genes_drop:
        variants = filter_drop(variants, genes_drop)
    variants = filter_igh(variants, igh)
    variants = filter_maf(variants)
    variants = filter_maf_cosmic(variants, mode) # Also include other historical cohorts than cosmic?
    variants = filter_nonpass(variants, mode)
    variants = filter_normals(variants)
    variants = filter_vaf(variants)
    variants = filter_dir(variants)

    ## OUTPUT
    name = namecore(infile)
    filter_export(variants, outdir, name, mode)
    print('Variant processing complete')
    return(variants) # Added this here - may be necessary for test?
