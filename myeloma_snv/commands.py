"""variants_process main command."""

import timeit
import re
from datetime import datetime
from os.path import join
import pandas as pd
import numpy as np
import pybedtools as pyb

pd.options.mode.chained_assignment = None  # default='warn'

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
    elif re.search('.tsv$', path):
        try:
            variants = pd.read_csv(
                filepath_or_buffer=path,
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
                        f'extension: try .csv, .tsv or .tsv.gz')

def determine_type(variants):
    len_var = variants['REF'].head(50)
    if len_var.str.len().mean() == 1:
        return('snv')
    return('indel')

## 1ST ANNOTATION
## Summaries of existing variables

def annotate_vaf_dep(variants):
    """
    Generate columns:
    TARGET_DEPTH: Mean TARGET_DEPTH from different callers
    REFERENCE_DEPTH: Mean REFERENCE_DEPTH from different callers
    """
    variants['TARGET_DEPTH'] = variants.filter(regex='TARGET_DEPTH').mean(axis=1)
    variants['TARGET_DEPTH'] = np.where(
        variants['TARGET_DEPTH'].isnull(), 0, variants['TARGET_DEPTH'])
    variants['REFERENCE_DEPTH'] = variants.filter(regex='REFERENCE_DEPTH').mean(axis=1)
    variants['REFERENCE_DEPTH'] = np.where(
        variants['REFERENCE_DEPTH'].isnull(), 0, variants['REFERENCE_DEPTH'])
    variants['DIRPROP'] = variants.filter(regex='DIRPROP').mean(axis=1)

    return(variants)

def annotate_maf(variants):
    """
    Generate column:
    MAX_MAF: Maximal MAF of variant in any normal database
    """
    variants['MAX_MAF'] = variants.filter(regex='_AF$|_AF_').max(axis=1)
    variants['MAX_MAF'] = np.where(
        variants['MAX_MAF'].isnull(), 0, variants['MAX_MAF'])
    return(variants)

## 1ST FILTERS
## Remove calls based on information present in input.

def filter_callers(variants):
    """
    Keep if passed by 1 or more callers.
    """
    vars_out = variants[variants['NUMBER_OF_CALLERS'] > 0]
    return(vars_out)

def filter_vartype(variants):
    """
    Keep coding variants (VAG_EFFECT)
    """

    # Variants considered as coding. Copied from leukgen/click_annotvcf 8.24.2018
    CODING_EFFECT = [
        "initiator_codon_change",
        "non_synonymous_codon",
        "splice_site_variant",
        "stop_gained",
        "stop_lost",
        "stop_retained_variant",
        "complex_change_in_transcript",
        "frameshift_variant",
        "inframe_codon_gain",
        "inframe_codon_loss",
        "inframe_variant"]

    vars_out = variants[variants["VAG_EFFECT"].isin(CODING_EFFECT)]

    return(vars_out)

def filter_vaf_dep(variants, vaf):
    """
    Keep if Target VAF > cut-off (default 2 %).
    Keep if Target Depth > 100
    """
    vars_out = variants[(variants['TARGET_VAF_MEAN'] > vaf) &
                        (variants['TARGET_DEPTH'] > 100)]
    return(vars_out)

def filter_maf(variants):
    """
    Keep if variant MAF < 0.05 % in exac/1000genomes
    Remove if between 0.025 and 0.05 %, if:
        Missense
        Not in BRCA1, BRCA2 or TP53
    """
    vars_out = variants[variants['MAX_MAF'] < 0.005]

    keepgenes = ['BRCA1', 'BRCA2', 'TP53']

    vars_out = vars_out[~((vars_out['MAX_MAF'] > 0.0025) & 
                          (vars_out['VAG_EFFECT'] == 'non_synonymous_codon') &
                          (~vars_out['VAG_GENE'].isin(keepgenes)))]
    return(vars_out)

def filter_qual(variants, mode):
    """
    Keep if reads on both strands
    Keep if 5 or more supporting reads
    Keep if no SE or SR flag
    """
    vars_out = variants[((variants['mutect_READS_FORWARD'].isnull()) |
                        (variants['mutect_READS_FORWARD'] >0)) &
                        ((variants['mutect_READS_REVERSE'].isnull()) |
                        (variants['mutect_READS_REVERSE'] > 0)) &
                        ((variants['mutect_READS_FORWARD'].isnull()) |
                        ((variants['mutect_READS_FORWARD'] + variants['mutect_READS_FORWARD']) > 4))]

    if mode == 'snv':
        vars_out = variants[(variants['FLAG_SE'] == 0) &
                            (variants['FLAG_SR'] == 0) &
                            ((variants['caveman_READS_FORWARD'].isnull()) |
                            (variants['caveman_READS_FORWARD'] >0)) &
                            ((variants['caveman_READS_REVERSE'].isnull()) |
                            (variants['caveman_READS_REVERSE'] >0)) &
                            ((variants['caveman_READS_FORWARD'].isnull()) |
                            ((variants['caveman_READS_FORWARD'] + variants['caveman_READS_FORWARD']) > 4))]
    return(vars_out)

def filter_panel(variants, genes_bed):
    """
    Remove if variant is not in BED file of regions to keep
    """
    # Ensure correct input data format
    variants["START"] = variants["START"].astype(np.int64)
    variants["END"] = variants["END"].astype(np.int64)
    # Generate bed-file from dataFrame
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
    
    vars_out = variants[variants.ID_VARIANT.isin(flaglist)]
    return(vars_out)

def filter_drop(variants, genes_drop):
    """
    Remove if variant is in list of genes to drop.
    """
    drop = pd.read_excel(io=genes_drop)['GENE']
    vars_out = variants[~variants.VAG_GENE.isin(drop)]
    return(vars_out)

def filter_igh(variants, igh_path):
    """
    Remove if variant in IGH locus
    """
    variants_bed = variants[["CHR", "START", "END", "ID_VARIANT"]]
    variants_bed = pyb.BedTool.from_dataframe(variants_bed)
    igh = pyb.BedTool(igh_path)
    variants_inter = variants_bed.intersect(igh, u=True)
    flaglist = []
    if not variants_inter.head(n=1, as_string=True) == '':
        flaglist = pyb.BedTool.to_dataframe(variants_inter)['name']
    
    vars_out = variants[~variants.ID_VARIANT.isin(flaglist)]
    return(vars_out)

def filter_normals(variants):
    """
    Remove if variant has genomic exact or pos in normals
    """
    match = ['genomic_exact', 'genomic_same_pos']
    vars_out = variants[~variants['NORMALS_Match_Class'].isin(match)]
    return(vars_out)

## 2ND ANNOTATION
## Final variant annotation

def annotate_cosmic(variants):
    """
    HEME_EXACT: Number of exact matches for hematopoietic and
    lymphoid tissue in cosmic.
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
    freqlist = freqlist[['GENE', 'MAX_MUTFREQ']].rename(
        index=str, columns={"GENE": "VAG_GENE"})
    variants = pd.merge(variants, freqlist, how='left')
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

## EXPORT
def namecore(infile):
    """
    Returns the "core" of the input file name, for use in output files.
    """
    name = infile.split('/')[-1]
    if re.search('.csv$', name):
        return(re.sub('.csv$', '', name))
    return(re.sub('.tsv.gz$', '', name))

def filter_export(allcalls, goodcalls, outdir, name, mode, vaf):
    """
    Function properties:
    1. Defines "good" and "bad" variants
    2. Writes files with good and bad variants.
    3. Creates processing summary report.
    """
    # Define badcalls
    badcalls = allcalls[~allcalls["ID_VARIANT"].isin(goodcalls["ID_VARIANT"])]

    # Define output names
    date = str(datetime.today()).split()[0].split("-")
    name = '_'.join([name, '_'.join(date)])
    goodname = join(outdir, name + '_goodcalls.csv')
    badname = join(outdir, name + '_badcalls.csv')
    textname = join(outdir, name + '_report.txt')

    # Export files
    goodcalls.to_csv(
        path_or_buf=goodname,
        index=False)
    badcalls.to_csv(
        path_or_buf=badname,
        index=False)

    # Summary report
    stop = timeit.default_timer()

    with open(textname, 'w') as f:
        f.write(
            f'Somatic variant processing for myTYPE\nv.1.0\n '
            f'Completed time: {str(datetime.today()).split(".")[0]}\n')
        f.write(f'Run time: {round(stop-START, 3)}\n')
        f.write(f'####\nMode: {mode}\n')
        f.write(f'VAF threshold: {vaf*100} %\n####\n')
        f.write(f'Imported calls: {allcalls.shape[0]}\n')
        f.write(f'Filtered out: {badcalls.shape[0]}\n')
        f.write(f'Passed: {goodcalls.shape[0]}\n')
    return()

# Main Function
def process(
        infile,
        outdir,
        vaf,
        genes,
        genes_drop,
        genes_bed,
        igh,
        mmrf,
        bolli,
        lohr,
        mytype):
    """Main function to process myTYPE SNV and indel output"""
    ## IMPORTING DATA
    variants = import_variants(infile)
    mode = determine_type(variants)
    ## 1ST ANNOTATIONS
    variants = annotate_maf(variants)
    variants = annotate_vaf_dep(variants)

    ## 1ST FILTERS
    filtering = filter_callers(variants)
    filtering = filter_vartype(filtering)
    filtering = filter_vaf_dep(filtering, vaf)
    filtering = filter_maf(filtering)
    # filtering = filter_qual(filtering, mode) - not using this filter until we can figure out strand support. 
    if genes_bed:
        filtering = filter_panel(filtering, genes_bed)
    if genes_drop:
        filtering = filter_drop(filtering, genes_drop)
    if igh:
        filtering = filter_igh(filtering, igh)
    
    if 'NORMALS_Match_Class' in filtering.columns:
        filtering = filter_normals(filtering)

    ## 2ND ANNOTATION
    ## Final variant annotation

    annotation = annotate_cosmic(filtering)

    if mmrf:
        annotation = annotate_mmrf(annotation, mmrf)
    if bolli:
        annotation = annotate_bolli(annotation, bolli)
    if lohr:
        annotation = annotate_lohr(annotation, lohr)
    if mytype:
        annotation = annotate_mytype(annotation, mytype)
    if genes:
        annotation = annotate_genefreq(annotation, genes)

    ## OUTPUT
    name = namecore(infile)
    filter_export(variants, annotation, outdir, name, mode, vaf)
    print('Processing complete')
    return(annotation) # Added this here - may be necessary for test?
