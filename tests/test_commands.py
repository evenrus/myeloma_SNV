"""myeloma_snv commands tests."""

import pandas as pd

from myeloma_snv import commands

## Testing import function
def test_import_variants_csv(ref_snv, args):
    """Test for import variants csv format"""
    test = commands.import_variants(args['infile_csv'],
                                    int(args['skiplines_def']))
    ref = ref_snv["ID_VARIANT"].tolist().sort()
    test = test["ID_VARIANT"].tolist().sort()
    assert ref == test

def test_import_variants_gz(ref_snv, args):
    """Test for import variants gz format"""
    test = commands.import_variants(args['infile_gz'],
                                    int(args['skiplines_def']))
    ref = ref_snv["ID_VARIANT"].tolist().sort()
    test = test["ID_VARIANT"].tolist().sort()
    assert ref == test

## Testing annotate functions for SNVs
def test_annotate_cosmic_snv(data_snv, ref_snv):
    """Test for annotation of heme exact in cosmic"""
    test = commands.annotate_cosmic(data_snv)
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["HEME_EXACT"].fillna(0).astype('int').tolist()
    test = test["HEME_EXACT"].fillna(0).astype('int').tolist()
    assert ref == test

def test_annotate_genefreq_csv(data_snv, ref_snv, args):
    """Test for annotate gene frequency for SNVs"""
    test = commands.annotate_genefreq(data_snv,
                                    args['genes'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["MAX_MUTFREQ"].fillna(0).astype('float').round(1).tolist()
    test = test["MAX_MUTFREQ"].fillna(0).astype('float').round(1).tolist()
    assert ref == test

def test_annotate_maf_csv(data_snv, ref_snv):
    """Test for annotate maf for SNVs"""
    test = commands.annotate_maf(data_snv)
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["MAX_MAF"].fillna(0).astype('float').round(3).tolist()
    test = test["MAX_MAF"].fillna(0).astype('float').round(3).tolist()
    assert ref == test

def test_annotate_normals_snv(data_snv, ref_snv, args):
    """Test for good normal annotation of SNVs"""
    test = commands.annotate_normals(data_snv,
                                     args['normals'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Normals_Frequency"].astype('str').tolist()
    test_freq = test["Normals_Frequency"].astype('str').tolist()
    ref_vaf = ref["Normals_median_VAF"].astype('float').round(2).tolist()
    test_vaf = test["Normals_median_VAF"].astype('float').round(2).tolist()
    assert (ref_freq == test_freq) & (ref_vaf == test_vaf)

def test_annotate_mmrf_snv(data_snv, ref_snv, args):
    """Test for mmrf annotation of SNVs"""
    test = commands.annotate_mmrf(data_snv,
                                  args['mmrf'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["MMRF_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["MMRF_Frequency"].fillna(0).astype('str').tolist()
    ref_class = ref["MMRF_Class"].fillna(0).tolist()
    test_class = test["MMRF_Class"].fillna(0).tolist()
    assert (ref_freq == test_freq) & (ref_class == test_class)

def test_annotate_bolli_snv(data_snv, ref_snv, args):
    """Test for bolli annotation of SNVs"""
    test = commands.annotate_bolli(data_snv,
                                  args['bolli'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Bolli_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["Bolli_Frequency"].fillna(0).astype('str').tolist()
    ref_annot = ref["Bolli_Annotation"].fillna(0).tolist()
    test_annot = test["Bolli_Annotation"].fillna(0).tolist()
    assert (ref_freq == test_freq) & (ref_annot == test_annot)
    
def test_annotate_lohr_snv(data_snv, ref_snv, args):
    """Test for lohr annotation of SNVs"""
    test = commands.annotate_lohr(data_snv,
                                  args['lohr'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Lohr_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["Lohr_Frequency"].fillna(0).astype('str').tolist()
    ref_class = ref["Lohr_Class"].fillna(0).tolist()
    test_class = test["Lohr_Class"].fillna(0).tolist()
    assert (ref_freq == test_freq) & (ref_class == test_class)

def test_annotate_mytype_snv(data_snv, ref_snv, args):
    """Test for mytype annotation of SNVs"""
    test = commands.annotate_mytype(data_snv,
                                  args['mytype'])
    ref = ref_snv.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["myTYPE_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["myTYPE_Frequency"].fillna(0).astype('str').tolist()
    ref_annot = ref["myTYPE_Annotation"].fillna(0).tolist()
    test_annot = test["myTYPE_Annotation"].fillna(0).tolist()
    assert (ref_freq == test_freq) & (ref_annot == test_annot)

def test_annotate_known_snv(ref_snv, args):
    """Test for annotate known variants for SNVs"""
    test = ref_snv.drop('KNOWN_MM', axis=1)
    test = commands.annotate_known(test,
                                   args['mytype'])
    ref = ref_snv.sort_values('ID_VARIANT')["KNOWN_MM"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')["KNOWN_MM"].astype('int').tolist()
    assert ref == test

## Testing filters for SNVs

def test_filter_panel_snv(ref_snv, args):
    """Test for filter panel function for SNVs"""
    test = ref_snv.drop('MFLAG_PANEL', axis=1)
    test = commands.filter_panel(test,
                                 args['genes_bed'])
    ref = ref_snv.sort_values('ID_VARIANT')
    ref = ref["MFLAG_PANEL"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_PANEL"].astype('int').tolist()
    assert ref == test

def test_filter_igh_snv(ref_snv, args):
    """Test for filter igh function for SNVs"""
    test = ref_snv.drop('MFLAG_IGH', axis=1)
    test = commands.filter_igh(test,
                               args['igh'])
    ref = ref_snv.sort_values('ID_VARIANT')
    ref = ref["MFLAG_IGH"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_IGH"].astype('int').tolist()
    assert ref == test

