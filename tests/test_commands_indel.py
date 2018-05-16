"""myeloma_snv commands tests."""

from myeloma_snv import commands

## Testing annotate functions for indels
def test_annotate_cosmic_indel(data_indel, ref_indel):
    """Test for annotation of heme exact in cosmic for indels"""
    test = commands.annotate_cosmic(data_indel)
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["HEME_EXACT"].fillna(0).astype('int').tolist()
    test = test["HEME_EXACT"].fillna(0).astype('int').tolist()
    assert ref == test

def test_annotate_genefreq_csv(data_indel, ref_indel, args):
    """Test for annotate gene frequency for indels"""
    test = commands.annotate_genefreq(data_indel,
                                      args['genes'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["MAX_MUTFREQ"].fillna(0).astype('float').round(1).tolist()
    test = test["MAX_MUTFREQ"].fillna(0).astype('float').round(1).tolist()
    assert ref == test

def test_annotate_maf_csv(data_indel, ref_indel):
    """Test for annotate maf for indels"""
    test = commands.annotate_maf(data_indel)
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref = ref["MAX_MAF"].fillna(0).astype('float').round(3).tolist()
    test = test["MAX_MAF"].fillna(0).astype('float').round(3).tolist()
    assert ref == test

def test_annotate_normals_indel(data_indel, ref_indel, args):
    """Test for good normal annotation of indels"""
    test = commands.annotate_normals(data_indel,
                                     args['normals_indel'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Normals_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["Normals_Frequency"].fillna(0)
    try:
        test_freq = test_freq.astype(int).astype('str').tolist()
    except ValueError:
        test_freq = test_freq.astype('str').tolist()
    ref_class = ref["Normals_Class"].fillna(0).tolist()
    test_class = test["Normals_Class"].fillna(0).tolist()
    assert ref_freq == test_freq
    assert ref_class == test_class

def test_annotate_mmrf_indel(data_indel, ref_indel, args):
    """Test for mmrf annotation of indels"""
    test = commands.annotate_mmrf(data_indel,
                                  args['mmrf'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["MMRF_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["MMRF_Frequency"].fillna(0)
    try:
        test_freq = test_freq.astype(int).astype('str').tolist()
    except ValueError:
        test_freq = test_freq.astype('str').tolist()
    ref_class = ref["MMRF_Class"].fillna(0).tolist()
    test_class = test["MMRF_Class"].fillna(0).tolist()
    assert ref_freq == test_freq
    assert ref_class == test_class

def test_annotate_bolli_indel(data_indel, ref_indel, args):
    """Test for bolli annotation of indels"""
    test = commands.annotate_bolli(data_indel,
                                   args['bolli'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Bolli_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["Bolli_Frequency"].fillna(0)
    try:
        test_freq = test_freq.astype(int).astype('str').tolist()
    except ValueError:
        test_freq = test_freq.astype('str').tolist()
    ref_annot = ref["Bolli_Annotation"].fillna(0).tolist()
    test_annot = test["Bolli_Annotation"].fillna(0).tolist()
    assert ref_freq == test_freq
    assert ref_annot == test_annot

def test_annotate_lohr_indel(data_indel, ref_indel, args):
    """Test for lohr annotation of indels"""
    test = commands.annotate_lohr(data_indel,
                                  args['lohr'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["Lohr_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["Lohr_Frequency"].fillna(0)
    try:
        test_freq = test_freq.astype(int).astype('str').tolist()
    except ValueError:
        test_freq = test_freq.astype('str').tolist()
    ref_class = ref["Lohr_Class"].fillna(0).tolist()
    test_class = test["Lohr_Class"].fillna(0).tolist()
    assert ref_freq == test_freq
    assert ref_class == test_class

def test_annotate_mytype_indel(data_indel, ref_indel, args):
    """Test for mytype annotation of indels"""
    test = commands.annotate_mytype(data_indel,
                                    args['mytype'])
    ref = ref_indel.sort_values('ID_VARIANT')
    test = test.sort_values('ID_VARIANT')
    ref_freq = ref["myTYPE_Frequency"].fillna(0).astype('str').tolist()
    test_freq = test["myTYPE_Frequency"].fillna(0)
    try:
        test_freq = test_freq.astype(int).astype('str').tolist()
    except ValueError:
        test_freq = test_freq.astype('str').tolist()
    ref_annot = ref["myTYPE_Annotation"].fillna(0).tolist()
    test_annot = test["myTYPE_Annotation"].fillna(0).tolist()
    assert ref_freq == test_freq
    assert ref_annot == test_annot

def test_annotate_known_indel(ref_indel, args):
    """Test for annotate known variants for indels"""
    test = ref_indel.drop('KNOWN_MM', axis=1)
    test = commands.annotate_known(test,
                                   args['mytype'])
    ref = ref_indel.sort_values('ID_VARIANT')["KNOWN_MM"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')["KNOWN_MM"].astype('int').tolist()
    assert ref == test

## Testing filters for indels

def test_filter_panel_indel(ref_indel, args):
    """Test for filter panel function for indels"""
    test = ref_indel.drop('MFLAG_PANEL', axis=1)
    test = commands.filter_panel(test,
                                 args['genes_bed'])
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_PANEL"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_PANEL"].astype('int').tolist()
    assert ref == test

def test_filter_drop_indel(ref_indel, args):
    """Test for filter drop function for indels"""
    test = ref_indel.drop('MFLAG_DROP', axis=1)
    test = commands.filter_drop(test,
                                args['genes_drop'])
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_DROP"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_DROP"].astype('int').tolist()
    assert ref == test

def test_filter_igh_indel(ref_indel, args):
    """Test for filter igh function for indels"""
    test = ref_indel.drop('MFLAG_IGH', axis=1)
    test = commands.filter_igh(test,
                               args['igh'])
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_IGH"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_IGH"].astype('int').tolist()
    assert ref == test

def test_filter_maf_indel(ref_indel):
    """Test for filter maf function for indels"""
    test = ref_indel.drop('MFLAG_MAF', axis=1)
    test['MAX_MAF'] = test['MAX_MAF'].astype('float')
    test = commands.filter_maf(test)
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_MAF"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_MAF"].astype('int').tolist()
    assert ref == test

def test_filter_maf_cosmic_indel(ref_indel, args):
    """Test for filter maf cosmic function for indels"""
    test = ref_indel.drop('MFLAG_MAFCOS', axis=1)
    test['MAX_MAF'] = test['MAX_MAF'].astype('float')
    test['ANY_EXACT_POS'] = test['ANY_EXACT_POS'].astype('float')
    test = commands.filter_maf_cosmic(test,
                                      args['mode_indel'])
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_MAFCOS"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_MAFCOS"].astype('int').tolist()
    assert ref == test

def test_filter_nonpass_indel(ref_indel, args):
    """Test for filter nonpass function for indels"""
    test = ref_indel.drop('MFLAG_NONPASS', axis=1)
    test['MAX_MAF'] = test['MAX_MAF'].astype('float')
    test['ANY_EXACT_POS'] = test['ANY_EXACT_POS'].astype('float')
    test['KNOWN_MM'] = test['KNOWN_MM'].astype('float')
    test = commands.filter_nonpass(test,
                                   args['mode_indel'])
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_NONPASS"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_NONPASS"].astype('int').tolist()
    assert ref == test

def test_filter_normals(ref_indel):
    """Test for filter normals function for indels"""
    test = ref_indel.drop('MFLAG_NORM', axis=1)
    test = commands.filter_normals(test)
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_NORM"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_NORM"].astype('int').tolist()
    assert ref == test

def test_filter_vaf(ref_indel):
    """Test for filter VAF function for indels"""
    test = ref_indel.drop('MFLAG_VAF', axis=1)
    test['TARGET_VAF'] = test['TARGET_VAF'].astype('float')
    test = commands.filter_vaf(test)
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_VAF"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_VAF"].astype('int').tolist()
    assert ref == test

def test_filter_bidir(ref_indel):
    """Test for filter bidir function for indels"""
    test = ref_indel.drop('MFLAG_BIDIR', axis=1)
    test['BIDIR'] = test['BIDIR'].astype('int')
    test = commands.filter_bidir(test)
    ref = ref_indel.sort_values('ID_VARIANT')
    ref = ref["MFLAG_BIDIR"].astype('int').tolist()
    test = test.sort_values('ID_VARIANT')
    test = test["MFLAG_BIDIR"].astype('int').tolist()
    assert ref == test
