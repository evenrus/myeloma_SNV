import pytest
import pandas as pd

@pytest.fixture()
def args():
    args = {}
    with open("../tests/test_args_dict.txt") as f:
        for line in f:
            (key, val) = line.split()
            args[key] = val
    return(args)

@pytest.fixture()
def data_snv(args):
    data = pd.read_csv(
        filepath_or_buffer=args['data_snv'],
        low_memory=False,
        na_values='NA')
    return(data)

@pytest.fixture()
def ref_snv(args):
    reference = pd.read_csv(
        filepath_or_buffer=args['reference_snv'],
        low_memory=False,
        na_values='NA',
        dtype = 'object')
    return(reference)

@pytest.fixture()
def data_indel(args):
    data = pd.read_csv(
        filepath_or_buffer=args['data_indel'],
        low_memory=False,
        na_values='NA')
    return(data)

@pytest.fixture()
def ref_indel(args):
    reference = pd.read_csv(
        filepath_or_buffer=args['reference_indel'],
        low_memory=False,
        na_values='NA',
        dtype = 'object')
    return(reference)