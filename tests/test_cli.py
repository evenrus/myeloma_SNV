"""myeloma_snv cli tests."""
'''
from datetime import datetime
from os.path import isfile, join
from click.testing import CliRunner

from myeloma_snv import cli

def test_main_snv(tmpdir):
    """Test for main command."""

    outdir = str(tmpdir)

    params = []
    with open("../tests/test_args_snv.txt") as f:
        params = f.readlines()
    params = [x.strip() for x in params]

    params.append(join('--outdir', outdir))

    runner = CliRunner()

    result = runner.invoke(cli.main, params)
    assert result.exit_code == 0 # Exit code is 2!!

    date = str(datetime.today()).split()[0].split("-")
    name = '/ID131074'
    name = '_'.join([name, '_'.join(date)])
    expected_goodcalls_csv = join(outdir, name + '_goodcalls.csv')
    expected_badcalls_csv = join(outdir, name + '_badcalls.csv')
    expected_report_txt = join(outdir, name + '_report.txt')

    assert isfile(expected_goodcalls_csv)
    assert isfile(expected_badcalls_csv)
    assert isfile(expected_report_txt)

    # Something is not working here:
    #       Exit code is not 0
    #       outdir does not seem to work
    #
    # Check if output file exists. Pass the tempdir as outdir to the function.
    #
    # Can make environmental variables in pytest.ini folder.
    # This can include dictionaries with variables to pass.
'''
