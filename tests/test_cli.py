"""myeloma_snv cli tests."""

from click.testing import CliRunner

from myeloma_snv import cli

def test_main_snv():
    """Test for main command."""
    params = []
    with open("../tests/test_args_snv.txt") as f:
        params = f.readlines()
    params = [x.strip() for x in params]

    runner = CliRunner()

    result = runner.invoke(cli.main, params)

    assert result.output
