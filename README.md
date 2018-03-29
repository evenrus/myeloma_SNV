# myeloma_snv

[![pypi badge][pypi_badge]][pypi_base]
[![travis badge][travis_badge]][travis_base]
[![codecov badge][codecov_badge]][codecov_base]
[![docker badge][docker_badge]][docker_base]
[![docker badge][automated_badge]][docker_base]

Program for post-processing of SNV data from myTYPE with myeloma-specific annotations and filters

## Usage

Run with [containers][docker_base]:

        # docker usage
        docker run evenrus/myeloma_snv --help

        # singularity usage
        singularity run docker://evenrus/myeloma_snv --help

See [docker2singularity] if you want to use a [singularity] image instead of using the `docker://` prefix.

## Contributing

Contributions are welcome, and they are greatly appreciated, check our [contributing guidelines](.github/CONTRIBUTING.md)!

## Credits

This package was created using [Cookiecutter] and the
[leukgen/cookiecutter-toil] project template.

<!-- References -->
[singularity]: http://singularity.lbl.gov/
[docker2singularity]: https://github.com/singularityware/docker2singularity
[cookiecutter]: https://github.com/audreyr/cookiecutter
[leukgen/cookiecutter-toil]: https://github.com/leukgen/cookiecutter-toil
[`--batchSystem`]: http://toil.readthedocs.io/en/latest/developingWorkflows/batchSystem.html?highlight=BatchSystem

<!-- Badges -->
[docker_base]: https://hub.docker.com/r/evenrus/myeloma_snv
[docker_badge]: https://img.shields.io/docker/build/evenrus/myeloma_snv.svg
[automated_badge]: https://img.shields.io/docker/automated/leukgen/myeloma_snv.svg
[codecov_badge]: https://codecov.io/gh/evenrus/myeloma_snv/branch/master/graph/badge.svg
[codecov_base]: https://codecov.io/gh/evenrus/myeloma_snv
[pypi_badge]: https://img.shields.io/pypi/v/myeloma_snv.svg
[pypi_base]: https://pypi.python.org/pypi/myeloma_snv
[travis_badge]: https://img.shields.io/travis/evenrus/myeloma_snv.svg
[travis_base]: https://travis-ci.org/evenrus/myeloma_snv
