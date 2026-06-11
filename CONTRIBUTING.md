# Welcome to the breakfast contributing guide <!-- omit in toc -->

In this guide you will get an overview of the contribution workflow from
setting up a development environment, testing your changes, submitting a pull
request, and performing a release.

## New contributor guide

To get an overview of the project itself, read the [README](README.md).

## Getting started

breakfast is written in Python and tries to follow the excellent packaging
guidelines ["Hypermodern Python" by Claudio
Jolowicz](https://cjolowicz.github.io/posts/hypermodern-python-01-setup/).
Most work on breakfast happens in an environment where administrator access is
not available, and the package should remain installable via
[conda](https://docs.conda.io/en/latest/index.html), [mamba](https://github.com/mamba-org/mamba),
and [bioconda](https://bioconda.github.io/).

### Setting up your development tools

Some tooling needs to be set up before you can work on breakfast. To install
this we use conda, and place the tools in their own environment:

```sh
conda env create -n breakfast-dev -f envs/breakfast-dev.yml
```

Then activate the environment:

```sh
conda activate breakfast-dev
```

The rest of this document assumes that you have the `breakfast-dev`
environment active. If the environment already exists and
`envs/breakfast-dev.yml` changed, update it with:

```sh
conda env update -n breakfast-dev -f envs/breakfast-dev.yml
```

After `just` is installed, the same update can be run as:

```sh
just env-update
```

### Installing the package

Install the project and development tools from the locked uv environment:

```sh
just sync
just --list
```

This runs `uv sync --locked --group dev`. If you prefer not to activate the
conda environment, run recipes through conda:

```sh
conda run -n breakfast-dev --live-stream just sync
```

### Testing

Run tests and quality checks with the local task runner:

```sh
just check
```

This runs the underlying uv commands:

```sh
uv run ruff format --check src tests
uv run ruff check src tests
uv run pytest --cov
```

If you prefer not to activate the conda environment, use:

```sh
conda run -n breakfast-dev --live-stream just check
```

CI runs tests on Python 3.11 and 3.12. If you need to test a specific Python
version locally, install it with uv and sync against it:

```sh
uv python install 3.11
uv sync --locked --group dev --python 3.11
just test
```

### Adding dependencies, updating dependency versions

Add runtime dependencies with uv:

```sh
uv add scikit-learn
uv add pandas
```

Add development dependencies to the `dev` group:

```sh
uv add --group dev pytest
```

Update the lockfile after changing dependency constraints:

```sh
uv lock
```

The project sets `tool.uv.exclude-newer = "1 week"` so uv does not lock very
fresh upstream releases immediately.

### Releasing a new version

Before starting a release, make sure the local checks and package build pass:

```sh
just check
just build
```

First update the version in `pyproject.toml` using uv:

```sh
uv version --bump patch
# <it will say the new version number here, e.g. 0.4.7>
git commit -am "Bump version"
git push
```

Then tag the commit with the same version number, using a `v` prefix, and push
the tag:

```sh
NEW_VERSION=v$(uv version --short)
git tag "$NEW_VERSION"
git push origin "$NEW_VERSION"
```

Now go to GitHub and create a release for that tag. Publishing the release
triggers CI, builds the package with `uv build`, and publishes it to PyPI with
`uv publish` if the tests pass.

If you followed the bioconda instructions above, the new PyPI package will be
detected automatically and trigger a GitHub pull request for the conda package
to be updated as well. After that is reviewed and approved, the bioconda
package will be updated. No interaction is necessary, although you might have
to do something extra if you change your package dependencies.

### Updating the Python version dependency

The supported Python range is declared in `pyproject.toml`:

```toml
[project]
requires-python = ">=3.11,<3.13"
```

After changing the range, refresh the lockfile and run the tests:

```sh
uv lock
just test
```

You might also need to update the version of Python in your conda development
environment.

### Updating the bioconda package when dependencies change

For package updates that do not add or remove dependencies, change dependency
version ranges, or change the allowed Python version, a normal release is
sufficient to automatically update both the PyPI and bioconda packages.
However, dependency or Python support changes may require updating the bioconda
`meta.yaml` file. This is explained in the
[bioconda docs](https://bioconda.github.io/contributor/updating.html).

TODO: document this process. I haven't had to do one of these updates yet.
