# Welcome to the breakfast contributing guide <!-- omit in toc -->

In this guide you will get an overview of the contribution workflow from setting up a development environment, testing your changes, submitting a pull request and performing a release.

Use the table of contents icon on the top left corner of this document to get to a specific section of this guide quickly.

## New contributor guide

To get an overview of the project itself, read the [README](README.md).

## Getting started

breakfast is written in Python and tries to follow the excellent packaging guidelines ["Hypermodern Python" by Claudio Jolowicz](https://cjolowicz.github.io/posts/hypermodern-python-01-setup/). Nevertheless, there are some places where breakfast differs from those guidelines, and we have tried to outline those differences here wherever relevant. The main differences are caused by most work on breakfast happening in an environment where administrator access is not available (a shared Linux HPC), and also because we want our package to be installable via [conda](https://docs.conda.io/en/latest/index.html) or [mamba](https://github.com/mamba-org/mamba), from the [bioconda](https://bioconda.github.io/) channel in particular.

### Setting up your development tools

Some tooling needs to be set up before you can work on breakfast. To install this we use conda, and place them in their own environment:

```sh
conda env create -n breakfast-dev -f envs/breakfast-dev.yml
```

Then when you want to work on the project, or at the very least if you want to use poetry commands or run tests, you need to switch to this environment:

```sh
conda activate breakfast-dev
```

The rest of this document assumes that you have the breakfast-dev environment active.

### Installing the package

As you're developing, you can install what you have developed using poetry install into your breakfast-dev conda environment:

```sh
poetry install
breakfast --version
```

### Testing

Some basic tests have been implemented. They can be run using nox:

```sh
nox
```

This will run the tests in the `tests/` folder, as well as linting using flake8 and code formatting checks using black.

See Hypermodern Python [part 2](https://cjolowicz.github.io/posts/hypermodern-python-02-testing/) and [part 3](https://cjolowicz.github.io/posts/hypermodern-python-03-linting/) for more details.

### Adding dependencies, updating dependency versions

You can add dependencies using poetry:

```sh
poetry add scikit-learn
poetry add pandas
```

You can automatically update the dependency to the newest minor or patch release like this:

```sh
poetry update pandas
```

and for major releases you have to be more explicit, assuming you're coming from 1.x to 2.x:

```sh
poetry update pandas^2.0
```

### Releasing a new version

First update the version in pyproject.toml using poetry:

```sh
poetry version patch
# <it will say the new version number here, e.g. 0.3.1>
git commit -am "Bump version"
git push
```

Then tag the commit with the same version number (note the "v" prefix), push the code and push the tag:

```sh
NEW_VERSION=v$(poetry version -s)
git tag "$NEW_VERSION"
git push origin "$NEW_VERSION"
```

Now go to github.com and do a release, selecting the version number tag you just pushed. This will automatically trigger the new version being tested and pushed to PyPI if the tests pass.

If you followed the bioconda instructions above, the new pypi package will be detected automatically and trigger a GitHub pull request for the conda package to be updated as well. After that is reviewed and approved, the bioconda package will be updated. No interaction is necessary, although you might have to do something extra if you change your package dependencies.

### Updating the python version dependency

Aside from updating package dependencies, it is also sometimes useful to update the dependency on python itself. One way to do this is to edit the pyproject.toml file and change the python version description. Versions can be specified using constraints that are documented in the [poetry docs](https://python-poetry.org/docs/dependency-specification/):

```
[tool.poetry.dependencies]
python = "^3.10"  # <-- this
```

Afterwards, you need to use poetry to update the poetry.lock file to reflect the change that you just made to the  pyproject.toml file. Be sure to use the `--no-update` flag to not update the locked versions of all dependency packages.

```sh
poetry lock --no-update
```

Then you need to run your tests to make sure everything is working, commit and push the changes.

You might also need to update/change the version of python in your conda environment, but I'm not certain about that.

### Updating the bioconda package when dependencies, dependency versions, or the python version has been changed

For package updates that don't lead to added/removed dependencies, changes to dependency versions, or changes to the allowed python version, a normal release (as above) is sufficient to automatically update both the PyPI and bioconda packages. However, for changes that do result in changes to dependencies it is necessary to update the bioconda meta.yml file. This is explained in [bioconda docs](https://bioconda.github.io/contributor/updating.html), and they also provide tools to help you with this.

TODO: document this process. I haven't had to do one of these updates yet.

### Conda and poetry virtual environment problems

When you first run poetry, it will detect some virtual environments on the system. In our case it should detect the conda env you currently have active and set that as the default for poetry. You can check that:

```
$ poetry env info

Virtualenv
Python:         3.10.4
Implementation: CPython
Path:           /home/<redacted>/.conda/envs/breakfast-dev
Valid:          True

System
Platform: linux
OS:       posix
Python:   /home/<redacted>/.conda/envs/breakfast-dev
```

If you decide to rename your conda development environment or have multiple projects and decide to use multiple conda environments, then you might have to switch the environment that poetry is using. This can be done by running the commands:

```
conda activate my-new-dev-env
poetry env use "$CONDA_PREFIX"
```

Note that the new environment needs to have poetry installed in this case.
