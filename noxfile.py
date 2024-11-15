"""Nox sessions."""

import tempfile
from typing import Any

import nox
from nox.sessions import Session


package = "breakfast"
nox.options.sessions = "lint", "tests"
nox.options.default_venv_backend = "conda"
locations = "src", "tests", "./noxfile.py"


def install_with_constraints(session: Session, *args: str, **kwargs: Any) -> None:
    """Install packages constrained by Poetry's lock file.
    This function is a wrapper for nox.sessions.Session.install. It
    invokes pip to install packages inside of the session's virtualenv.
    Additionally, pip is passed a constraints file generated from
    Poetry's lock file, to ensure that the packages are pinned to the
    versions specified in poetry.lock. This allows you to manage the
    packages as Poetry development dependencies.
    Arguments:
        session: The Session object.
        args: Command-line arguments for pip.
        kwargs: Additional keyword arguments for Session.install.
    """
    with tempfile.NamedTemporaryFile() as requirements:
        # Hide a warning that I have already handled
        session.run(
            "poetry",
            "config",
            "warnings.export",
            "false",
            external=True,
        )
        session.run(
            "poetry",
            "export",
            "--with",
            "dev",
            "--format=requirements.txt",
            "--without-hashes",
            f"--output={requirements.name}",
            external=True,
        )
        session.install("-r", requirements.name, *args, **kwargs)


@nox.session(python="3.12")
def black(session: Session) -> None:
    """Run black code formatter."""
    args = session.posargs or locations
    install_with_constraints(session)
    session.run("black", *args)


@nox.session(python="3.12")
def lint(session: Session) -> None:
    """Lint using flake8."""
    args = session.posargs or locations
    install_with_constraints(
        session,
        "flake8",
        "flake8-annotations",
        "flake8-bandit",
        "flake8-black",
        "flake8-bugbear",
        "flake8-docstrings",
        "flake8-import-order",
        "darglint",
    )
    session.run("flake8", *args)


@nox.session(python=["3.11", "3.12"])
def tests(session):
    args = session.posargs or ["--cov"]
    session.run("poetry", "install", "--without", "dev", external=True)
    install_with_constraints(session)
    session.run("pytest", *args)
