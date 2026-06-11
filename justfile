export UV_CACHE_DIR := ".uv-cache"

default:
    @just --list

sync:
    uv sync --locked --group dev

format:
    uv run ruff format src tests

lint:
    uv run ruff format --check src tests
    uv run ruff check src tests

test:
    uv run pytest --cov

check: lint test

build:
    uv build

env-update:
    conda env update -n breakfast-dev -f envs/breakfast-dev.yml

clean:
    rm -rf dist .coverage .pytest_cache .ruff_cache .uv-cache .venv src/*.egg-info
