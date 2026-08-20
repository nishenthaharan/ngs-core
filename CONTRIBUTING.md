# Contributing

Thank you for helping improve NGS Core.

## Development setup

```bash
git clone https://github.com/nishenthaharan/ngs-core.git
cd ngs-core
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Before submitting a change, run:

```bash
ruff check .
ruff format --check .
pytest
python -m build
```

## Pull requests

- Keep each pull request focused on one behaviour or tightly related group of changes.
- Explain the biological or operational reason for changed defaults.
- Include tests for valid inputs, malformed inputs, and important boundary conditions.
- Update the README and schema documentation when the CLI or output changes.
- Do not commit identifiable human sequence data, production datasets, credentials, or
  large binary fixtures. Use small synthetic reads in tests.

## Scientific correctness

Changes that alter QC calculations or filtering semantics should include a worked example
or comparison with an established implementation. Platform- or assay-specific assumptions
must be explicit rather than encoded as universal defaults.
