# Documentation & CI

## Local doc build
```bash
python -m venv venv
source venv/bin/activate
pip install -r docs/requirements.txt

# Generate Doxygen XML/HTML and build the Sphinx site
doxygen docs/Doxyfile
sphinx-build -b html docs/sphinx docs/_build/html
```
Open `docs/_build/html/index.html` to view the site locally. The raw Doxygen HTML remains in `docs/doxygen/html/` if you need it.

## GitHub Pages via Actions
The repository includes a workflow (`.github/workflows/docs.yml`) that:
1. Installs Doxygen and Python doc dependencies.
2. Runs `doxygen docs/Doxyfile` (XML + HTML).
3. Builds the Sphinx site (`sphinx-build -b html docs/sphinx docs/_build/html`).
4. Publishes `docs/_build/html/` to GitHub Pages.

To test this from a fork on a branch (e.g., `docs-playground`):
1. Push your docs branch.
2. Ensure repository Settings → Pages → Source is set to “GitHub Actions.”
3. Check the workflow run for the published URL (typically `https://<username>.github.io/qasmtrans/`).
