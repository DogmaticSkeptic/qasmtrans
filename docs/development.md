# Documentation & CI

## Local doc build
```bash
python -m venv venv
source venv/bin/activate
pip install -r docs/requirements.txt

# Generate Doxygen HTML/XML and serve the MkDocs site
doxygen docs/Doxyfile
mkdocs serve
```
The Doxygen HTML will be available under `reference/html/` inside the MkDocs site.

## GitHub Pages via Actions
The repository includes a workflow (`.github/workflows/docs.yml`) that:
1. Installs Doxygen and MkDocs dependencies.
2. Runs `doxygen docs/Doxyfile`.
3. Builds the MkDocs site (`mkdocs build`), which copies the generated Doxygen HTML.
4. Publishes the `site/` directory to GitHub Pages.

To test this from a fork on a branch (e.g., `docs-playground`):
1. Push your docs branch.
2. Ensure repository Settings → Pages → Source is set to “GitHub Actions.”
3. Check the workflow run for the published URL (typically `https://<username>.github.io/qasmtrans/`).
