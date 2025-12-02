# Development Notes

- Doxygen XML is expected at `docs/doxygen/xml/` (set in `docs/Doxyfile`).
- Sphinx configuration lives in `docs/sphinx/conf.py`. Theme: Furo, with Breathe + MyST Markdown enabled.
- Build the site locally with:
  ```bash
  pip install -r docs/requirements.txt
  doxygen docs/Doxyfile
  sphinx-build -b html docs/sphinx docs/_build/html
  ```
- Use `sphinx-autobuild` (optional) for live reload during editing:
  ```bash
  pip install sphinx-autobuild
  sphinx-autobuild docs/sphinx docs/_build/html
  ```
