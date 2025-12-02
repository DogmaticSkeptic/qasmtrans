# API Reference

The C++ API reference is generated with Doxygen and rendered by Sphinx via Breathe.

## Generate locally
```bash
# From the repository root
doxygen docs/Doxyfile
sphinx-build -b html docs/sphinx docs/_build/html
```
Open `docs/_build/html/index.html` in your browser to view the site (API is under “API Reference”).

## Outputs
- Doxygen XML: `docs/doxygen/xml/` (consumed by Breathe).
- Doxygen HTML: `docs/doxygen/html/` (optional raw Doxygen output).
- Sphinx HTML site: `docs/_build/html/` (what GitHub Pages publishes).
