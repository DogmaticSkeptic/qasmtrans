# API Reference

The C++ API reference is generated with Doxygen from the sources under `include/` and `src/`.

## Generate locally
```bash
# From the repository root
doxygen docs/Doxyfile
mkdocs serve
```
After running Doxygen, open `http://127.0.0.1:8000/reference/html/index.html` while `mkdocs serve` is running.

## What gets generated
- XML output: `docs/reference/xml/` (for future integration with MkDocs plugins if desired).
- HTML output: `docs/reference/html/`, which is copied verbatim into the MkDocs site so it is available under `/reference/html/` on GitHub Pages.

If the API section appears empty on the published site, make sure the Doxygen step ran before `mkdocs build`.
