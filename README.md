# research-analysis-tools

## Set of tools for scientific data analysis.

As of march 2026 only scripts and libraries for hopping transport analysis are added.

### hopping project structure

```
research-analysis-tools/
├─ hopping/
│  ├─ src/
│  │  ├─ database/
│  │  ├─ datafiles/
│  │  ├─ lib.py
│  │  ├─ TL10.py
│  ├─ output/
│  │  ├─ Plots in pdf or png format
│  ├─ Jupyter notebooks for analysis
```

* `database/` folder contains databases in form of SQLite and JSON,
* `datafiles/` conductance measurements in form of CSV files,
* `lib.py` library with useful functions,
* `TL10.py` module for managing database through **pandas** library.

Analysis is done through jupyter notebooks.

### TODO

* implement SQLite database located in `src/database/`, now database in form of JSON is used,
* add more jupyter notebooks to *hopping* project,
* add tools for *contactless electroreflectance (CER) spectrometry* and *UV-Vis spectrophotometry* data analysis.
