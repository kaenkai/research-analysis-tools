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
│  ├─ output/
│  │  ├─ Plots in pdf or png format
│  ├─ Jupyter notebooks for analysis
├─ playground/
```

* `hopping`: a project for electron transport analisys
  * `database/` databases in form of SQLite and JSON (depreciated),
  * `datafiles/` conductance measurements in form of CSV files,
  * `lib.py` library with many useful functions for analysis and `TL10` class for managing database through **pandas** library.
* `playground`: a project mainly used as a sandbox for learning EDA and ML.

Analysis is done through jupyter notebooks.

### TODO

* add more jupyter notebooks to *hopping* project,
* add tools for *contactless electroreflectance (CER) spectrometry* and *UV-Vis spectrophotometry* data analysis.
