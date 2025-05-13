# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.4] 13-05-2025

## Changed

- Improved README for publication
- Switched installer option from `hatch` to `uv`

## [0.6.3] 16-04-2025

## Changed

- Changed key names in params json file/schema to better align with output parameter designation

## [0.6.2] 12-04-2025

## Fixed

- Fixed sorting error when fermo and mzmine ion adduct annotations are mixed

## [0.6.1] 10-04-2025

## Added

- Retrieve library and ion identity annotation from mzmine3/4 file

## [0.6.0] 10-04-2025

## Changed

- Flattened parameter .json format
- Module logic: all modules must be specified explicitly, no module runs by default
- Removed all default values from data models to omit redundancy
- Spectral library now accepts a dir of files (e.g. multiple mgf files), not a single file
- Removed automatic timeout of long-running functions

## [0.5.1] 05-04-2025

## Fixed

- Fixed publishing to PyPI

## [0.5.0] 05-04-2025

## Added

- Implemented compatibility with mzmine4 output

## Changed

- QuantPhenotype: reimplemented, allowing for flexible FDR
- SpecLibMgfParser: implemented lazy loading to better deal with large libraries, filter out invalid spectra

## Fixed 

- Fixed bug in mzmine parser hanging on "NaN to integer"


## [0.4.3] 22-07-2024

### Fixed:

- PhenotypeManager: prevented Pearson calculation on constant or NaN-containing arrays
- GeneralParser: fixed error-handling on malformed input files.
- MS2DeepScoreNetworker: fixed MS2 spectra filtering for ms2deepscore algorithm

## [0.4.2] 16-06-2024

### Fixed

- Fixed bug in SummaryWriter: a nonexisting function was referenced, leading to premature exit of module.

## [0.4.1] 16-06-2024

### Fixed

- Versioning

## [0.4.0] 15-06-2024

### Removed

- [Breaking change] Removed MS2Query de novo annotation after observation of process instability (unforeseen process termination by system with SIGKILL (9))

## [0.3.3] 06-06-2024

### Fixed

- Added exception for fail of Sample Scores calculation in case of lack of spectral similarity networking data (e.g. no MS/MS data provided)

## [0.3.2] 05-06-2024

### Fixed

- Corrected erroneous assignment of 'True' to 'module_passed' for FragmentAnnotator and NeutralLossAnnotator if no MS/MS information was provided

## [0.3.1] 05-06-2024

### Changed

- Loosened typing restrictions for Feature and Sample object attributes: area and height (intensity) now accept float values.

### Removed

- [Breaking change] Removed toggle 'nonbiological' from 'FragmentAnnotator' and from parameters file; 'nonbiological' fragment annotation is now performed automatically

### Security

- Instead of full file paths, only filenames are now written to the 'out.fermo.session.json' file

## [0.3.0] 03-06-2024

### Changed

- [Breaking change] Parameter settings for `additional_modules/feature_filtering` were changed from a list of ranges to a dictionary with explicit values.
- Reworked score assignment for qualitative phenotype data: phenotype-associated features now always receive a score of 1.0, and non-associated ones a score of 0.0.
- For all modules with runtime restriction, the 'maximum_runtime' parameter was set to a default of '0' (unlimited runtime). Therefore, runtime restriction must now be specified explicitly.
- Added a 'module_passed' parameter to all modules. This allows a more accurate description via the SummaryWriter (e.g. module was activated but timed out, and lack of e.g. annotation is due to premature ending and not because there were no hits).

## [0.2.2] 27-05-2024

### Changed

- Removed 'phenotypes' as separate Feature attribute: write to Annotation object instead.
- Implemented sorting of annotation entries in descending order

## [0.2.1] 26-05-2024

### Fixed

- Fixed bug in "SummaryWriter": implemented error catching.

## [0.2.0] 26-05-2024

### Changed

- Replaced global logger with logger specific for `main_cli()`. `main()` now needs an argument `logger`
- Reworked output file naming: all output files now start with `out.fermo.` and a suffix specifying their type
- Removed output directory selection: the output directory is now always `results` located in the directory in which the peaktable resides.
- Features now always have default result values (previously, some Features could have an empty dictionary)
- MS2Query assignment now uses temporary directories for data reading/writing

### Fixed

- Removed `sys.exit(0)` in case of successful run to fix compatibility issue with Celery task manager (`fermo_gui`)

## [0.1.1] 24-05-2024

### Changed

- Implemented logging verboseness parameter for CLI interface

## [0.1.0] 19-05-2024

- First public release of `fermo_core`