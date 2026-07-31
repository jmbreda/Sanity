## [2.0] - 2026-07-29
### Added
- New methods for calculating v_g (`-v_m` option):
    - MAP (**default**) maximum a posteriori estimate of v_g
    - EAP expected a posteriori estimate of v_g
    - MLE maximum likelihood estimate of v_g
    - MARG marginal likelihood estimate of v_g
- Lazy reading of the input file, and lazy output writing which significantly reduces memory usage for large datasets.
- Support for gzipped input files (`.gz` extension).
- Timestamped run command is now saved by default in the output directory (`sanity_command.txt`).

### Changed
- By default, *Sanity* will use the MAP method for calculating v_g, if not specified otherwise.
- No `*_vmax` output files. For all methods the output files are named the same.
