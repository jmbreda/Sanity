## [2.0.0] - 2026-07-29
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
- Fixed a sign error in the trigamma asymptotic series, which slightly affected the error bars on mu (and hence `d_mu.txt` and `ltq_error_bars.txt`) for genes with very low total counts, by under 1e-4 relative.
