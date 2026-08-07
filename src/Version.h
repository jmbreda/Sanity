#ifndef _Version_h_
#define _Version_h_

// Single source of truth for the version reported by Sanity, Sanity_distance and
// Sanity_gene_correlation. The version also appears in Dockerfile and CHANGELOG.md,
// which are not generated from this header and must be changed by hand each release.
constexpr char SANITY_VERSION[] = "2.0.0";

#endif
