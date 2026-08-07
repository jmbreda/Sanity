# Contributing to Sanity

Thank you for considering contributing to Sanity!

## Code of Conduct

This project and everyone participating in it is governed by our
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## How can I contribute?

### Reporting bugs

Please check the existing issues first to avoid duplicates. A useful bug report includes:

- A clear, descriptive title
- The exact command line you ran, and the input it ran on
- What you expected to happen, and what happened instead
- Environment details: operating system, compiler version, and the output of `./bin/Sanity -v`

### Suggesting features

Feature requests are welcome. It helps to describe the problem you are trying to solve, not only
the solution you have in mind, along with any alternatives you considered.

### Improving documentation

Documentation improvements are always welcome: fixing typos, adding examples, or clarifying
sections that were confusing when you first read them.

### Submitting code

Look for issues labelled `good first issue` or `help wanted` for places to start.

## Development setup

### Prerequisites

- A C++ compiler with **C++17** support (the Makefile builds with `-std=c++17`, and the code uses
  `<filesystem>` and structured bindings, so this is a hard requirement)
- `make`
- Git
- OpenMP
- zlib
- Python 3 with `numpy` (`pip3 install numpy`), to run `tests/compare.py`

For how to install the OpenMP and zlib dependencies on Linux and macOS, follow the
[Installation section of the README](README.md#installation) rather than duplicating the steps
here — note in particular that on macOS Apple's `clang++` does not support OpenMP, so a real GCC
is required.

### Getting started

```bash
# 1. Fork the repository on GitHub, then clone your fork
git clone https://github.com/<your-username>/Sanity.git
cd Sanity

# 2. Add the upstream remote
git remote add upstream https://github.com/jmbreda/Sanity.git

# 3. Create a branch for your changes
git checkout -b feature/your-feature-name

# 4. Make your changes, then compile
cd src
make clean
make

# 5. Run the tests
cd ../tests
python3 compare.py
```

## Pull request process

1. Update your branch with the latest upstream changes:

   ```bash
   git fetch upstream
   git rebase upstream/master
   ```

2. Rebuild and re-run the tests, and confirm they still pass:

   ```bash
   cd src && make clean && make
   cd ../tests && python3 compare.py
   ```

3. Update the documentation if you changed the command-line options or added a feature.

4. Push your branch to your fork and open a pull request against the `master` branch. There is no
   pull request template; a clear description of what changed and why is enough.

### Checklist

- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new compiler warnings
- [ ] `python3 compare.py` passes for all four methods

## Style guide

### Commit messages

Please write clear, descriptive commit messages: a short summary line saying what changed, and a
body explaining why when the reason is not obvious. Sanity does not follow a formal commit message
convention.

### Code style

Sanity does not currently follow a single standard code style. Please match the style of the file
you are editing, and keep formatting changes out of commits that also change behaviour — mixing
the two makes a change very hard to review.

### Testing

Sanity does not have a comprehensive test suite. `tests/compare.py` checks the output of a run
against checked-in reference outputs for all four `-v_m` methods; see [tests/README.md](tests/README.md).
Any change that is meant to be behaviour-preserving should leave those outputs identical. If you
add a feature or fix a bug, adding a test is very welcome.

## Community

- [Discussions](https://github.com/jmbreda/Sanity/discussions) — ask questions

---

Thank you for contributing!
