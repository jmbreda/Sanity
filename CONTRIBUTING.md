# Contributing to Sanity

First off, thank you for considering contributing to Sanity! It's 
people like you that make Sanity such a great tool.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Pull Request Process](#pull-request-process)
- [Style Guide](#style-guide)
- [Community](#community)

## Code of Conduct

This project and everyone participating in it is governed by our 
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are 
expected to uphold this code.

## How Can I Contribute?

### 🐛 Reporting Bugs

Before creating bug reports, please check existing issues to avoid 
duplicates. When you create a bug report, include as many details as 
possible.

**Great bug reports include:**
- A clear, descriptive title
- Steps to reproduce the behavior
- Expected behavior vs actual behavior
- Screenshots (if applicable)
- Environment details (OS, browser, version)

### 💡 Suggesting Features

Feature requests are welcome!

**Great feature requests include:**
- Clear problem statement: "I'm frustrated when..."
- Proposed solution
- Alternative solutions you've considered
- Additional context

### 📝 Improving Documentation

Documentation improvements are always welcome! This includes:
- Fixing typos
- Adding examples
- Clarifying confusing sections
- Translating documentation

### 🔧 Submitting Code

Look for issues labeled `good first issue` or `help wanted` for 
great places to start.

## Development Setup

### Prerequisites

- C++ compiler with C++11 support
- make
- Git
- libomp
- zlib

### Getting Started

```bash
# 1. Fork the repository on GitHub

# 2. Clone your fork locally
git clone https://github.com/breda/Sanity.git
cd Sanity

# 3. Add upstream remote
git remote add upstream https://github.com/jmbreda/Sanity.git

# 4. Create a branch for your changes
git checkout -b feature/your-feature-name

# 5. Code your changes
 ...

# Compile the project
cd src
make clean
make

# 7. Run tests to verify your changes
cd tests
python compare.py
```

## Pull Request Process

### Before Submitting

1. **Update your branch** with the latest upstream changes:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Run provided simple tests** and ensure all tests pass:
   ```bash
    cd tests
    python compare.py
   ```

3. **Update documentation** if you've changed APIs or added features.

### Submitting

1. Push your branch to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

2. Open a Pull Request against the `main` branch.

3. Fill out the PR template completely.

4. Wait for review.

### PR Checklist

- [ ] My code follows the project's style guidelines
- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes

## Style Guide

### Commit Messages

We follow [Conventional Commits](https://conventionalcommits.org/):

```
(): 

[optional body]

[optional footer]
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation only
- `style`: Formatting, missing semicolons, etc.
- `refactor`: Code change that neither fixes a bug nor adds a feature
- `test`: Adding missing tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(auth): add OAuth2 support
fix(api): handle null response from payment provider
docs(readme): update installation instructions
```

### Code Style

- While currently Sanity is not following any standard code style, we encourage contributors to follow a consistent style within their own contributions. 

### Testing

- Currentlty, Sanity does not have a comprehensive test suite but rather a simple test script to verify that results are matching the correct ones. However, we encourage contributors to add tests for new features and bug fixes.

## Community

- [Discussions](https://github.com/breda/Sanity/discussions) - Ask questions

## Recognition

Contributors are added to our [CONTRIBUTORS.md](CONTRIBUTORS.md) file 
and mentioned in release notes for significant contributions.

---

Thank you for contributing! 🎉