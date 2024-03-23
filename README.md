npd
===========


[![](https://travis-ci.org/LostFan123/npd.svg?branch=master)](https://travis-ci.org/LostFan123/npd "Travis CI")
[![](https://dev.azure.com/skorobogatov/npd/_apis/build/status/LostFan123.npd?branchName=master)](https://dev.azure.com/skorobogatov/npd/_build/latest?definitionId=2&branchName=master "Azure Pipelines")
[![](https://codecov.io/gh/LostFan123/npd/branch/master/graph/badge.svg)](https://codecov.io/gh/LostFan123/npd "Codecov")
[![](https://img.shields.io/github/license/LostFan123/npd.svg)](https://github.com/LostFan123/npd/blob/master/LICENSE "License")
[![](https://badge.fury.io/py/npd.svg)](https://badge.fury.io/py/npd "PyPI")

Summary
-------

Npd (stands for **N**on-convex **p**olygon **d**ecomposition) is a Python library that implements an algorithm for
non-convex polygon decomposition into separate parts depending on the area 
requirements.

---

In what follows, `python` is an alias for `python3.8` or any later
version.

Installation
------------

Install the latest `pip` & `setuptools` packages versions:
  ```bash
  python -m pip install --upgrade pip setuptools
  ```

### User

Download and install the latest stable version from `PyPI` repository:
  ```bash
  python -m pip install --upgrade npd
  ```

### Developer

Download the latest version from `GitHub` repository
```bash
git clone https://github.com/LostFan123/npd.git
cd npd
```

Install:
  ```bash
  poetry install
  ```

Development
-----------

### Bumping version

#### Preparation

Install
[bump-my-version](https://github.com/callowayproject/bump-my-version/tree/master?tab=readme-ov-file#installation).

#### Pre-release

Choose which version number category to bump following [semver
specification](http://semver.org/).

Test bumping version
```bash
bump-my-version bump --dry-run --verbose $CATEGORY
```

where `$CATEGORY` is the target version number category name, possible
values are `patch`/`minor`/`major`.

Bump version
```bash
bump-my-version bump --verbose $CATEGORY
```

This will set version to `major.minor.patch-alpha`. 

#### Release

Test bumping version
```bash
bump-my-version bump --dry-run --verbose release
```

Bump version
```bash
bump-my-version bump --verbose release
```

This will set version to `major.minor.patch`.

#### Notes

To avoid inconsistency between branches and pull requests,
bumping version should be merged into `master` branch 
as separate pull request.

### Running tests

Plain:
  ```bash
  pytest
  ```

Inside `Docker` container:
  ```bash
  docker-compose run --entrypoint pytest npd-cpython
  ```
