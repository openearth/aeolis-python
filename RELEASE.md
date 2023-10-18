# Making a new release of AeoLiS

We publish the package on PyPI. Only administrators can create new releases, and releases are made using the **main** branch.

## 1. Prepare Release
All packaging instructions are in the [pyproject.toml](pyproject.toml)

A. Before generating a new release, make sure that all  relevant code changes have been merge to the **main** branch.

B. Bump the version by updating the `pyproject.toml`. For consistency, you also have to update the version on `CITATION.cff` and `README.md` (citation section).
<br> (For pre-release versioning [see](https://packaging.python.org/en/latest/guides/distributing-packages-using-setuptools/#pre-release-versioning))

C. To test if the Python source package (`.tar.gz`) and the binary package (`.whl`) can be built without problem in the `dist/` directory, install `build` and build the package:

```bash
# install build
pip install build

# build package
python -m build
```

D. Document the changes in the new release by updating the `docs/whatsnew.rst` file. Document changes following the format of previous releases.

E. Make any relavant updates or changes to the package documentation. Then, push and merge changes to the `main` branch.

## 2. Publish Release

A. On GitHub, create a new release. Use `v.x.x.x` as **tag** and **title** for the new release. In the description of the release use the same text we used in `docs/whatsnew.rst`. Mind about the MarkDown formatting.

B. Then, publish the relaese. This will trigger the `Upload PyPI release` CI workflow which builds and uploads the package to PyPI. This will take some time to complete.

