# Making a new release of AeoLiS

We publish the package on PyPI. Only administrators can create new releases.

## 1. Prepare Release
 All of the Python
packaging instructions are in the the following files:

* `aoelis/VERSION`
* `setup.cfg`
* `setup.py`

A. Before generating a package, we first need to install `build`.

```bash
pip install build 
```

B. Bump the version by updating the `aeolis/VERISON`. For consistency, you have to update the version on `CITATION.cff` and `README.md` (citation section).

C. To test if the Python source package (`.tar.gz`) and the binary package (`.whl`) can be built without problem in the `dist/` directory, do:

```bash
python -m build
```

D. Document the changes in the new release by updating the `docs/whatsnew.rst` file. Document changes following the format of previous releases.

E. Make any relavant updates or changes to the package documentation. Then, push and merge changes to the `main` branch.

## 2. Publish Release

A. On GitHub, create a new release. Use `v.x.x.x` as tag and title for the new release. In the description of the release use the same text we used in `docs/whatsnew.rst`. Mind about the MarkDown formatting.

C. Then, publish the relaese. This will trigger the `Upload PyPI release` workflow which builds and uploads the package to PyPI. This will take some time to complete.

