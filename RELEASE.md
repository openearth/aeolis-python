# Making a new release of AeoLiS

The package is published to `PyPI`

## Manual release
 All of the Python
packaging instructions are in the the following files:

* `aoelis/VERSION`
* `setup.cfg`
* `setup.py`

Before generating a package, we first need to install `build`.

```bash
pip install build 
```

Bump the version using by updating the `aeolis/VERISON`. For consistency, you have to update the version on `CITATION.cff` and `README.md` (citation section).

To create a Python source package (`.tar.gz`) and the binary package (`.whl`) in the `dist/` directory, do:

```bash
python -m build
```

Document the changes in the new version by updating the `docs/whatsnew.rst` file. Document changes following format of previous releases.

Push and merge changes to the `main` branch.

## Publishing 

On GitHub, create a new release. Use `v.x.x.x` as title and tag for the new release, and in the description of the release we use the same text we used in `docs/whatsnew.rst`. Mind about the MarkDown formatting.

Then, publish the relaese. This will trigger the `Upload PyPI release` workflow which builds and uploads the package to PyPI. This will take some time to complete.




