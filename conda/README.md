# FastNLR conda recipe

This directory contains a recipe suitable for local conda builds and for
submission to conda-forge or bioconda.

For a bioinformatics command-line tool, the usual primary channel is
**bioconda**. Submitting the same package name to both conda-forge and
bioconda is discouraged because users may resolve different builds from the
two channels. Choose one channel as the canonical home; the recipe files are
compatible with either.

## Local build

```bash
conda install -y conda-build
conda build conda/fastnlr
```

## conda-forge submission

1. Fork `conda-forge/staged-recipes`.
2. Create a directory `recipes/fastnlr`.
3. Copy `meta.yaml`, `build.sh`, and `bld.bat` from this directory.
4. Verify the `sha256` value in `meta.yaml` matches the final release tarball.
5. Open a pull request.

## bioconda submission

1. Fork `bioconda/bioconda-recipes`.
2. Create a directory `recipes/fastnlr`.
3. Copy the same recipe files.
4. Verify the version, source URL, and `sha256`.
5. Open a pull request and follow the bioconda review checklist.

If the release tag changes, update `version` and `sha256` in `meta.yaml`.
