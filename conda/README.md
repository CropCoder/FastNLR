# FastNLR conda recipe

This directory holds the Bioconda recipe for `fastnlr` and acts as the canonical copy of the
files submitted to [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).

Bioconda is the only conda channel FastNLR targets. The package is deliberately **not**
submitted to conda-forge: publishing the same package name to both channels lets users resolve
different builds depending on channel order, so there is a single canonical home.

## Contents

| File | Purpose |
|---|---|
| `fastnlr/meta.yaml` | recipe metadata: source URL, sha256, requirements, tests, license |
| `fastnlr/build.sh` | Unix build script: bundle third-party licenses, install the `fastnlr` binary |
| `fastnlr/bld.bat` | retained for reference only; Bioconda does not build Windows packages |

## Local build

```bash
conda install -y conda-build
conda build conda/fastnlr
```

## Bioconda submission

The published recipe lives in `recipes/fastnlr/` inside the `bioconda/bioconda-recipes`
repository and is kept in sync with the files here.

1. Fork `bioconda/bioconda-recipes`.
2. Create or update `recipes/fastnlr/`.
3. Copy `meta.yaml` and `build.sh` from this directory.
4. Confirm the `version` and `sha256` values against the release tarball.
5. Open a pull request and follow the Bioconda review checklist.

Recipe notes:

- `run_exports` is required for every Bioconda recipe, command-line tools included, and is
  enforced by the `missing_run_exports` lint.
- `build.sh` runs `cargo install --locked --bins --path crates/nlr-cli`, installing the
  `fastnlr` binary into `$PREFIX/bin`.
- `cargo-bundle-licenses` generates `THIRDPARTY.yml`, which is shipped as a license file and
  copied under `$PREFIX/licenses`.

## Release checklist

For every new tag:

1. Set `version` in `meta.yaml` to the tag without the leading `v`.
2. Recompute the hash:

   ```bash
   curl -L https://github.com/CropCoder/FastNLR/archive/refs/tags/v<version>.tar.gz | sha256sum
   ```

3. Reset `build: number` to `0`.
4. Apply the identical change to `bioconda/bioconda-recipes`.
