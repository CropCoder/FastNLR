# Repository Guidelines

FastNLR is a Rust workspace for high-speed NLR immune-receptor locus annotation. Keep changes focused, tested, and compatible with the existing crate dependency graph.

## Project Structure & Module Organization

- `crates/` — Cargo workspace crates. Each crate has one responsibility; source lives in `src/`, integration tests in `tests/`.
- `crates/nlr-cli/` — CLI entry point and pipeline orchestration; `crates/nlr-cli/data/` contains the embedded motif profiles.
- `data/` — optional motif libraries, generation scripts, and their documentation.
- `conda/` — Bioconda packaging; `.github/workflows/` — CI and release automation.
- `docs/` — GitHub Pages content. `genome/`, `target/`, `/dist`, and runtime checkpoints are ignored by Git.

## Build, Test, and Development Commands

```bash
cargo build --workspace --locked   # build all crates with locked dependencies
cargo build --release              # optimized binary at target/release/fastnlr
cargo test --workspace --locked    # run all unit and integration tests
cargo fmt --all -- --check         # verify rustfmt compliance
cargo clippy --workspace --all-targets -- -D warnings  # lint with warnings as errors
```

Run locally with `./target/release/fastnlr -i genome.fasta -p out -t 8`.

## Coding Style & Naming Conventions

- Rust 2021 edition; use standard `rustfmt` formatting and four-space indentation.
- Keep the layered dependency direction: `nlr-core` depends on nothing, `nlr-cli` orchestrates the rest.
- Use `snake_case` for crates, modules, functions, and variables; `CamelCase` for types; `SCREAMING_SNAKE_CASE` for constants.
- Keep modules single-purpose, avoid unnecessary abstractions, and add doc comments for public items.

## Testing Guidelines

- Put unit tests in `#[cfg(test)]` modules; integration tests go under `crates/*/tests/`.
- Name integration files by behavior, such as `scan.rs`, `output.rs`, or `e2e.rs`.
- Cover sequence translation, scanning, assembly, output formats, and CLI end-to-end behavior.
- New bug fixes or features should include regression or integration tests, and all tests must pass before merge.

## Commit & Pull Request Guidelines

- Use conventional commit prefixes from existing history: `feat:`, `fix:`, `perf:`, `docs:`, `test:`, `build:`, `refactor:`, `cli:`.
- Write short, imperative subjects. Do not include agent or tooling metadata in commit messages.
- Target `main`, keep PRs scoped, and explain the change, rationale, and verification commands.
- Link related issues when applicable. Include screenshots only for output or plot changes.
