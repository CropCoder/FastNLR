# FastNLR Desktop

Tauri 2 + React desktop application for genome and protein NLR annotation. The backend calls
`nlr-service` directly and does not start a local HTTP server.

## Supported targets

- Windows 10/11 x64: NSIS installer.
- macOS 12+ Apple Silicon: DMG installer.
- Linux: development and CI compile target only; the published Linux artifact is the CLI.

The current unsigned builds are intended for evaluation. Public warning-free distribution requires
Windows code signing and Apple Developer ID signing plus notarization.

## Requirements

- Node.js 22+
- Rust 1.90+
- Windows: WebView2
- macOS: Xcode command-line tools
- Ubuntu/Debian:

```bash
sudo apt install libwebkit2gtk-4.1-dev librsvg2-dev \
  patchelf libdbus-1-dev pkg-config
```

## Development

```bash
cd app
npm ci
npm run dev
```

Run the native desktop shell:

```bash
npm run tauri dev
```

## Validation

```bash
npm run check:versions
npm run typecheck
npm test
npm run build
cargo test --manifest-path src-tauri/Cargo.toml --locked
cargo clippy --manifest-path src-tauri/Cargo.toml --all-targets -- -D warnings
```

## Packaging

```bash
# Windows x64
npm run tauri build -- --bundles nsis

# macOS Apple Silicon
npm run tauri build -- --target aarch64-apple-darwin --bundles dmg
```

Tagged releases are built by `.github/workflows/release.yml`. Branch and pull-request artifacts
are produced by `.github/workflows/desktop-ci.yml`.

## Signing credentials

The current release uses an ad-hoc macOS signature and no Windows certificate. When credentials
are available, add them as GitHub Actions secrets and configure Tauri signing rather than
committing certificate files:

- Apple: Developer ID Application certificate, certificate password, signing identity, Apple ID
  or App Store Connect API key, team ID, and notarization settings.
- Windows: Authenticode certificate or a cloud signing command supported by Tauri.

## Project layout

- `src/`: React pages, workflow state, result tables, and charts.
- `src-tauri/src/lib.rs`: Tauri commands, parameter validation, job control, and result adapters.
- `src-tauri/capabilities/default.json`: minimal dialog, opener, and store permissions.
