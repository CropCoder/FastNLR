// 检查桌面端、CLI 与 Release 标签版本是否一致。
// Author: Jiwen Zhao (https://github.com/CropCoder)
// 用法：node scripts/S1-Check-Versions-V1.mjs [vX.Y.Z]

import fs from "node:fs";
import path from "node:path";
import process from "node:process";
import { fileURLToPath } from "node:url";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(scriptDir, "../..");

function read(relativePath) {
  return fs.readFileSync(path.join(repoRoot, relativePath), "utf8");
}

function cargoVersion(relativePath) {
  const text = read(relativePath);
  const match = text.match(/^version\s*=\s*"([^"]+)"/m);
  if (!match) throw new Error(`Cannot find package version in ${relativePath}`);
  return match[1];
}

const jsonVersions = {
  "app/package.json": JSON.parse(read("app/package.json")).version,
  "app/package-lock.json": JSON.parse(read("app/package-lock.json")).version,
  "app/src-tauri/tauri.conf.json": JSON.parse(
    read("app/src-tauri/tauri.conf.json"),
  ).version,
};

const cargoVersions = {
  "Cargo.toml": cargoVersion("Cargo.toml"),
  "app/src-tauri/Cargo.toml": cargoVersion("app/src-tauri/Cargo.toml"),
};

const versions = { ...jsonVersions, ...cargoVersions };
const expected = process.argv[2]?.replace(/^v/, "");
const unique = new Set(Object.values(versions));

if (unique.size !== 1 || (expected && !unique.has(expected))) {
  console.error("Version mismatch detected:");
  for (const [file, version] of Object.entries(versions)) {
    console.error(`  ${file}: ${version}`);
  }
  if (expected) console.error(`  release tag: ${expected}`);
  process.exit(1);
}

console.log(`Version check passed: ${expected ?? [...unique][0]}`);
