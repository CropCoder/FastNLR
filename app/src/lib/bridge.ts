/**
 * Typed Tauri bridge and browser-preview fallbacks.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import { openPath } from "@tauri-apps/plugin-opener";
import { getCurrentWebviewWindow } from "@tauri-apps/api/webviewWindow";
import type {
  InputInfo,
  InputMode,
  JobSnapshot,
  RunDraft,
  RunResults,
  RuntimeInfo,
} from "../types";

function isTauriRuntime() {
  return typeof window !== "undefined" && "__TAURI_INTERNALS__" in window;
}

export async function getRuntimeInfo(): Promise<RuntimeInfo> {
  if (!isTauriRuntime()) {
    return {
      version: "1.9.2",
      os: "browser",
      arch: "web",
      logicalCpus: navigator.hardwareConcurrency || 1,
    };
  }
  return invoke<RuntimeInfo>("get_runtime_info");
}

export async function inspectInput(
  path: string,
  _mode: InputMode = "genome",
): Promise<InputInfo> {
  if (!isTauriRuntime()) {
    const name = path.split(/[\\/]/).pop() || "input.fasta";
    const prefix = name.replace(/\.gz$/i, "").replace(/\.(fa|fasta|fna|fas|ffn)$/i, "");
    return { path, name, prefix, sizeBytes: 0, compressed: /\.gz$/i.test(name) };
  }
  return invoke<InputInfo>("inspect_input", { path, mode: _mode });
}

export async function pickGenome(): Promise<string | null> {
  if (!isTauriRuntime()) return null;
  const result = await open({
    multiple: false,
    directory: false,
    title: "Select genome FASTA",
    filters: [
      {
        name: "Genome FASTA",
        extensions: ["fa", "fasta", "fna", "fas", "ffn", "gz"],
      },
    ],
  });
  return typeof result === "string" ? result : null;
}

export async function pickGff(): Promise<string | null> {
  if (!isTauriRuntime()) return null;
  const result = await open({
    multiple: false,
    directory: false,
    title: "Select GFF3 annotation",
    filters: [{ name: "GFF3", extensions: ["gff", "gff3", "gz"] }],
  });
  return typeof result === "string" ? result : null;
}

export async function pickFiles(title: string): Promise<string[]> {
  if (!isTauriRuntime()) return [];
  const result = await open({
    multiple: true,
    directory: false,
    title,
    filters: [{ name: "Model files", extensions: ["txt", "hmm", "dat"] }],
  });
  if (Array.isArray(result)) return result;
  return typeof result === "string" ? [result] : [];
}

export async function pickDirectory(title: string): Promise<string | null> {
  if (!isTauriRuntime()) return null;
  const result = await open({ directory: true, multiple: false, title });
  return typeof result === "string" ? result : null;
}

export async function startRun(draft: RunDraft): Promise<number> {
  if (!isTauriRuntime()) throw new Error("FastNLR jobs can only run in the desktop app.");
  return invoke<number>(
    draft.mode === "genome" ? "start_genome_job" : "start_annotation_job",
    { request: toBackendRequest(draft) },
  );
}

export function toBackendRequest(draft: RunDraft) {
  if (draft.mode === "genome") {
    return {
      inputPath: draft.inputPath,
      outputDir: draft.outputDir,
      threads: draft.threads,
      flank: draft.flank,
      seqsPerThread: draft.seqsPerThread,
      checkpointDir: draft.checkpointDir || null,
      motifAcceptP: Number(draft.motifAcceptP),
      motifPrelimP: Number(draft.motifPrelimP),
      tirAcceptP: Number(draft.tirAcceptP),
      flexibleSeed: draft.flexibleSeed ? Number(draft.flexibleSeed) : null,
      noFlexibleSeed: draft.noFlexibleSeed,
      requirePloop: draft.requirePloop,
      relaxedSeed: draft.relaxedSeed ? Number(draft.relaxedSeed) : null,
      keepPartialLoci: draft.keepPartialLoci,
      motFile: draft.motFile || null,
      storeFile: draft.storeFile || null,
      motifCategories: draft.motifCategories ? [draft.motifCategories] : [],
      extraSeeds: draft.extraSeeds ? [draft.extraSeeds] : [],
      extraSignatures: draft.extraSignatures ? [draft.extraSignatures] : [],
    };
  }
  return {
    inputPath: draft.inputPath,
    gffPath: draft.gffPath,
    outputDir: draft.outputDir,
    threads: draft.threads,
    batch: draft.batch,
    hmmFiles: draft.hmmFiles,
  };
}

export async function getJob(jobId: number): Promise<JobSnapshot> {
  return invoke<JobSnapshot>("get_job", { jobId });
}

export async function cancelJob(jobId: number): Promise<boolean> {
  return invoke<boolean>("cancel_job", { jobId });
}

export async function hasExistingResults(draft: RunDraft): Promise<boolean> {
  if (!isTauriRuntime()) return false;
  const input = await inspectInput(draft.inputPath, draft.mode);
  return invoke<boolean>("has_existing_results", {
    outputDir: draft.outputDir,
    prefix: input.prefix,
    mode: draft.mode,
  });
}

export async function loadResults(entry: {
  mode: InputMode;
  outputDir: string;
  prefix: string;
}): Promise<RunResults> {
  return invoke<RunResults>("load_results", { request: entry });
}

export async function revealPath(path: string): Promise<void> {
  if (!isTauriRuntime()) return;
  await openPath(path);
}

export async function listenForFileDrop(
  handler: (path: string) => void,
): Promise<() => void> {
  if (!isTauriRuntime()) return () => undefined;
  return getCurrentWebviewWindow().onDragDropEvent((event) => {
    if (event.payload.type === "drop" && event.payload.paths[0]) {
      handler(event.payload.paths[0]);
    }
  });
}
