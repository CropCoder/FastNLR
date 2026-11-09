/**
 * Desktop state initialization and interruption tests.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { beforeEach, describe, expect, it, vi } from "vitest";
import { useAppStore } from "./appStore";
import type { RunHistoryEntry } from "../types";

const loadPersistedState = vi.fn();
const savePersistedState = vi.fn();
const getRuntimeInfo = vi.fn();

vi.mock("../lib/persistence", () => ({
  loadPersistedState: (...args: unknown[]) => loadPersistedState(...args),
  savePersistedState: (...args: unknown[]) => savePersistedState(...args),
}));

vi.mock("../lib/bridge", () => ({
  getRuntimeInfo: (...args: unknown[]) => getRuntimeInfo(...args),
}));

const interruptedEntry = {
  id: "run-1",
  mode: "genome",
  inputPath: "/data/genome.fa",
  inputName: "genome.fa",
  outputDir: "/data/results",
  prefix: "genome",
  status: "Running",
  phase: "Scanning",
  startedAt: "2026-11-09T08:00:00Z",
  finishedAt: null,
  scanned: 3,
  estimatedTotal: 10,
  resultCount: null,
  error: null,
  draft: {
    mode: "genome",
    inputPath: "/data/genome.fa",
    outputDir: "/data/results",
    threads: 4,
    flank: 2000,
    seqsPerThread: 1000,
    checkpointDir: "",
    motifAcceptP: "0.00001",
    motifPrelimP: "0.0001",
    tirAcceptP: "0.0000001",
    flexibleSeed: "3",
    noFlexibleSeed: false,
    requirePloop: false,
    relaxedSeed: "",
    keepPartialLoci: false,
    motFile: "",
    storeFile: "",
    motifCategories: "",
    extraSeeds: "",
    extraSignatures: "",
  },
} satisfies RunHistoryEntry;

describe("appStore initialization", () => {
  beforeEach(() => {
    getRuntimeInfo.mockResolvedValue({
      version: "1.8.0",
      os: "windows",
      arch: "x86_64",
      logicalCpus: 12,
    });
    loadPersistedState.mockResolvedValue({
      settings: { defaultThreads: 8, defaultOutputDir: "/data" },
      history: [interruptedEntry],
    });
    useAppStore.setState({
      initialized: false,
      runtime: null,
      settings: { defaultThreads: 1, defaultOutputDir: "" },
      history: [],
      activeRun: null,
      pendingDraft: null,
    });
  });

  it("marks unfinished persisted jobs as interrupted", async () => {
    await useAppStore.getState().initialize();
    const value = useAppStore.getState();
    expect(value.initialized).toBe(true);
    expect(value.runtime?.logicalCpus).toBe(12);
    expect(value.settings.defaultThreads).toBe(8);
    expect(value.history[0].status).toBe("Interrupted");
    expect(value.history[0].phase).toBe("Interrupted");
    expect(savePersistedState).toHaveBeenCalled();
  });
});
