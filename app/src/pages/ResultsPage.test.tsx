/**
 * Desktop result loading and table tests.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, Route, Routes } from "react-router-dom";
import { beforeEach, describe, expect, it, vi } from "vitest";
import { ResultsPage } from "./ResultsPage";
import { useAppStore } from "../store/appStore";
import type { RunHistoryEntry, RunResults } from "../types";

const loadResults = vi.fn();

vi.mock("../lib/bridge", () => ({
  loadResults: (...args: unknown[]) => loadResults(...args),
  revealPath: vi.fn(),
}));

const entry: RunHistoryEntry = {
  id: "run-1",
  mode: "genome",
  inputPath: "/data/genome.fa",
  inputName: "genome.fa",
  outputDir: "/data/results",
  prefix: "genome",
  status: "Completed",
  phase: "Completed",
  startedAt: "2026-11-09T08:00:00Z",
  finishedAt: "2026-11-09T08:01:00Z",
  scanned: 10,
  estimatedTotal: 10,
  resultCount: 2,
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
};

const results: RunResults = {
  mode: "genome",
  prefix: "genome",
  outputDir: "/data/results",
  outputs: [
    {
      name: "genome.summary.json",
      path: "/data/results/genome.summary.json",
      sizeBytes: 1024,
      isDir: false,
    },
  ],
  genome: {
    sequence_count: 1,
    motif_total: 9,
    nlr_total: 2,
    per_chromosome: { chr1: { motifs: 9, nlrs: 2 } },
    motif_counts: { "1": 2, "6": 2, "4": 2 },
    nlr_type_counts: { "CC-NBARC-LRR": 2 },
    nlr_family_counts: { CNL: 2 },
    loci: [
      {
        name: "nlr_1",
        seq: "chr1",
        start: 100,
        end: 900,
        strand: "+",
        class: "CC-NBARC-LRR",
        family: "CNL",
        has_stop: false,
        motif_count: 5,
      },
      {
        name: "nlr_2",
        seq: "chr1",
        start: 1200,
        end: 2100,
        strand: "-",
        class: "CC-NBARC-LRR",
        family: "CNL",
        has_stop: false,
        motif_count: 4,
      },
    ],
  },
  protein: null,
};

describe("ResultsPage", () => {
  beforeEach(() => {
    loadResults.mockResolvedValue(results);
    useAppStore.setState({
      initialized: true,
      history: [entry],
      activeRun: null,
    });
  });

  it("loads result metrics and switches to the locus table", async () => {
    render(
      <MemoryRouter initialEntries={["/results/run-1"]}>
        <Routes>
          <Route path="/results/:historyId" element={<ResultsPage />} />
        </Routes>
      </MemoryRouter>,
    );

    await waitFor(() => expect(loadResults).toHaveBeenCalled());
    expect(screen.getByText(/2 NLR candidates/i)).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /NLR loci/i }));
    expect(await screen.findByText("nlr_1")).toBeInTheDocument();
    expect(screen.getAllByText("chr1").length).toBeGreaterThan(0);
  });
});
