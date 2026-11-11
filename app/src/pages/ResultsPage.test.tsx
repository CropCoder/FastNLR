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
    per_chromosome: { chr1: { motifs: 9, nlrs: 2, length: 250_000_000 } },
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
  annotation: null,
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
    expect(
      screen.getByRole("button", { name: /Chromosome map/i }),
    ).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /NLR loci/i }));
    expect(await screen.findByText("nlr_1")).toBeInTheDocument();
    expect(screen.getAllByText("chr1").length).toBeGreaterThan(0);
  });

  it("renders GFF3 annotation gene and domain tables", async () => {
    const annotationEntry: RunHistoryEntry = {
      ...entry,
      id: "run-gff",
      mode: "annotation",
      inputName: "genome.fa",
      draft: {
        mode: "annotation",
        inputPath: "/data/genome.fa",
        gffPath: "/data/annotation.gff3",
        outputDir: "/data/results",
        threads: 4,
        batch: 1000,
        hmmFiles: [],
      },
    };
    useAppStore.setState({ history: [annotationEntry] });
    loadResults.mockResolvedValue({
      mode: "annotation",
      prefix: "genome",
      outputDir: "/data/results",
      outputs: [],
      genome: null,
      annotation: {
        proteinsScanned: 1,
        nlrTotal: 1,
        completeNlrTotal: 1,
        genesWithCds: 1,
        skippedGenes: 0,
        elapsed: "1.00s",
        modelLabels: ["PF00931(NB-ARC)"],
        modelLibraryVersion: "FastNLR-Pfam-HMM-23-2026-11-11",
        modelCount: 23,
        typeCounts: { "CC-NBARC-LRR": 1 },
        completeTypeCounts: { "CC-NBARC-LRR": 1 },
        familyCounts: { CNL: 1 },
        records: [
          {
            geneId: "g1",
            geneName: "RPM1",
            geneBiotype: "protein_coding",
            transcriptId: "t1",
            transcriptName: "",
            seqid: "chr1",
            geneStart: 1,
            geneEnd: 100,
            strand: "+",
            cdsCount: 1,
            cdsLengthBp: 90,
            proteinLength: 30,
            class: "CC-NBARC-LRR",
            nterm: "CC",
            hasNbarc: true,
            hasTir: false,
            hasLrrAny: true,
            hasCtermLrr: true,
            complete: true,
            isNlr: true,
            domainCount: 3,
            domains: "PF00931:5-20",
            family: "CNL",
            internalStopCount: 0,
            hasInternalStop: false,
            trailingPartialBases: 0,
            phaseConsistent: true,
          },
        ],
        domains: [
          {
            geneId: "g1",
            transcriptId: "t1",
            seqid: "chr1",
            strand: "+",
            source: "hmm",
            model: "PF00931",
            modelName: "NB-ARC",
            category: "NBARC",
            proteinStart: 5,
            proteinEnd: 20,
            cdsStart: 13,
            cdsEnd: 60,
            genomicStart: 13,
            genomicEnd: 60,
            genomicSegments: "chr1:13-60",
            bitscore: 100,
            pvalue: 1e-20,
            passed: true,
          },
        ],
      },
    });

    render(
      <MemoryRouter initialEntries={["/results/run-gff"]}>
        <Routes>
          <Route path="/results/:historyId" element={<ResultsPage />} />
        </Routes>
      </MemoryRouter>,
    );

    await waitFor(() => expect(loadResults).toHaveBeenCalledWith({
      mode: "annotation",
      outputDir: "/data/results",
      prefix: "genome",
    }));
    expect(screen.getByText("Complete NLR")).toBeInTheDocument();
    expect(screen.getByText(/23 models/i)).toBeInTheDocument();
    fireEvent.click(await screen.findByRole("button", { name: /NLR genes/i }));
    expect(await screen.findByText("g1")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /Domain hits/i }));
    expect(await screen.findByText("PF00931")).toBeInTheDocument();
  });
});
