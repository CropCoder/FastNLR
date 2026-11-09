/**
 * Annotation form behavior tests.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { fireEvent, render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { beforeEach, describe, expect, it, vi } from "vitest";
import { NewRunPage, newRunTestHelpers } from "./NewRunPage";
import { useAppStore } from "../store/appStore";

vi.mock("../lib/bridge", () => ({
  inspectInput: vi.fn(),
  listenForFileDrop: vi.fn().mockResolvedValue(() => undefined),
  pickDirectory: vi.fn(),
  pickFasta: vi.fn(),
  pickFiles: vi.fn(),
  hasExistingResults: vi.fn(),
  startRun: vi.fn(),
}));

describe("NewRunPage", () => {
  beforeEach(() => {
    useAppStore.setState({
      initialized: true,
      runtime: {
        version: "1.8.0",
        os: "linux",
        arch: "x86_64",
        logicalCpus: 8,
      },
      settings: { defaultThreads: 8, defaultOutputDir: "" },
      history: [],
      activeRun: null,
      pendingDraft: null,
    });
  });

  it("switches to protein mode and reveals HMM controls", () => {
    render(
      <MemoryRouter>
        <NewRunPage />
      </MemoryRouter>,
    );

    expect(screen.getByText("Launch a new NLR analysis")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("tab", { name: /Protein/i }));
    fireEvent.click(screen.getByRole("button", { name: /Advanced parameters/i }));

    expect(screen.getByText("Custom HMM files")).toBeInTheDocument();
    expect(screen.queryByText("Keep partial loci")).not.toBeInTheDocument();
  });

  it("maps genome form values to the backend draft contract", () => {
    const values = {
      ...newRunTestHelpers.baseDefaults(),
      inputPath: "/data/genome.fa",
      outputDir: "/data/results",
      threads: 8,
      motifCategories: "21=CC",
      extraSeeds: "21,4",
    };
    const draft = newRunTestHelpers.formToDraft(values);
    expect(draft.mode).toBe("genome");
    if (draft.mode === "genome") {
      expect(draft.inputPath).toBe("/data/genome.fa");
      expect(draft.threads).toBe(8);
      expect(draft.extraSeeds).toBe("21,4");
    }
  });
});
