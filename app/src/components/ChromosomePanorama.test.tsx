/**
 * Chromosome panorama interaction and helper tests.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import {
  ChromosomePanorama,
  chromosomePanoramaTestHelpers,
} from "./ChromosomePanorama";
import type { GenomeLocus } from "../types";

const loci: GenomeLocus[] = [
  {
    name: "chr1_nlr1",
    seq: "chr1",
    start: 10_000_000,
    end: 11_000_000,
    strand: "+",
    class: "CC-NBARC-LRR",
    family: "CNL",
    has_stop: false,
    motif_count: 5,
  },
  {
    name: "chr1_nlr2",
    seq: "chr1",
    start: 70_000_000,
    end: 71_000_000,
    strand: "-",
    class: "TIR-NBARC-LRR",
    family: "TNL",
    has_stop: false,
    motif_count: 6,
  },
];

describe("ChromosomePanorama", () => {
  it("opens a locus detail dialog and displays physical coordinates", () => {
    render(
      <ChromosomePanorama
        chromosomes={{
          chr1: { motifs: 8, nlrs: 2, length: 100_000_000 },
        }}
        loci={loci}
      />,
    );

    fireEvent.click(
      screen.getByRole("button", { name: /chr1_nlr2.*TNL/i }),
    );
    expect(screen.getByRole("dialog")).toBeInTheDocument();
    expect(screen.getByText("TIR-NBARC-LRR")).toBeInTheDocument();
    expect(screen.getByText("70.0–71.0 Mb")).toBeInTheDocument();

    fireEvent.click(
      screen.getByRole("button", { name: /close locus details/i }),
    );
    expect(screen.queryByRole("dialog")).not.toBeInTheDocument();
  });

  it("uses stable nice axis steps", () => {
    expect(chromosomePanoramaTestHelpers.niceStep(37)).toBe(50);
    expect(chromosomePanoramaTestHelpers.niceStep(4.2)).toBe(5);
    expect(chromosomePanoramaTestHelpers.niceStep(0.12)).toBe(0.2);
  });
});
