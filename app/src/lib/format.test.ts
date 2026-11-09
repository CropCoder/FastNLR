/**
 * Unit tests for display formatting helpers.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { describe, expect, it } from "vitest";
import { countEntries, formatBytes, toSeries } from "./format";

describe("format helpers", () => {
  it("formats byte sizes", () => {
    expect(formatBytes(512)).toBe("512 B");
    expect(formatBytes(1536)).toBe("1.50 KB");
    expect(formatBytes(20 * 1024 * 1024)).toBe("20.0 MB");
  });

  it("sorts chart series by count and name", () => {
    expect(
      toSeries({ TNL: 4, CNL: 8, RNL: 4 }),
    ).toEqual([
      { name: "CNL", value: 8 },
      { name: "RNL", value: 4 },
      { name: "TNL", value: 4 },
    ]);
    expect(countEntries({ TNL: 4, CNL: 8 })).toBe(12);
  });
});
