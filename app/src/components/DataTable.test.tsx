/**
 * Result table behavior tests.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import type { ColumnDef } from "@tanstack/react-table";
import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { DataTable } from "./DataTable";

interface Row {
  name: string;
  family: string;
}

const columns: ColumnDef<Row>[] = [
  { accessorKey: "name", header: "Name" },
  { accessorKey: "family", header: "Family" },
];

describe("DataTable", () => {
  it("filters result rows", () => {
    render(
      <DataTable
        data={[
          { name: "locus_1", family: "CNL" },
          { name: "locus_2", family: "TNL" },
        ]}
        columns={columns}
      />,
    );
    expect(screen.getByText("locus_1")).toBeInTheDocument();
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "TNL" } });
    expect(screen.queryByText("locus_1")).not.toBeInTheDocument();
    expect(screen.getByText("locus_2")).toBeInTheDocument();
  });
});
