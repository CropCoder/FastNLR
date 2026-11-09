/**
 * Genome and protein result overview, tables, and output browser.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import type { ColumnDef } from "@tanstack/react-table";
import {
  AlertCircle,
  Boxes,
  Dna,
  FileJson,
  Files,
  FolderOpen,
  Layers3,
  RefreshCw,
  RotateCcw,
  Table2,
  type LucideIcon,
} from "lucide-react";
import { useEffect, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { DataTable } from "../components/DataTable";
import { DistributionChart } from "../components/DistributionChart";
import { EmptyState } from "../components/EmptyState";
import { MetricCard } from "../components/MetricCard";
import { PageHeader } from "../components/PageHeader";
import { loadResults, revealPath } from "../lib/bridge";
import { formatBytes, formatDate } from "../lib/format";
import { useAppStore } from "../store/appStore";
import type {
  GenomeLocus,
  ProteinDomain,
  ProteinRecord,
  RunResults,
} from "../types";

const locusColumns: ColumnDef<GenomeLocus>[] = [
  { accessorKey: "name", header: "Locus" },
  { accessorKey: "seq", header: "Sequence" },
  {
    accessorKey: "start",
    header: "Start",
    cell: ({ getValue }) => Number(getValue()).toLocaleString(),
  },
  {
    accessorKey: "end",
    header: "End",
    cell: ({ getValue }) => Number(getValue()).toLocaleString(),
  },
  { accessorKey: "strand", header: "Strand" },
  {
    accessorKey: "class",
    header: "Architecture",
    cell: ({ getValue }) => <span className="table-code">{String(getValue())}</span>,
  },
  { accessorKey: "family", header: "Family" },
  { accessorKey: "motif_count", header: "Motifs" },
  {
    accessorKey: "has_stop",
    header: "Stop",
    cell: ({ getValue }) => (getValue() ? "Yes" : "No"),
  },
];

const proteinColumns: ColumnDef<ProteinRecord>[] = [
  { accessorKey: "proteinId", header: "Protein" },
  {
    accessorKey: "length",
    header: "Length",
    cell: ({ getValue }) => Number(getValue()).toLocaleString(),
  },
  {
    accessorKey: "class",
    header: "Architecture",
    cell: ({ getValue }) => <span className="table-code">{String(getValue())}</span>,
  },
  { accessorKey: "family", header: "Family" },
  { accessorKey: "nterm", header: "N-terminus" },
  { accessorKey: "domainCount", header: "Domains" },
  {
    accessorKey: "domains",
    header: "Domain layout",
    cell: ({ getValue }) => (
      <span className="domain-layout" title={String(getValue())}>
        {String(getValue())}
      </span>
    ),
  },
];

const domainColumns: ColumnDef<ProteinDomain>[] = [
  { accessorKey: "proteinId", header: "Protein" },
  { accessorKey: "source", header: "Source" },
  {
    accessorKey: "model",
    header: "Model",
    cell: ({ row }) => (
      <span>
        {row.original.model}
        <small className="cell-subtitle">{row.original.modelName}</small>
      </span>
    ),
  },
  { accessorKey: "category", header: "Category" },
  { accessorKey: "start", header: "Start" },
  { accessorKey: "end", header: "End" },
  {
    accessorKey: "bitscore",
    header: "Bitscore",
    cell: ({ getValue }) =>
      getValue() == null ? "—" : Number(getValue()).toFixed(1),
  },
  {
    accessorKey: "pvalue",
    header: "p-value",
    cell: ({ getValue }) =>
      getValue() == null ? "—" : Number(getValue()).toExponential(2),
  },
  {
    accessorKey: "passed",
    header: "Passed",
    cell: ({ getValue }) => (
      <span className={`pass-mark ${getValue() ? "yes" : "no"}`}>
        {getValue() ? "Yes" : "No"}
      </span>
    ),
  },
];

type ResultTab = "overview" | "loci" | "domains" | "files";

const resultTabs = (mode: "genome" | "protein"): Array<{
  value: ResultTab;
  label: string;
  icon: LucideIcon;
}> => [
  { value: "overview", label: "Overview", icon: Layers3 },
  {
    value: "loci",
    label: mode === "genome" ? "NLR loci" : "NLR proteins",
    icon: Table2,
  },
  ...(mode === "protein"
    ? [{ value: "domains" as const, label: "Domain hits", icon: Boxes }]
    : []),
  { value: "files", label: "Output files", icon: Files },
];

export function ResultsPage() {
  const navigate = useNavigate();
  const { historyId } = useParams();
  const history = useAppStore((state) => state.history);
  const activeRun = useAppStore((state) => state.activeRun);
  const updateHistory = useAppStore((state) => state.updateHistory);
  const setPendingDraft = useAppStore((state) => state.setPendingDraft);
  const setActiveRun = useAppStore((state) => state.setActiveRun);
  const entry = history.find((item) => item.id === historyId);
  const [results, setResults] = useState<RunResults | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [tab, setTab] = useState<ResultTab>("overview");
  const [reloadKey, setReloadKey] = useState(0);

  useEffect(() => {
    if (activeRun?.historyId === historyId) setActiveRun(null);
  }, [activeRun?.historyId, historyId, setActiveRun]);

  useEffect(() => {
    if (!entry) return;
    let disposed = false;
    setLoading(true);
    setError(null);
    void loadResults({
      mode: entry.mode,
      outputDir: entry.outputDir,
      prefix: entry.prefix,
    })
      .then((value) => {
        if (disposed) return;
        setResults(value);
        const count = value.genome?.nlr_total ?? value.protein?.nlrTotal ?? 0;
        updateHistory(entry.id, { resultCount: count });
      })
      .catch((reason) => {
        if (!disposed) setError(String(reason));
      })
      .finally(() => {
        if (!disposed) setLoading(false);
      });
    return () => {
      disposed = true;
    };
  }, [entry?.id, reloadKey]);

  if (!entry) {
    return (
      <div className="page">
        <PageHeader
          eyebrow="Results"
          title="Run record not found"
          description="The local history record may have been removed."
        />
        <section className="panel">
          <EmptyState
            icon={<AlertCircle size={30} />}
            title="No matching run"
            description="Choose a completed run from History."
            action={
              <button className="btn primary" onClick={() => navigate("/history")}>
                Open History
              </button>
            }
          />
        </section>
      </div>
    );
  }

  if (loading) {
    return (
      <div className="page">
        <PageHeader
          eyebrow="Results"
          title={entry.inputName}
          description="Loading the machine-readable result set..."
        />
        <div className="result-skeleton">
          <div />
          <div />
          <div />
        </div>
      </div>
    );
  }

  if (error || !results) {
    return (
      <div className="page">
        <PageHeader
          eyebrow="Results"
          title={entry.inputName}
          description={`Output directory: ${entry.outputDir}`}
        />
        <section className="panel">
          <EmptyState
            icon={<AlertCircle size={30} />}
            title="Results could not be loaded"
            description={error || "The generated summary files are missing or incomplete."}
            action={
              <div className="inline-actions">
                <button className="btn secondary" onClick={() => setReloadKey((value) => value + 1)}>
                  <RefreshCw size={16} />
                  Retry
                </button>
                <button className="btn ghost" onClick={() => void revealPath(entry.outputDir)}>
                  <FolderOpen size={16} />
                  Open folder
                </button>
              </div>
            }
          />
        </section>
      </div>
    );
  }

  const genome = results.genome;
  const protein = results.protein;
  const resultCount = genome?.nlr_total ?? protein?.nlrTotal ?? 0;
  const chromosomeNlrs = genome
    ? Object.fromEntries(
        Object.entries(genome.per_chromosome)
          .map(([name, value]) => [name, value.nlrs] as const)
          .sort((left, right) => right[1] - left[1])
          .slice(0, 20),
      )
    : {};

  return (
    <div className="page results-page">
      <PageHeader
        eyebrow={`${entry.mode === "genome" ? "Genome" : "Protein"} results`}
        title={entry.inputName}
        description={`${resultCount.toLocaleString()} NLR candidates · completed ${formatDate(entry.finishedAt || entry.startedAt)}`}
        actions={
          <div className="inline-actions">
            <button className="btn secondary" onClick={() => void revealPath(entry.outputDir)}>
              <FolderOpen size={16} />
              Open output folder
            </button>
            <button
              className="btn ghost"
              onClick={() => {
                setPendingDraft(entry.draft);
                navigate("/new");
              }}
            >
              <RotateCcw size={16} />
              Run again
            </button>
          </div>
        }
      />

      <nav className="tab-bar" aria-label="Result views">
        {resultTabs(entry.mode).map(({ value, label, icon: Icon }) => (
          <button
            key={value}
            type="button"
            className={tab === value ? "active" : ""}
            onClick={() => setTab(value as ResultTab)}
          >
            <Icon size={16} />
            {label}
          </button>
        ))}
      </nav>

      {tab === "overview" ? (
        <>
          <div className="metrics-grid">
            <MetricCard
              label="NLR candidates"
              value={resultCount.toLocaleString()}
              detail={entry.mode === "genome" ? "assembled loci" : "classified proteins"}
              icon={<Dna size={19} />}
              accent="mint"
            />
            {genome ? (
              <>
                <MetricCard
                  label="Motif hits"
                  value={genome.motif_total.toLocaleString()}
                  detail="accepted sequence motifs"
                  icon={<Layers3 size={19} />}
                  accent="iris"
                />
                <MetricCard
                  label="Sequences"
                  value={genome.sequence_count.toLocaleString()}
                  detail="chromosomes and contigs"
                  icon={<Table2 size={19} />}
                  accent="sun"
                />
                <MetricCard
                  label="Families"
                  value={Object.keys(genome.nlr_family_counts).length}
                  detail="NLR family classes"
                  icon={<Boxes size={19} />}
                  accent="coral"
                />
              </>
            ) : (
              <>
                <MetricCard
                  label="Proteins scanned"
                  value={(protein?.proteinsScanned ?? 0).toLocaleString()}
                  detail="input sequences"
                  icon={<Layers3 size={19} />}
                  accent="iris"
                />
                <MetricCard
                  label="Domain hits"
                  value={(protein?.domains.length ?? 0).toLocaleString()}
                  detail="HMM and LRR evidence"
                  icon={<Boxes size={19} />}
                  accent="sun"
                />
                <MetricCard
                  label="Elapsed"
                  value={protein?.elapsed ?? "—"}
                  detail={`${protein?.modelLabels.length ?? 0} domain models`}
                  icon={<RefreshCw size={19} />}
                  accent="coral"
                />
              </>
            )}
          </div>

          <div className="charts-grid">
            <DistributionChart
              title="NLR families"
              subtitle="Biological family classification"
              values={
                genome?.nlr_family_counts ?? protein?.familyCounts ?? {}
              }
            />
            <DistributionChart
              title="Domain architectures"
              subtitle="Detected domain combinations"
              values={genome?.nlr_type_counts ?? protein?.typeCounts ?? {}}
            />
            {genome ? (
              <DistributionChart
                title="NLR loci by sequence"
                subtitle="Top chromosomes and contigs"
                values={chromosomeNlrs}
                horizontal
              />
            ) : (
              <article className="panel chart-card result-summary-card">
                <div className="chart-heading">
                  <div>
                    <h3>Domain model set</h3>
                    <p>Models used for protein classification</p>
                  </div>
                </div>
                <div className="model-tags">
                  {protein?.modelLabels.map((label) => (
                    <span key={label}>{label}</span>
                  ))}
                </div>
              </article>
            )}
          </div>
        </>
      ) : null}

      {tab === "loci" ? (
        <section className="panel table-panel">
          <div className="section-heading">
            <span className="section-number">01</span>
            <div>
              <h2>{entry.mode === "genome" ? "Assembled NLR loci" : "Classified NLR proteins"}</h2>
              <p>Sort and filter the complete result table.</p>
            </div>
          </div>
          {entry.mode === "genome" && genome ? (
            <DataTable
              data={genome.loci}
              columns={locusColumns}
              searchPlaceholder="Filter by locus, sequence, family, or architecture"
            />
          ) : protein ? (
            <DataTable
              data={protein.records}
              columns={proteinColumns}
              searchPlaceholder="Filter by protein, family, or architecture"
            />
          ) : null}
        </section>
      ) : null}

      {tab === "domains" && protein ? (
        <section className="panel table-panel">
          <div className="section-heading">
            <span className="section-number">02</span>
            <div>
              <h2>Domain evidence</h2>
              <p>All HMM and PWM LRR hits reported by the protein pipeline.</p>
            </div>
          </div>
          <DataTable
            data={protein.domains}
            columns={domainColumns}
            searchPlaceholder="Filter by protein, model, category, or source"
          />
        </section>
      ) : null}

      {tab === "files" ? (
        <section className="panel files-panel">
          <div className="section-heading">
            <span className="section-number">03</span>
            <div>
              <h2>Generated files</h2>
              <p>Open individual outputs or reveal the complete result directory.</p>
            </div>
          </div>
          <div className="output-file-list">
            {results.outputs.map((file) => (
              <button
                className="output-file"
                key={file.path}
                onClick={() => void revealPath(file.path)}
              >
                <span className="file-icon">
                  {file.name.endsWith(".json") ? <FileJson size={19} /> : <Files size={19} />}
                </span>
                <span>
                  <strong>{file.name}</strong>
                  <small title={file.path}>{file.path}</small>
                </span>
                <span>{file.isDir ? "Folder" : formatBytes(file.sizeBytes)}</span>
              </button>
            ))}
          </div>
        </section>
      ) : null}
    </div>
  );
}
