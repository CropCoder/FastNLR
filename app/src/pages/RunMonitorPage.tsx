/**
 * Active annotation progress and cancellation page.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Activity,
  Check,
  CircleStop,
  Clock3,
  Database,
  FileOutput,
  LoaderCircle,
  TriangleAlert,
} from "lucide-react";
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { PageHeader } from "../components/PageHeader";
import { StatusBadge } from "../components/StatusBadge";
import { cancelJob } from "../lib/bridge";
import { formatElapsed } from "../lib/format";
import { useAppStore } from "../store/appStore";

const phaseOrder = ["Pending", "ParsingAnnotation", "Scanning", "Writing", "Completed"] as const;

const phaseCopy: Record<string, string> = {
  Pending: "Waiting for the execution thread",
  ParsingAnnotation: "Parsing GFF3 and selecting the longest CDS transcript per gene",
  Scanning: "Scanning sequence fragments and assembling loci",
  Writing: "Writing output tables, sequences, and summaries",
  Completed: "All output files are ready",
  Failed: "The annotation pipeline failed",
  Canceled: "The run was canceled",
};

export function RunMonitorPage() {
  const navigate = useNavigate();
  const activeRun = useAppStore((state) => state.activeRun);
  const history = useAppStore((state) => state.history);
  const setActiveRun = useAppStore((state) => state.setActiveRun);
  const [clock, setClock] = useState(0);
  const [cancelRequested, setCancelRequested] = useState(false);

  const entry = activeRun
    ? history.find((item) => item.id === activeRun.historyId)
    : null;

  useEffect(() => {
    const timer = window.setInterval(() => setClock((value) => value + 1), 1000);
    return () => window.clearInterval(timer);
  }, []);

  useEffect(() => {
    if (!activeRun || entry?.status !== "Completed") return;
    setActiveRun(null);
    navigate(`/results/${activeRun.historyId}`, { replace: true });
  }, [activeRun, entry?.status, navigate, setActiveRun]);

  if (!activeRun || !entry) {
    return (
      <div className="page">
        <PageHeader
          eyebrow="Run monitor"
          title="No active annotation run"
          description="Start a new analysis or open a completed run from History."
        />
        <section className="panel centered-panel">
          <Activity size={40} />
          <h2>Nothing is running</h2>
          <button className="btn primary" onClick={() => navigate("/new")}>
            Start a new run
          </button>
        </section>
      </div>
    );
  }

  const phase = entry.phase;
  const status = entry.status;
  const scanned = entry.scanned;
  const total = entry.estimatedTotal;
  const progress = total > 0 ? Math.min(100, (scanned / total) * 100) : null;
  const terminal = status === "Failed" || status === "Canceled";
  const phaseIndex = Math.max(0, phaseOrder.indexOf(phase as (typeof phaseOrder)[number]));

  async function requestCancel() {
    if (!activeRun || cancelRequested) return;
    const confirmed = window.confirm(
      "Cancel this annotation run? Partial scan results will not be written.",
    );
    if (!confirmed) return;
    setCancelRequested(await cancelJob(activeRun.jobId));
  }

  return (
    <div className="page">
      <PageHeader
        eyebrow="Run monitor"
        title={entry.inputName}
        description={`${activeRun.mode === "genome" ? "Genome annotation" : "GFF3 protein annotation"} writing to ${activeRun.outputDir}`}
        actions={<StatusBadge status={status} phase={phase} />}
      />

      <section className="monitor-hero panel">
        <div className="monitor-orbit" data-state={status.toLowerCase()}>
          {status === "Completed" ? (
            <Check size={38} />
          ) : terminal ? (
            <TriangleAlert size={38} />
          ) : (
            <LoaderCircle className="spin" size={42} />
          )}
        </div>
        <div className="monitor-copy">
          <span className="eyebrow">Current phase</span>
          <h2>{phase}</h2>
          <p>{phaseCopy[phase] ?? "Processing the annotation pipeline"}</p>
        </div>
        <div className="monitor-time">
          <Clock3 size={18} />
          <span>Elapsed</span>
          <strong>{formatElapsed(activeRun.startedAt)}</strong>
        </div>
      </section>

      <section className="panel progress-panel">
        <div className="progress-heading">
          <div>
            <h2>Pipeline progress</h2>
            <p>
              {scanned.toLocaleString()} {entry.mode === "genome" ? "fragments" : "genes"} processed
              {total > 0 ? ` of approximately ${total.toLocaleString()}` : ""}
            </p>
          </div>
          <strong>{progress == null ? "—" : `${progress.toFixed(1)}%`}</strong>
        </div>
        <div className={`progress-track${progress == null ? " indeterminate" : ""}`}>
          <span style={{ width: progress == null ? undefined : `${progress}%` }} />
        </div>

        <div className="phase-rail">
          {phaseOrder.map((item, index) => (
            <div
              key={item}
              className={`phase-step${index < phaseIndex ? " done" : ""}${item === phase ? " current" : ""}`}
            >
              <span>{index < phaseIndex ? <Check size={14} /> : index + 1}</span>
              <div>
                <strong>{item}</strong>
                <small>{phaseCopy[item]}</small>
              </div>
            </div>
          ))}
        </div>
      </section>

      <div className="monitor-meta-grid">
        <article className="panel compact-panel">
          <Database size={20} />
          <span>Input</span>
          <strong title={entry.inputPath}>{entry.inputName}</strong>
        </article>
        <article className="panel compact-panel">
          <FileOutput size={20} />
          <span>Output prefix</span>
          <strong>{entry.prefix}</strong>
        </article>
        <article className="panel compact-panel">
          <Activity size={20} />
          <span>Phase</span>
          <strong>{phase}</strong>
        </article>
      </div>

      {entry.error ? <div className="alert error">{entry.error}</div> : null}

      <div className="monitor-actions">
        <button className="btn ghost" onClick={() => navigate("/new")}>
          New analysis
        </button>
        {!terminal && status !== "Completed" ? (
          <button
            className="btn danger"
            onClick={requestCancel}
            disabled={cancelRequested}
          >
            <CircleStop size={18} />
            {cancelRequested ? "Cancellation requested" : "Cancel run"}
          </button>
        ) : status === "Completed" ? (
          <button
            className="btn primary"
            onClick={() => navigate(`/results/${activeRun.historyId}`)}
          >
            View results
          </button>
        ) : (
          <button
            className="btn secondary"
            onClick={() => {
              setActiveRun(null);
              navigate("/history");
            }}
          >
            Open history
          </button>
        )}
      </div>
    </div>
  );
}
