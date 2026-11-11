/**
 * Local annotation run history.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Clock3,
  Dna,
  ExternalLink,
  FolderOpen,
  History,
  Play,
  RotateCcw,
  Trash2,
} from "lucide-react";
import { useNavigate } from "react-router-dom";
import { EmptyState } from "../components/EmptyState";
import { PageHeader } from "../components/PageHeader";
import { StatusBadge } from "../components/StatusBadge";
import { revealPath } from "../lib/bridge";
import { formatDate } from "../lib/format";
import { useAppStore } from "../store/appStore";

export function HistoryPage() {
  const navigate = useNavigate();
  const history = useAppStore((state) => state.history);
  const clearHistory = useAppStore((state) => state.clearHistory);
  const setPendingDraft = useAppStore((state) => state.setPendingDraft);
  const activeRun = useAppStore((state) => state.activeRun);

  function clearAll() {
    if (
      window.confirm(
        "Remove all local run records? Generated result folders and files will not be deleted.",
      )
    ) {
      clearHistory();
    }
  }

  return (
    <div className="page">
      <PageHeader
        eyebrow="Local run archive"
        title="Annotation history"
        description="Recent runs are stored locally. Output folders remain untouched when a history record is removed."
        actions={
          history.length ? (
            <button className="btn ghost danger-text" onClick={clearAll}>
              <Trash2 size={16} />
              Clear history
            </button>
          ) : null
        }
      />

      {history.length === 0 ? (
        <section className="panel">
          <EmptyState
            icon={<History size={32} />}
            title="No annotation runs yet"
            description="Start a genome or GFF3-guided analysis and its local summary will appear here."
            action={
              <button className="btn primary" onClick={() => navigate("/new")}>
                <Play size={16} />
                Start first run
              </button>
            }
          />
        </section>
      ) : (
        <div className="history-list">
          {history.map((entry) => {
            const isActive = activeRun?.historyId === entry.id;
            const canViewResults =
              entry.status === "Completed" || entry.status === "Interrupted";
            return (
              <article className="panel history-card" key={entry.id}>
                <div className={`history-mode ${entry.mode}`}>
                  <Dna size={21} />
                </div>
                <div className="history-main">
                  <div className="history-title-row">
                    <div>
                      <h3>{entry.inputName}</h3>
                      <p title={entry.inputPath}>{entry.inputPath}</p>
                    </div>
                    <StatusBadge status={entry.status} phase={entry.phase} />
                  </div>
                  <div className="history-meta">
                    <span>
                      <Clock3 size={14} />
                      {formatDate(entry.startedAt)}
                    </span>
                    <span>{entry.mode === "genome" ? "Genome mode" : "GFF3 mode"}</span>
                    <span>
                      {entry.resultCount == null
                        ? `${entry.scanned.toLocaleString()} processed`
                        : `${entry.resultCount.toLocaleString()} NLR candidates`}
                    </span>
                  </div>
                  {entry.error ? <p className="history-error">{entry.error}</p> : null}
                </div>
                <div className="history-actions">
                  {isActive ? (
                    <button className="btn primary small" onClick={() => navigate("/run")}>
                      <ExternalLink size={15} />
                      Monitor
                    </button>
                  ) : canViewResults ? (
                    <button
                      className="btn primary small"
                      onClick={() => navigate(`/results/${entry.id}`)}
                    >
                      <ExternalLink size={15} />
                      Results
                    </button>
                  ) : null}
                  <button
                    className="btn secondary small"
                    onClick={() => {
                      setPendingDraft(entry.draft);
                      navigate("/new");
                    }}
                  >
                    <RotateCcw size={15} />
                    Reuse
                  </button>
                  <button
                    className="icon-button"
                    aria-label="Open output directory"
                    onClick={() => void revealPath(entry.outputDir)}
                  >
                    <FolderOpen size={17} />
                  </button>
                </div>
              </article>
            );
          })}
        </div>
      )}
    </div>
  );
}
