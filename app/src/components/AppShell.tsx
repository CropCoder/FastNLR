/**
 * Desktop navigation shell and global active-job polling.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Activity,
  AlertCircle,
  Check,
  Dna,
  History,
  Info,
  Play,
  Settings,
} from "lucide-react";
import { useEffect } from "react";
import { NavLink, Outlet } from "react-router-dom";
import { getJob, loadResults } from "../lib/bridge";
import { useAppStore } from "../store/appStore";

const navigation = [
  { to: "/new", label: "New Run", icon: Play },
  { to: "/history", label: "History", icon: History },
  { to: "/settings", label: "Settings", icon: Settings },
  { to: "/about", label: "About", icon: Info },
];

export function AppShell() {
  const runtime = useAppStore((state) => state.runtime);
  const activeRun = useAppStore((state) => state.activeRun);
  const history = useAppStore((state) => state.history);
  const updateHistory = useAppStore((state) => state.updateHistory);
  const activeEntry = activeRun
    ? history.find((entry) => entry.id === activeRun.historyId)
    : null;
  const runComplete = activeEntry?.status === "Completed";
  const runEnded =
    activeEntry?.status === "Failed" || activeEntry?.status === "Canceled";

  useEffect(() => {
    if (!activeRun) return;
    let disposed = false;
    let timer: number | undefined;

    const poll = async () => {
      try {
        const snapshot = await getJob(activeRun.jobId);
        if (disposed) return;
        updateHistory(activeRun.historyId, {
          status: snapshot.status,
          phase: snapshot.phase,
          scanned: snapshot.scanned,
          estimatedTotal: snapshot.estimated_total,
          error: snapshot.error,
          finishedAt:
            snapshot.status === "Completed" ||
            snapshot.status === "Failed" ||
            snapshot.status === "Canceled"
              ? new Date().toISOString()
              : null,
        });

        if (snapshot.status === "Completed") {
          try {
            const results = await loadResults({
              mode: activeRun.mode,
              outputDir: activeRun.outputDir,
              prefix: activeRun.prefix,
            });
            if (disposed) return;
            updateHistory(activeRun.historyId, {
              resultCount:
                results.genome?.nlr_total ?? results.protein?.nlrTotal ?? 0,
            });
          } catch {
            // 结果页会给出更具体的读取错误，历史状态仍保留为 Completed。
          }
          return;
        }

        if (snapshot.status === "Failed" || snapshot.status === "Canceled") {
          return;
        }

        timer = window.setTimeout(poll, 700);
      } catch (error) {
        if (disposed) return;
        updateHistory(activeRun.historyId, {
          status: "Failed",
          phase: "Failed",
          error: String(error),
          finishedAt: new Date().toISOString(),
        });
      }
    };

    void poll();
    return () => {
      disposed = true;
      if (timer) window.clearTimeout(timer);
    };
  }, [
    activeRun?.jobId,
    activeRun?.historyId,
    activeRun?.mode,
    activeRun?.outputDir,
    activeRun?.prefix,
    updateHistory,
  ]);

  return (
    <div className="app-shell">
      <aside className="sidebar">
        <div className="sidebar-brand">
          <div className="brand-mark">
            <Dna size={22} strokeWidth={2.4} />
          </div>
          <div>
            <strong>FastNLR</strong>
            <span>Desktop</span>
          </div>
        </div>

        <nav className="sidebar-nav" aria-label="Main navigation">
          {navigation.map(({ to, label, icon: Icon }) => (
            <NavLink
              key={to}
              to={to}
              className={({ isActive }) =>
                `nav-item${isActive ? " active" : ""}`
              }
            >
              <Icon size={18} />
              <span>{label}</span>
              {to === "/new" && activeRun ? (
                <span className="nav-pulse" title="A run is active" />
              ) : null}
            </NavLink>
          ))}
        </nav>

        {activeRun && activeEntry ? (
          <NavLink
            to={runComplete ? `/results/${activeRun.historyId}` : "/run"}
            className={`active-run-chip${runComplete ? " complete" : ""}${runEnded ? " ended" : ""}`}
          >
            {runComplete ? (
              <Check size={15} />
            ) : runEnded ? (
              <AlertCircle size={15} />
            ) : (
              <Activity size={15} />
            )}
            <span>
              <strong>
                {runComplete
                  ? "Run complete"
                  : runEnded
                    ? "Run ended"
                    : "Run in progress"}
              </strong>
              <small>{activeRun.inputName}</small>
            </span>
          </NavLink>
        ) : (
          <div className="sidebar-ready">
            <span className="ready-dot" />
            Ready for annotation
          </div>
        )}

        <div className="sidebar-footer">
          <span>v{runtime?.version ?? "1.8.0"}</span>
          <span>
            {runtime?.os ?? "desktop"} / {runtime?.arch ?? "native"}
          </span>
        </div>
      </aside>
      <main className="main-content">
        <Outlet />
      </main>
    </div>
  );
}
