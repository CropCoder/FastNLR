/**
 * Job and history status badge.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import type { HistoryStatus, JobPhase } from "../types";

export function StatusBadge({
  status,
  phase,
}: {
  status: HistoryStatus;
  phase?: JobPhase | "Interrupted";
}) {
  const label = status === "Running" && phase ? phase : status;
  const tone = status.toLowerCase();
  return <span className={`status-badge ${tone}`}>{label}</span>;
}
