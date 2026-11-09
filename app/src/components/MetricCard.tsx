/**
 * Result overview metric card.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import type { ReactNode } from "react";

export function MetricCard({
  label,
  value,
  detail,
  icon,
  accent = "mint",
}: {
  label: string;
  value: string | number;
  detail?: string;
  icon?: ReactNode;
  accent?: "mint" | "sun" | "iris" | "coral";
}) {
  return (
    <article className={`metric-card ${accent}`}>
      <div className="metric-top">
        <span>{label}</span>
        {icon ? <span className="metric-icon">{icon}</span> : null}
      </div>
      <strong>{value}</strong>
      {detail ? <small>{detail}</small> : null}
    </article>
  );
}
