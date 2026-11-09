/**
 * Result distribution chart using the FastNLR scientific palette.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import { toSeries } from "../lib/format";

const palette = ["#5cc5b2", "#9c88ff", "#eccc68", "#ff6b6b", "#55a7d8", "#8dbb75"];

export function DistributionChart({
  title,
  subtitle,
  values,
  horizontal = false,
}: {
  title: string;
  subtitle: string;
  values: Record<string, number>;
  horizontal?: boolean;
}) {
  const data = toSeries(values);
  const height = horizontal ? Math.max(250, data.length * 34 + 60) : 280;
  return (
    <article className="panel chart-card">
      <div className="chart-heading">
        <div>
          <h3>{title}</h3>
          <p>{subtitle}</p>
        </div>
        <span>{data.reduce((sum, item) => sum + item.value, 0)}</span>
      </div>
      {data.length ? (
        <div style={{ height }}>
          <ResponsiveContainer width="100%" height="100%">
            <BarChart
              data={data}
              layout={horizontal ? "vertical" : "horizontal"}
              margin={{ top: 8, right: 18, bottom: 8, left: horizontal ? 14 : 0 }}
            >
              <CartesianGrid
                stroke="#e7ecea"
                strokeDasharray="3 3"
                horizontal={!horizontal}
                vertical={horizontal}
              />
              {horizontal ? (
                <>
                  <XAxis
                    type="number"
                    allowDecimals={false}
                    tick={{ fontFamily: "Arial", fontSize: 11, fill: "#6b7c86" }}
                  />
                  <YAxis
                    type="category"
                    dataKey="name"
                    width={115}
                    tick={{ fontFamily: "Arial", fontSize: 11, fill: "#334852" }}
                  />
                </>
              ) : (
                <>
                  <XAxis
                    dataKey="name"
                    interval={0}
                    angle={data.length > 4 ? -24 : 0}
                    textAnchor={data.length > 4 ? "end" : "middle"}
                    height={data.length > 4 ? 62 : 36}
                    tick={{ fontFamily: "Arial", fontSize: 11, fill: "#334852" }}
                  />
                  <YAxis
                    allowDecimals={false}
                    tick={{ fontFamily: "Arial", fontSize: 11, fill: "#6b7c86" }}
                  />
                </>
              )}
              <Tooltip
                cursor={{ fill: "rgba(92, 197, 178, 0.08)" }}
                contentStyle={{
                  fontFamily: "Arial",
                  fontSize: 12,
                  borderRadius: 10,
                  border: "1px solid #dfe7e4",
                }}
              />
              <Bar dataKey="value" radius={horizontal ? [0, 5, 5, 0] : [5, 5, 0, 0]}>
                {data.map((entry, index) => (
                  <Cell
                    key={entry.name}
                    fill={palette[index % palette.length]}
                  />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>
      ) : (
        <div className="chart-empty">No data available for this chart.</div>
      )}
    </article>
  );
}
