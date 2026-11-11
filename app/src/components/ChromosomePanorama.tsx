/**
 * Interactive chromosome panorama with locus blocks and smoothed NLR density.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Check,
  Copy,
  Map as MapIcon,
  MousePointer2,
  RotateCcw,
  X,
} from "lucide-react";
import {
  useEffect,
  useMemo,
  useRef,
  useState,
  type KeyboardEvent,
} from "react";
import type { ChromosomeStats, GenomeLocus } from "../types";

const FAMILY_COLORS: Record<string, string> = {
  TNL: "#2a9d82",
  CNL: "#7b6ce4",
  RNL: "#e36565",
  "TIR-only": "#c89a20",
  "NB-only": "#4e8fbd",
  Other: "#7d8785",
};

const FALLBACK_COLORS = [
  "#2a9d82",
  "#7b6ce4",
  "#e36565",
  "#c89a20",
  "#4e8fbd",
  "#8d7b69",
];

interface PanoramaChromosome {
  name: string;
  length: number;
  loci: GenomeLocus[];
}

interface DensityPoint {
  x: number;
  value: number;
}

interface PathPoint {
  x: number;
  y: number;
}

export function ChromosomePanorama({
  chromosomes,
  loci,
}: {
  chromosomes: Record<string, ChromosomeStats>;
  loci: GenomeLocus[];
}) {
  const hostRef = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(1080);
  const [hiddenFamilies, setHiddenFamilies] = useState<Set<string>>(new Set());
  const [hoveredLocus, setHoveredLocus] = useState<{
    locus: GenomeLocus;
    x: number;
    y: number;
  } | null>(null);
  const [selectedLocus, setSelectedLocus] = useState<GenomeLocus | null>(null);
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    const observer = new ResizeObserver((entries) => {
      const next = entries[0]?.contentRect.width;
      if (next) setWidth(Math.max(620, Math.floor(next)));
    });
    observer.observe(host);
    return () => observer.disconnect();
  }, []);

  useEffect(() => {
    if (!selectedLocus) return;
    const onKeyDown = (event: globalThis.KeyboardEvent) => {
      if (event.key === "Escape") setSelectedLocus(null);
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [selectedLocus]);

  const families = useMemo(() => {
    const counts = new Map<string, number>();
    for (const locus of loci) {
      counts.set(locus.family, (counts.get(locus.family) ?? 0) + 1);
    }
    return [...counts.entries()].sort(
      (left, right) => right[1] - left[1] || left[0].localeCompare(right[0]),
    );
  }, [loci]);

  const resolved = useMemo(() => {
    const bySequence = new Map<string, GenomeLocus[]>();
    for (const locus of loci) {
      const group = bySequence.get(locus.seq) ?? [];
      group.push(locus);
      bySequence.set(locus.seq, group);
    }
    const names = new Set([
      ...Object.keys(chromosomes),
      ...bySequence.keys(),
    ]);
    return [...names]
      .map((name): PanoramaChromosome => {
        const sequenceLoci = bySequence.get(name) ?? [];
        const observedEnd = sequenceLoci.reduce(
          (maximum, locus) => Math.max(maximum, locus.end),
          0,
        );
        const declaredLength = chromosomes[name]?.length ?? 0;
        const length =
          declaredLength > 0
            ? declaredLength
            : observedEnd > 0
              ? Math.ceil(observedEnd * 1.04)
              : 1_000_000;
        return { name, length, loci: sequenceLoci };
      })
      .sort(
        (left, right) =>
          right.length - left.length || left.name.localeCompare(right.name),
      );
  }, [chromosomes, loci]);

  const visibleLoci = useMemo(
    () => loci.filter((locus) => !hiddenFamilies.has(locus.family)),
    [hiddenFamilies, loci],
  );

  const visibleBySequence = useMemo(() => {
    const bySequence = new Map<string, GenomeLocus[]>();
    for (const locus of visibleLoci) {
      const group = bySequence.get(locus.seq) ?? [];
      group.push(locus);
      bySequence.set(locus.seq, group);
    }
    return bySequence;
  }, [visibleLoci]);

  const plot = useMemo(() => {
    const left = 148;
    const right = 22;
    const top = 58;
    const bottom = 46;
    const rowHeight = 116;
    const densityHeight = 32;
    const plotWidth = Math.max(360, width - left - right);
    const plotHeight = Math.max(rowHeight, resolved.length * rowHeight);
    const height = top + plotHeight + bottom;
    const maxLength = Math.max(1, ...resolved.map((item) => item.length));
    const maxMb = maxLength / 1_000_000;
    const stepMb = niceStep(maxMb / 6);
    const ticks: number[] = [];
    for (let value = 0; value <= maxMb + stepMb * 0.25; value += stepMb) {
      ticks.push(value);
    }
    const xForBp = (bp: number) =>
      left + (Math.max(0, Math.min(bp, maxLength)) / maxLength) * plotWidth;
    const rawDensity = resolved.map((chromosome) =>
      buildDensity(
        visibleBySequence.get(chromosome.name) ?? [],
        chromosome.length,
      ),
    );
    const maxDensity = Math.max(
      0.001,
      ...rawDensity.flatMap((values) => values),
    );
    return {
      left,
      right,
      top,
      bottom,
      rowHeight,
      densityHeight,
      plotWidth,
      height,
      maxLength,
      ticks,
      xForBp,
      rawDensity,
      maxDensity,
    };
  }, [resolved, visibleBySequence, width]);

  function toggleFamily(family: string) {
    setHiddenFamilies((current) => {
      const next = new Set(current);
      if (next.has(family)) next.delete(family);
      else next.add(family);
      return next;
    });
  }

  function openLocus(locus: GenomeLocus) {
    setSelectedLocus(locus);
    setCopied(false);
  }

  async function copyCoordinates() {
    if (!selectedLocus) return;
    const coordinates = `${selectedLocus.seq}:${selectedLocus.start}-${selectedLocus.end}`;
    try {
      await navigator.clipboard.writeText(coordinates);
      setCopied(true);
    } catch {
      setCopied(false);
    }
  }

  function onLocusKeyDown(
    event: KeyboardEvent<SVGGElement>,
    locus: GenomeLocus,
  ) {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      openLocus(locus);
    }
  }

  if (!resolved.length) {
    return (
      <div className="panorama-empty">
        <MapIcon size={26} />
        <strong>No chromosome coordinates available</strong>
        <span>Run genome mode to generate the physical-position panorama.</span>
      </div>
    );
  }

  return (
    <div className="chromosome-panorama" ref={hostRef}>
      <div className="panorama-toolbar">
        <div className="panorama-hint">
          <MousePointer2 size={15} />
          <span>
            Smooth curves show local NLR density. Select a locus block for full details.
          </span>
        </div>
        <div className="family-filter" aria-label="NLR family filters">
          {families.map(([family, count]) => {
            const hidden = hiddenFamilies.has(family);
            return (
              <button
                key={family}
                type="button"
                className={hidden ? "muted" : ""}
                onClick={() => toggleFamily(family)}
                title={hidden ? `Show ${family}` : `Hide ${family}`}
                aria-pressed={!hidden}
              >
                <span style={{ background: familyColor(family) }} />
                {family}
                <small>{count}</small>
              </button>
            );
          })}
          {hiddenFamilies.size ? (
            <button
              type="button"
              className="filter-reset"
              onClick={() => setHiddenFamilies(new Set())}
              title="Show all NLR families"
            >
              <RotateCcw size={13} />
              Reset
            </button>
          ) : null}
        </div>
      </div>

      <div className="panorama-scroll">
        <svg
          className="panorama-svg"
          width={width}
          height={plot.height}
          role="img"
          aria-label="Chromosome panorama of NLR loci and density"
        >
          <g className="panorama-axis">
            {plot.ticks.map((tick) => {
              const x = plot.xForBp(tick * 1_000_000);
              return (
                <g key={tick}>
                  <line
                    x1={x}
                    x2={x}
                    y1={plot.top - 20}
                    y2={plot.height - plot.bottom}
                  />
                  <text
                    x={x}
                    y={plot.top - 28}
                    textAnchor="middle"
                  >
                    {formatTick(tick)}
                  </text>
                </g>
              );
            })}
            <text
              className="axis-title"
              x={plot.left + plot.plotWidth / 2}
              y={plot.height - 12}
              textAnchor="middle"
            >
              Physical position (Mb)
            </text>
          </g>

          {resolved.map((chromosome, rowIndex) => {
            const rowTop = plot.top + rowIndex * plot.rowHeight;
            const trackY = rowTop + 72;
            const densityBaseline = rowTop + 44;
            const chromosomeLoci =
              visibleBySequence.get(chromosome.name) ?? [];
            const density = plot.rawDensity[rowIndex] ?? [];
            const points = density.map((value, index): DensityPoint => ({
              x:
                plot.left +
                (index / Math.max(1, density.length - 1)) * plot.plotWidth,
              value,
            }));
            const linePath = smoothPath(
              points.map((point) => ({
                x: point.x,
                y:
                  densityBaseline -
                  (point.value / plot.maxDensity) * plot.densityHeight,
              })),
            );
            const areaPath = `${linePath} L ${
              plot.left + plot.plotWidth
            } ${densityBaseline} L ${plot.left} ${densityBaseline} Z`;
            const trackWidth = plot.xForBp(chromosome.length) - plot.left;

            return (
              <g
                className="panorama-row"
                key={chromosome.name}
                style={{ animationDelay: `${Math.min(rowIndex * 35, 280)}ms` }}
              >
                <line
                  className="row-separator"
                  x1={0}
                  x2={width}
                  y1={rowTop}
                  y2={rowTop}
                />
                <text
                  className="chromosome-name"
                  x={plot.left - 14}
                  y={trackY + 1}
                  textAnchor="end"
                >
                  {truncateLabel(chromosome.name, 18)}
                  <title>{chromosome.name}</title>
                </text>
                <text
                  className="chromosome-length"
                  x={plot.left - 14}
                  y={trackY + 17}
                  textAnchor="end"
                >
                  {formatMb(chromosome.length)} Mb
                </text>

                <path className="density-area" d={areaPath} />
                <path className="density-line" d={linePath} />
                <rect
                  className="chromosome-track"
                  x={plot.left}
                  y={trackY}
                  width={Math.max(2, trackWidth)}
                  height={14}
                  rx={2}
                />

                {chromosomeLoci.map((locus) => {
                  const startX = plot.xForBp(locus.start);
                  const endX = plot.xForBp(locus.end);
                  const blockWidth = Math.max(3, endX - startX);
                  const blockX = Math.min(
                    plot.left + trackWidth - 1,
                    Math.max(plot.left, startX),
                  );
                  const color = familyColor(locus.family);
                  return (
                    <g
                      key={locus.name}
                      className="locus-block"
                      role="button"
                      tabIndex={0}
                      aria-label={`${locus.name}, ${locus.family}, ${formatMb(
                        locus.start,
                      )} to ${formatMb(locus.end)} Mb`}
                      onClick={() => openLocus(locus)}
                      onKeyDown={(event) => onLocusKeyDown(event, locus)}
                      onMouseEnter={() =>
                        setHoveredLocus({
                          locus,
                          x: Math.min(
                            width - 210,
                            Math.max(plot.left, blockX - 80),
                          ),
                          y: rowTop + 7,
                        })
                      }
                      onMouseLeave={() => setHoveredLocus(null)}
                      onFocus={() =>
                        setHoveredLocus({
                          locus,
                          x: Math.min(
                            width - 210,
                            Math.max(plot.left, blockX - 80),
                          ),
                          y: rowTop + 7,
                        })
                      }
                      onBlur={() => setHoveredLocus(null)}
                    >
                      <rect
                        x={blockX}
                        y={trackY - 2}
                        width={Math.min(blockWidth, trackWidth)}
                        height={18}
                        rx={2}
                        fill={color}
                      >
                        <title>
                          {locus.name} · {locus.family} ·{" "}
                          {formatMb(locus.start)}-{formatMb(locus.end)} Mb
                        </title>
                      </rect>
                    </g>
                  );
                })}

                <text
                  className="density-label"
                  x={plot.left + plot.plotWidth}
                  y={rowTop + 17}
                  textAnchor="end"
                >
                  {chromosomeLoci.length} loci
                </text>
              </g>
            );
          })}

          {hoveredLocus ? (
            <g
              className="panorama-hover"
              transform={`translate(${hoveredLocus.x}, ${hoveredLocus.y})`}
              pointerEvents="none"
            >
              <rect width={202} height={48} rx={4} />
              <text x={10} y={18}>
                {truncateLabel(hoveredLocus.locus.name, 28)}
              </text>
              <circle
                cx={14}
                cy={33}
                r={3}
                fill={familyColor(hoveredLocus.locus.family)}
              />
              <text className="hover-meta" x={23} y={37}>
                {hoveredLocus.locus.family} · {hoveredLocus.locus.class}
              </text>
            </g>
          ) : null}
        </svg>
      </div>

      {selectedLocus ? (
        <div
          className="locus-modal-backdrop"
          role="presentation"
          onMouseDown={() => setSelectedLocus(null)}
        >
          <section
            className="locus-modal"
            role="dialog"
            aria-modal="true"
            aria-labelledby="locus-modal-title"
            onMouseDown={(event) => event.stopPropagation()}
          >
            <div className="locus-modal-head">
              <div>
                <span className="eyebrow">NLR locus detail</span>
                <h2 id="locus-modal-title">{selectedLocus.name}</h2>
              </div>
              <button
                type="button"
                className="icon-button"
                onClick={() => setSelectedLocus(null)}
                aria-label="Close locus details"
              >
                <X size={17} />
              </button>
            </div>

            <div className="locus-modal-family">
              <span style={{ background: familyColor(selectedLocus.family) }} />
              <strong>{selectedLocus.family}</strong>
              <span>{selectedLocus.class}</span>
            </div>

            <div className="locus-detail-grid">
              <Detail label="Chromosome" value={selectedLocus.seq} />
              <Detail
                label="Coordinates"
                value={`${formatBp(selectedLocus.start)}–${formatBp(
                  selectedLocus.end,
                )} bp`}
              />
              <Detail
                label="Physical position"
                value={`${formatMb(selectedLocus.start)}–${formatMb(
                  selectedLocus.end,
                )} Mb`}
              />
              <Detail
                label="Strand"
                value={selectedLocus.strand === "-" ? "Reverse (−)" : "Forward (+)"}
              />
              <Detail label="Motif count" value={String(selectedLocus.motif_count)} />
              <Detail label="Stop codon" value={selectedLocus.has_stop ? "Present" : "Absent"} />
            </div>

            <div className="locus-coordinate-box">
              <span>{selectedLocus.seq}:{selectedLocus.start}-{selectedLocus.end}</span>
              <button type="button" onClick={() => void copyCoordinates()}>
                {copied ? <Check size={14} /> : <Copy size={14} />}
                {copied ? "Copied" : "Copy"}
              </button>
            </div>

            <p className="locus-modal-note">
              Coordinates use the unflanked core locus interval from the genome assembly.
            </p>
          </section>
        </div>
      ) : null}
    </div>
  );
}

function Detail({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function familyColor(family: string) {
  if (FAMILY_COLORS[family]) return FAMILY_COLORS[family];
  let hash = 0;
  for (const character of family) hash = (hash * 31 + character.charCodeAt(0)) | 0;
  return FALLBACK_COLORS[Math.abs(hash) % FALLBACK_COLORS.length];
}

function buildDensity(loci: GenomeLocus[], length: number, bins = 96) {
  const values = Array.from({ length: bins }, () => 0);
  if (!loci.length || length <= 0) return values;
  const sigma = Math.max(1.2, bins * 0.025);
  for (const locus of loci) {
    const midpoint = (locus.start + locus.end) / 2;
    const center = (midpoint / length) * (bins - 1);
    for (let index = 0; index < bins; index += 1) {
      const distance = (index - center) / sigma;
      values[index] += Math.exp(-0.5 * distance * distance);
    }
  }
  return values;
}

function smoothPath(points: PathPoint[]) {
  if (!points.length) return "";
  if (points.length === 1) return `M ${points[0].x} ${points[0].y}`;
  let path = `M ${points[0].x.toFixed(2)} ${points[0].y.toFixed(2)}`;
  for (let index = 0; index < points.length - 1; index += 1) {
    const previous = points[Math.max(0, index - 1)];
    const current = points[index];
    const next = points[index + 1];
    const after = points[Math.min(points.length - 1, index + 2)];
    const cp1x = current.x + (next.x - previous.x) / 6;
    const cp1y = current.y + (next.y - previous.y) / 6;
    const cp2x = next.x - (after.x - current.x) / 6;
    const cp2y = next.y - (after.y - current.y) / 6;
    path += ` C ${cp1x.toFixed(2)} ${cp1y.toFixed(2)}, ${cp2x.toFixed(
      2,
    )} ${cp2y.toFixed(2)}, ${next.x.toFixed(2)} ${next.y.toFixed(2)}`;
  }
  return path;
}

function niceStep(value: number) {
  if (!Number.isFinite(value) || value <= 0) return 1;
  const exponent = Math.floor(Math.log10(value));
  const fraction = value / 10 ** exponent;
  const niceFraction =
    fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
  return niceFraction * 10 ** exponent;
}

function formatTick(mb: number) {
  if (mb >= 100) return Math.round(mb).toString();
  if (mb >= 10) return mb.toFixed(0);
  return mb.toFixed(1);
}

function formatMb(bp: number) {
  return (bp / 1_000_000).toFixed(bp >= 10_000_000 ? 1 : 2);
}

function formatBp(bp: number) {
  return Math.round(bp).toLocaleString();
}

function truncateLabel(value: string, maxLength: number) {
  return value.length > maxLength
    ? `${value.slice(0, maxLength - 1)}…`
    : value;
}

export const chromosomePanoramaTestHelpers = {
  buildDensity,
  familyColor,
  niceStep,
  smoothPath,
};
