/**
 * Shared TypeScript contracts for Tauri commands and desktop state.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
export type InputMode = "genome" | "protein";

export type JobStatus =
  | "Pending"
  | "Running"
  | "Completed"
  | "Failed"
  | "Canceled";

export type JobPhase =
  | "Pending"
  | "Scanning"
  | "Writing"
  | "Completed"
  | "Failed"
  | "Canceled";

export type HistoryStatus = JobStatus | "Interrupted";

export interface RuntimeInfo {
  version: string;
  os: string;
  arch: string;
  logicalCpus: number;
}

export interface InputInfo {
  path: string;
  name: string;
  prefix: string;
  sizeBytes: number;
  compressed: boolean;
}

export interface JobSnapshot {
  id: number;
  kind: "Genome" | "Protein";
  status: JobStatus;
  phase: JobPhase;
  scanned: number;
  estimated_total: number;
  error: string | null;
}

export interface OutputFile {
  name: string;
  path: string;
  sizeBytes: number;
  isDir: boolean;
}

export interface GenomeLocus {
  name: string;
  seq: string;
  start: number;
  end: number;
  strand: string;
  class: string;
  family: string;
  has_stop: boolean;
  motif_count: number;
}

export interface ChromosomeStats {
  motifs: number;
  nlrs: number;
}

export interface GenomeSummary {
  sequence_count: number;
  motif_total: number;
  nlr_total: number;
  per_chromosome: Record<string, ChromosomeStats>;
  motif_counts: Record<string, number>;
  nlr_type_counts: Record<string, number>;
  nlr_family_counts: Record<string, number>;
  loci: GenomeLocus[];
}

export interface ProteinRecord {
  proteinId: string;
  length: number;
  class: string;
  nterm: string;
  hasNbarc: boolean;
  hasTir: boolean;
  hasLrr: boolean;
  domainCount: number;
  domains: string;
  family: string;
}

export interface ProteinDomain {
  proteinId: string;
  source: string;
  model: string;
  modelName: string;
  category: string;
  start: number;
  end: number;
  bitscore: number | null;
  pvalue: number | null;
  passed: boolean;
}

export interface ProteinSummary {
  proteinsScanned: number;
  nlrTotal: number;
  elapsed: string | null;
  modelLabels: string[];
  typeCounts: Record<string, number>;
  familyCounts: Record<string, number>;
  records: ProteinRecord[];
  domains: ProteinDomain[];
}

export interface RunResults {
  mode: InputMode;
  prefix: string;
  outputDir: string;
  outputs: OutputFile[];
  genome: GenomeSummary | null;
  protein: ProteinSummary | null;
}

export interface GenomeDraft {
  mode: "genome";
  inputPath: string;
  outputDir: string;
  threads: number;
  flank: number;
  seqsPerThread: number;
  checkpointDir: string;
  motifAcceptP: string;
  motifPrelimP: string;
  tirAcceptP: string;
  flexibleSeed: string;
  noFlexibleSeed: boolean;
  requirePloop: boolean;
  relaxedSeed: string;
  keepPartialLoci: boolean;
  motFile: string;
  storeFile: string;
  motifCategories: string;
  extraSeeds: string;
  extraSignatures: string;
}

export interface ProteinDraft {
  mode: "protein";
  inputPath: string;
  outputDir: string;
  threads: number;
  batch: number;
  motifAcceptP: string;
  motifPrelimP: string;
  hmmFiles: string[];
  motFile: string;
  storeFile: string;
}

export type RunDraft = GenomeDraft | ProteinDraft;

export interface AppSettings {
  defaultThreads: number;
  defaultOutputDir: string;
}

export interface RunHistoryEntry {
  id: string;
  mode: InputMode;
  inputPath: string;
  inputName: string;
  outputDir: string;
  prefix: string;
  status: HistoryStatus;
  phase: JobPhase | "Interrupted";
  startedAt: string;
  finishedAt: string | null;
  scanned: number;
  estimatedTotal: number;
  resultCount: number | null;
  error: string | null;
  draft: RunDraft;
}

export interface ActiveRun {
  jobId: number;
  historyId: string;
  mode: InputMode;
  inputName: string;
  prefix: string;
  outputDir: string;
  startedAt: string;
}

export interface PersistedAppState {
  settings: AppSettings;
  history: RunHistoryEntry[];
}
