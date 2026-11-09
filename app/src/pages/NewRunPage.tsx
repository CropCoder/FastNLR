/**
 * Genome and protein annotation setup form.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { zodResolver } from "@hookform/resolvers/zod";
import {
  ChevronDown,
  ChevronRight,
  FileSearch,
  FolderOpen,
  Gauge,
  Play,
  Plus,
  ServerCog,
  Trash2,
  UploadCloud,
} from "lucide-react";
import { useEffect, useState, type ReactNode } from "react";
import { useForm } from "react-hook-form";
import { useNavigate } from "react-router-dom";
import { z } from "zod";
import { PageHeader } from "../components/PageHeader";
import {
  hasExistingResults,
  inspectInput,
  listenForFileDrop,
  pickDirectory,
  pickFasta,
  pickFiles,
  startRun,
} from "../lib/bridge";
import { formatBytes } from "../lib/format";
import { useAppStore } from "../store/appStore";
import type {
  GenomeDraft,
  InputInfo,
  InputMode,
  ProteinDraft,
  RunDraft,
  RunHistoryEntry,
} from "../types";

interface FormValues {
  mode: InputMode;
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
  batch: number;
  hmmFiles: string[];
}

const probability = z
  .string()
  .trim()
  .refine((value) => {
    const parsed = Number(value);
    return Number.isFinite(parsed) && parsed > 0 && parsed <= 1;
  }, "Enter a probability greater than 0 and no more than 1");

const formSchema = z
  .object({
    mode: z.enum(["genome", "protein"]),
    inputPath: z.string().min(1, "Select an input FASTA"),
    outputDir: z.string().min(1, "Select an output directory"),
    threads: z.number().int().min(1).max(512),
    flank: z.number().int().min(0).max(1_000_000),
    seqsPerThread: z.number().int().min(1),
    checkpointDir: z.string(),
    motifAcceptP: probability,
    motifPrelimP: probability,
    tirAcceptP: probability,
    flexibleSeed: z.string(),
    noFlexibleSeed: z.boolean(),
    requirePloop: z.boolean(),
    relaxedSeed: z.string(),
    keepPartialLoci: z.boolean(),
    motFile: z.string(),
    storeFile: z.string(),
    motifCategories: z.string(),
    extraSeeds: z.string(),
    extraSignatures: z.string(),
    batch: z.number().int().min(1),
    hmmFiles: z.array(z.string()),
  })
  .superRefine((values, context) => {
    if (Number(values.motifPrelimP) < Number(values.motifAcceptP)) {
      context.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["motifPrelimP"],
        message: "Prefilter probability must be greater than or equal to the accept probability",
      });
    }
    if (values.relaxedSeed && !values.noFlexibleSeed) {
      context.addIssue({
        code: z.ZodIssueCode.custom,
        path: ["relaxedSeed"],
        message: "Relaxed seeding requires flexible seeding to be disabled",
      });
    }
  });

function baseDefaults(): FormValues {
  return {
    mode: "genome",
    inputPath: "",
    outputDir: "",
    threads: 1,
    flank: 2000,
    seqsPerThread: 1000,
    checkpointDir: "",
    motifAcceptP: "0.00001",
    motifPrelimP: "0.0001",
    tirAcceptP: "0.0000001",
    flexibleSeed: "3",
    noFlexibleSeed: false,
    requirePloop: false,
    relaxedSeed: "",
    keepPartialLoci: false,
    motFile: "",
    storeFile: "",
    motifCategories: "",
    extraSeeds: "",
    extraSignatures: "",
    batch: 256,
    hmmFiles: [],
  };
}

function draftToForm(draft: RunDraft): FormValues {
  const defaults = baseDefaults();
  if (draft.mode === "protein") {
    return {
      ...defaults,
      mode: "protein",
      inputPath: draft.inputPath,
      outputDir: draft.outputDir,
      threads: draft.threads,
      batch: draft.batch,
      motifAcceptP: draft.motifAcceptP,
      motifPrelimP: draft.motifPrelimP,
      hmmFiles: draft.hmmFiles,
      motFile: draft.motFile,
      storeFile: draft.storeFile,
    };
  }
  return {
    ...defaults,
    ...draft,
  };
}

function formToDraft(values: FormValues): RunDraft {
  if (values.mode === "protein") {
    const draft: ProteinDraft = {
      mode: "protein",
      inputPath: values.inputPath,
      outputDir: values.outputDir,
      threads: values.threads,
      batch: values.batch,
      motifAcceptP: values.motifAcceptP,
      motifPrelimP: values.motifPrelimP,
      hmmFiles: values.hmmFiles,
      motFile: values.motFile,
      storeFile: values.storeFile,
    };
    return draft;
  }
  const draft: GenomeDraft = {
    mode: "genome",
    inputPath: values.inputPath,
    outputDir: values.outputDir,
    threads: values.threads,
    flank: values.flank,
    seqsPerThread: values.seqsPerThread,
    checkpointDir: values.checkpointDir,
    motifAcceptP: values.motifAcceptP,
    motifPrelimP: values.motifPrelimP,
    tirAcceptP: values.tirAcceptP,
    flexibleSeed: values.flexibleSeed,
    noFlexibleSeed: values.noFlexibleSeed,
    requirePloop: values.requirePloop,
    relaxedSeed: values.relaxedSeed,
    keepPartialLoci: values.keepPartialLoci,
    motFile: values.motFile,
    storeFile: values.storeFile,
    motifCategories: values.motifCategories,
    extraSeeds: values.extraSeeds,
    extraSignatures: values.extraSignatures,
  };
  return draft;
}

function makeId() {
  return typeof crypto.randomUUID === "function"
    ? crypto.randomUUID()
    : `${Date.now()}-${Math.random().toString(16).slice(2)}`;
}

function defaultOutputPath(inputPath: string, prefix: string) {
  const separator = inputPath.includes("\\") ? "\\" : "/";
  const index = Math.max(inputPath.lastIndexOf("/"), inputPath.lastIndexOf("\\"));
  const parent = index >= 0 ? inputPath.slice(0, index) : ".";
  return `${parent}${separator}${prefix}-fastnlr`;
}

export function NewRunPage() {
  const navigate = useNavigate();
  const runtime = useAppStore((state) => state.runtime);
  const settings = useAppStore((state) => state.settings);
  const activeRun = useAppStore((state) => state.activeRun);
  const history = useAppStore((state) => state.history);
  const pendingDraft = useAppStore((state) => state.pendingDraft);
  const setPendingDraft = useAppStore((state) => state.setPendingDraft);
  const setActiveRun = useAppStore((state) => state.setActiveRun);
  const addHistory = useAppStore((state) => state.addHistory);
  const [inputInfo, setInputInfo] = useState<InputInfo | null>(null);
  const [advanced, setAdvanced] = useState(Boolean(pendingDraft));
  const [submitError, setSubmitError] = useState<string | null>(null);
  const [starting, setStarting] = useState(false);
  const activeEntry = activeRun
    ? history.find((entry) => entry.id === activeRun.historyId)
    : null;
  const runBusy =
    activeEntry?.status === "Pending" || activeEntry?.status === "Running";

  const defaults = pendingDraft
    ? draftToForm(pendingDraft)
    : {
        ...baseDefaults(),
        outputDir: settings.defaultOutputDir,
        threads: settings.defaultThreads || runtime?.logicalCpus || 1,
      };
  const {
    register,
    handleSubmit,
    watch,
    setValue,
    getValues,
    formState: { errors },
  } = useForm<FormValues>({
    resolver: zodResolver(formSchema),
    defaultValues: defaults,
  });

  const mode = watch("mode");
  const hmmFiles = watch("hmmFiles");
  const noFlexibleSeed = watch("noFlexibleSeed");
  const inputPath = watch("inputPath");

  useEffect(() => {
    if (pendingDraft) {
      setPendingDraft(null);
      if (pendingDraft.inputPath) {
        void inspectInput(pendingDraft.inputPath, pendingDraft.mode)
          .then(setInputInfo)
          .catch(() => setInputInfo(null));
      }
    }
  }, [pendingDraft, setPendingDraft]);

  useEffect(() => {
    let dispose: () => void = () => undefined;
    void listenForFileDrop((path) => {
      void applyInput(path, mode);
    }).then((unlisten) => {
      dispose = unlisten;
    });
    return () => dispose();
  }, [mode]);

  async function applyInput(path: string, nextMode = mode) {
    try {
      const info = await inspectInput(path, nextMode);
      setInputInfo(info);
      setValue("inputPath", info.path, { shouldValidate: true });
      if (!getValues("outputDir")) {
        setValue(
          "outputDir",
          settings.defaultOutputDir || defaultOutputPath(info.path, info.prefix),
          { shouldValidate: true },
        );
      }
      setSubmitError(null);
    } catch (error) {
      setSubmitError(String(error));
    }
  }

  async function chooseInput() {
    const path = await pickFasta(mode);
    if (path) await applyInput(path);
  }

  async function chooseOutput() {
    const path = await pickDirectory("Select output directory");
    if (path) setValue("outputDir", path, { shouldValidate: true });
  }

  async function chooseModelFiles(field: "motFile" | "storeFile") {
    const files = await pickFiles(field === "motFile" ? "Select mot.txt" : "Select store.txt");
    if (files[0]) setValue(field, files[0], { shouldValidate: true });
  }

  async function chooseHmmFiles() {
    const files = await pickFiles("Select HMM model files");
    if (!files.length) return;
    setValue("hmmFiles", Array.from(new Set([...getValues("hmmFiles"), ...files])), {
      shouldValidate: true,
    });
  }

  function switchMode(nextMode: InputMode) {
    setValue("mode", nextMode);
    setValue("inputPath", "");
    setInputInfo(null);
    setSubmitError(null);
  }

  async function submit(values: FormValues) {
    setStarting(true);
    setSubmitError(null);
    try {
      const info = await inspectInput(values.inputPath, values.mode);
      setInputInfo(info);
      const draft = formToDraft(values);
      if (await hasExistingResults(draft)) {
        const proceed = window.confirm(
          `Results for "${info.prefix}" already exist in this output directory. Continue and replace the generated files?`,
        );
        if (!proceed) return;
      }
      const jobId = await startRun(draft);
      const startedAt = new Date().toISOString();
      const historyId = makeId();
      const entry: RunHistoryEntry = {
        id: historyId,
        mode: values.mode,
        inputPath: values.inputPath,
        inputName: info.name,
        outputDir: values.outputDir,
        prefix: info.prefix,
        status: "Pending",
        phase: "Pending",
        startedAt,
        finishedAt: null,
        scanned: 0,
        estimatedTotal: 0,
        resultCount: null,
        error: null,
        draft,
      };
      addHistory(entry);
      setActiveRun({
        jobId,
        historyId,
        mode: values.mode,
        inputName: info.name,
        prefix: info.prefix,
        outputDir: values.outputDir,
        startedAt,
      });
      navigate("/run");
    } catch (error) {
      setSubmitError(String(error));
    } finally {
      setStarting(false);
    }
  }

  return (
    <div className="page">
      <PageHeader
        eyebrow="Annotation workspace"
        title="Launch a new NLR analysis"
        description="Run the same high-speed Rust pipeline as the FastNLR CLI, with validated parameters and live progress."
        actions={
          <div className="runtime-chip">
            <ServerCog size={16} />
            {runtime?.logicalCpus ?? 1} logical CPUs
          </div>
        }
      />

      <form onSubmit={handleSubmit(submit)} className="run-form">
        {runBusy ? (
          <div className="alert warning">
            An annotation run is already active. Open the run monitor or wait for it to
            finish before starting another analysis.
          </div>
        ) : null}

        <section className="panel mode-panel">
          <div className="mode-switch" role="tablist" aria-label="Input mode">
            <button
              type="button"
              role="tab"
              aria-selected={mode === "genome"}
              className={mode === "genome" ? "selected" : ""}
              onClick={() => switchMode("genome")}
            >
              <span className="mode-index">01</span>
              <span>
                <strong>Genome</strong>
                <small>Six-frame DNA annotation</small>
              </span>
            </button>
            <button
              type="button"
              role="tab"
              aria-selected={mode === "protein"}
              className={mode === "protein" ? "selected" : ""}
              onClick={() => switchMode("protein")}
            >
              <span className="mode-index">02</span>
              <span>
                <strong>Protein</strong>
                <small>HMM domain classification</small>
              </span>
            </button>
          </div>
        </section>

        <section className="panel">
          <div className="section-heading">
            <span className="section-number">1</span>
            <div>
              <h2>Input and destination</h2>
              <p>Select the FASTA file and the folder that will receive all generated files.</p>
            </div>
          </div>

          <div className="form-grid two-column">
            <div className="field full">
              <label htmlFor="input-path">Input FASTA</label>
              <button type="button" className="file-drop" onClick={chooseInput}>
                <UploadCloud size={25} />
                <span>
                  <strong>
                    {inputInfo?.name || (inputPath ? inputPath.split(/[\\/]/).pop() : "Choose or drop a FASTA file")}
                  </strong>
                  <small>
                    {inputInfo
                      ? `${formatBytes(inputInfo.sizeBytes)}${inputInfo.compressed ? " · gzip" : ""} · prefix ${inputInfo.prefix}`
                      : mode === "genome"
                        ? ".fa, .fasta, .fna, .fas, .ffn, or .gz"
                        : ".fa, .fasta, .faa, .fas, or .gz"}
                  </small>
                </span>
                <ChevronRight size={18} />
              </button>
              <input type="hidden" {...register("inputPath")} />
              {errors.inputPath ? <p className="field-error">{errors.inputPath.message}</p> : null}
            </div>

            <div className="field full">
              <label htmlFor="output-dir">Output directory</label>
              <div className="input-with-button">
                <input
                  id="output-dir"
                  placeholder="Choose a directory for FastNLR results"
                  {...register("outputDir")}
                />
                <button type="button" className="btn secondary" onClick={chooseOutput}>
                  <FolderOpen size={17} />
                  Browse
                </button>
              </div>
              {errors.outputDir ? <p className="field-error">{errors.outputDir.message}</p> : null}
            </div>

            <div className="field">
              <label htmlFor="threads">Threads</label>
              <input
                id="threads"
                type="number"
                min={1}
                max={512}
                {...register("threads", { valueAsNumber: true })}
              />
              <small>Use no more than the physical CPU core count for best throughput.</small>
              {errors.threads ? <p className="field-error">{errors.threads.message}</p> : null}
            </div>

            <div className="field">
              <label>Output prefix</label>
              <div className="readonly-field">
                <FileSearch size={16} />
                {inputInfo?.prefix || "Derived from the input filename"}
              </div>
            </div>
          </div>
        </section>

        <section className="panel">
          <button
            type="button"
            className="advanced-toggle"
            onClick={() => setAdvanced((value) => !value)}
            aria-expanded={advanced}
          >
            <span className="section-heading compact">
              <span className="section-number">2</span>
              <span>
                <h2>Advanced parameters</h2>
                <p>Motif thresholds, seeding, checkpoints, and external model files.</p>
              </span>
            </span>
            <ChevronDown className={advanced ? "rotate" : ""} size={20} />
          </button>

          {advanced ? (
            <div className="advanced-content">
              <div className="form-grid three-column">
                {mode === "genome" ? (
                  <>
                    <Field label="Flank bases" error={errors.flank?.message}>
                      <input type="number" min={0} {...register("flank", { valueAsNumber: true })} />
                    </Field>
                    <Field label="Fragments per batch" error={errors.seqsPerThread?.message}>
                      <input
                        type="number"
                        min={1}
                        {...register("seqsPerThread", { valueAsNumber: true })}
                      />
                    </Field>
                    <Field label="Checkpoint directory">
                      <input placeholder="Optional" {...register("checkpointDir")} />
                    </Field>
                    <Field label="Motif accept p" error={errors.motifAcceptP?.message}>
                      <input {...register("motifAcceptP")} />
                    </Field>
                    <Field label="Motif prefilter p" error={errors.motifPrelimP?.message}>
                      <input {...register("motifPrelimP")} />
                    </Field>
                    <Field label="TIR accept p" error={errors.tirAcceptP?.message}>
                      <input {...register("tirAcceptP")} />
                    </Field>
                    <Field label="Flexible seed minimum">
                      <input type="number" min={1} {...register("flexibleSeed")} />
                    </Field>
                    <Field label="Relaxed seed minimum" error={errors.relaxedSeed?.message}>
                      <input
                        type="number"
                        min={1}
                        placeholder="Optional"
                        disabled={!noFlexibleSeed}
                        {...register("relaxedSeed")}
                      />
                    </Field>
                  </>
                ) : (
                  <>
                    <Field label="Proteins per batch" error={errors.batch?.message}>
                      <input type="number" min={1} {...register("batch", { valueAsNumber: true })} />
                    </Field>
                    <Field label="Motif accept p" error={errors.motifAcceptP?.message}>
                      <input {...register("motifAcceptP")} />
                    </Field>
                    <Field label="Motif prefilter p" error={errors.motifPrelimP?.message}>
                      <input {...register("motifPrelimP")} />
                    </Field>
                  </>
                )}
              </div>

              {mode === "genome" ? (
                <div className="check-grid">
                  <label className="check-card">
                    <input type="checkbox" {...register("noFlexibleSeed")} />
                    <span>
                      <strong>Disable flexible seeding</strong>
                      <small>Use only exact motif-id seed combinations.</small>
                    </span>
                  </label>
                  <label className="check-card">
                    <input type="checkbox" {...register("requirePloop")} />
                    <span>
                      <strong>Require P-loop</strong>
                      <small>Higher precision for flexible seeds.</small>
                    </span>
                  </label>
                  <label className="check-card">
                    <input type="checkbox" {...register("keepPartialLoci")} />
                    <span>
                      <strong>Keep partial loci</strong>
                      <small>Restore legacy permissive filtering.</small>
                    </span>
                  </label>
                </div>
              ) : (
                <div className="hmm-section">
                  <div className="inline-heading">
                    <div>
                      <h3>Custom HMM files</h3>
                      <p>Optional models are added to the five embedded Pfam profiles.</p>
                    </div>
                    <button type="button" className="btn secondary small" onClick={chooseHmmFiles}>
                      <Plus size={16} />
                      Add HMM
                    </button>
                  </div>
                  {hmmFiles.length ? (
                    <div className="file-list">
                      {hmmFiles.map((file) => (
                        <div className="file-row" key={file}>
                          <span title={file}>{file}</span>
                          <button
                            type="button"
                            aria-label={`Remove ${file}`}
                            onClick={() =>
                              setValue(
                                "hmmFiles",
                                hmmFiles.filter((value) => value !== file),
                              )
                            }
                          >
                            <Trash2 size={15} />
                          </button>
                        </div>
                      ))}
                    </div>
                  ) : (
                    <p className="muted-copy">Using the embedded NB-ARC, TIR, CC, and LRR profiles.</p>
                  )}
                </div>
              )}

              <div className="form-grid two-column advanced-files">
                <Field label="Custom mot.txt">
                  <div className="input-with-button">
                    <input placeholder="Use embedded profile" {...register("motFile")} />
                    <button type="button" className="btn ghost" onClick={() => chooseModelFiles("motFile")}>
                      Browse
                    </button>
                  </div>
                </Field>
                <Field label="Custom store.txt">
                  <div className="input-with-button">
                    <input placeholder="Use embedded profile" {...register("storeFile")} />
                    <button type="button" className="btn ghost" onClick={() => chooseModelFiles("storeFile")}>
                      Browse
                    </button>
                  </div>
                </Field>
                {mode === "genome" ? (
                  <>
                    <Field label="Motif categories">
                      <input placeholder="21=CC,22=CC,25=TIR" {...register("motifCategories")} />
                    </Field>
                    <Field label="Extra seed combinations">
                      <input placeholder="21,4;22,4" {...register("extraSeeds")} />
                    </Field>
                    <Field label="Extra signatures">
                      <input placeholder="21,4;22,4" {...register("extraSignatures")} />
                    </Field>
                  </>
                ) : null}
              </div>
            </div>
          ) : null}
        </section>

        {submitError ? <div className="alert error">{submitError}</div> : null}

        <div className="launch-bar">
          <div>
            <Gauge size={18} />
            <span>
              <strong>{mode === "genome" ? "Genome annotation" : "Protein classification"}</strong>
              <small>
                {mode === "genome"
                  ? "Outputs include loci FASTA, GFF, BED, statistics, plots, and summary JSON."
                  : "Outputs include domain hits, NLR candidates, and a run summary."}
              </small>
            </span>
          </div>
          <button
            type="submit"
            className="btn primary large"
            disabled={starting || runBusy}
          >
            <Play size={18} fill="currentColor" />
            {starting ? "Starting..." : "Start annotation"}
          </button>
        </div>
      </form>
    </div>
  );
}

function Field({
  label,
  error,
  children,
}: {
  label: string;
  error?: string;
  children: ReactNode;
}) {
  return (
    <div className="field">
      <label>{label}</label>
      {children}
      {error ? <p className="field-error">{error}</p> : null}
    </div>
  );
}

export const newRunTestHelpers = {
  baseDefaults,
  draftToForm,
  formToDraft,
};
