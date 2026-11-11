/**
 * Desktop defaults and local data settings.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  Check,
  Cpu,
  FolderOpen,
  HardDrive,
  RotateCcw,
  Save,
  Settings,
  Trash2,
} from "lucide-react";
import { useEffect, useState } from "react";
import { PageHeader } from "../components/PageHeader";
import { InfoHint } from "../components/InfoHint";
import { pickDirectory } from "../lib/bridge";
import { useAppStore } from "../store/appStore";

export function SettingsPage() {
  const runtime = useAppStore((state) => state.runtime);
  const settings = useAppStore((state) => state.settings);
  const saveSettings = useAppStore((state) => state.saveSettings);
  const history = useAppStore((state) => state.history);
  const clearHistory = useAppStore((state) => state.clearHistory);
  const [threads, setThreads] = useState(settings.defaultThreads);
  const [outputDir, setOutputDir] = useState(settings.defaultOutputDir);
  const [saved, setSaved] = useState(false);

  useEffect(() => {
    setThreads(settings.defaultThreads);
    setOutputDir(settings.defaultOutputDir);
  }, [settings]);

  function save() {
    saveSettings({
      defaultThreads: Math.max(1, Math.trunc(threads || 1)),
      defaultOutputDir: outputDir.trim(),
    });
    setSaved(true);
    window.setTimeout(() => setSaved(false), 1800);
  }

  async function chooseDefaultOutput() {
    const path = await pickDirectory("Choose the default output directory");
    if (path) setOutputDir(path);
  }

  return (
    <div className="page">
      <PageHeader
        eyebrow="Preferences"
        title="Desktop settings"
        description="Defaults are applied to new runs. Existing run records keep their original parameters."
      />

      <div className="settings-grid">
        <section className="panel">
          <div className="section-heading">
            <span className="section-number">01</span>
            <div>
              <h2>Run defaults</h2>
              <p>Choose the initial thread count and destination for new analyses.</p>
            </div>
          </div>
          <div className="settings-form">
            <div className="field">
              <label htmlFor="default-threads">
                Default threads
                <InfoHint text="Applied to new runs. Existing runs keep their saved thread count." />
              </label>
              <input
                id="default-threads"
                type="number"
                min={1}
                max={512}
                value={threads}
                onChange={(event) => setThreads(Number(event.target.value))}
              />
              <small>This system reports {runtime?.logicalCpus ?? 1} logical CPUs.</small>
            </div>
            <div className="field">
              <label htmlFor="default-output">
                Default output directory
                <InfoHint text="Leave blank to create a result folder beside the selected input FASTA." />
              </label>
              <div className="input-with-button">
                <input
                  id="default-output"
                  value={outputDir}
                  onChange={(event) => setOutputDir(event.target.value)}
                  placeholder="Leave blank to create a folder beside the input"
                />
                <button className="btn secondary" onClick={chooseDefaultOutput}>
                  <FolderOpen size={16} />
                  Browse
                </button>
              </div>
            </div>
            <div className="settings-actions">
              <button
                className="btn ghost"
                onClick={() => {
                  setThreads(runtime?.logicalCpus ?? 1);
                  setOutputDir("");
                }}
              >
                <RotateCcw size={16} />
                Reset
              </button>
              <button className="btn primary" onClick={save}>
                {saved ? <Check size={17} /> : <Save size={17} />}
                {saved ? "Saved" : "Save settings"}
              </button>
            </div>
          </div>
        </section>

        <section className="panel">
          <div className="section-heading">
            <span className="section-number">02</span>
            <div>
              <h2>Application</h2>
              <p>Runtime information and local data controls.</p>
            </div>
          </div>
          <div className="settings-info">
            <div>
              <Cpu size={18} />
              <span>Target</span>
              <strong>
                {runtime?.os} / {runtime?.arch}
              </strong>
            </div>
            <div>
              <HardDrive size={18} />
              <span>FastNLR version</span>
              <strong>{runtime?.version}</strong>
            </div>
            <div>
              <Settings size={18} />
              <span>Stored run records</span>
              <strong>{history.length} / 50</strong>
            </div>
          </div>
          <div className="danger-zone">
            <div>
              <strong>Clear local run history</strong>
              <p>This removes records only. Output folders and files are never deleted.</p>
            </div>
            <button
              className="btn secondary danger-text"
              disabled={!history.length}
              onClick={() => {
                if (window.confirm("Remove all local run records?")) clearHistory();
              }}
            >
              <Trash2 size={16} />
              Clear
            </button>
          </div>
        </section>
      </div>
    </div>
  );
}
