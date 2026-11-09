/**
 * Zustand application state for jobs, settings, and local history.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { create } from "zustand";
import { getRuntimeInfo } from "../lib/bridge";
import { loadPersistedState, savePersistedState } from "../lib/persistence";
import type {
  ActiveRun,
  AppSettings,
  RunDraft,
  RunHistoryEntry,
  RuntimeInfo,
} from "../types";

const MAX_HISTORY = 50;

interface AppStore {
  initialized: boolean;
  runtime: RuntimeInfo | null;
  settings: AppSettings;
  history: RunHistoryEntry[];
  activeRun: ActiveRun | null;
  pendingDraft: RunDraft | null;
  initialize: () => Promise<void>;
  saveSettings: (settings: AppSettings) => void;
  addHistory: (entry: RunHistoryEntry) => void;
  updateHistory: (id: string, patch: Partial<RunHistoryEntry>) => void;
  clearHistory: () => void;
  setActiveRun: (run: ActiveRun | null) => void;
  setPendingDraft: (draft: RunDraft | null) => void;
}

const initialSettings: AppSettings = {
  defaultThreads: 1,
  defaultOutputDir: "",
};

let persistTimer: number | null = null;
let pendingPersist: { settings: AppSettings; history: RunHistoryEntry[] } | null =
  null;

function persistFromStore(get: () => AppStore, immediate = false) {
  const { settings, history } = get();
  pendingPersist = { settings, history };
  if (immediate) {
    if (persistTimer) window.clearTimeout(persistTimer);
    persistTimer = null;
    const payload = pendingPersist;
    pendingPersist = null;
    void savePersistedState(payload);
    return;
  }
  if (persistTimer) return;
  persistTimer = window.setTimeout(() => {
    persistTimer = null;
    const payload = pendingPersist;
    pendingPersist = null;
    if (payload) void savePersistedState(payload);
  }, 2000);
}

export const useAppStore = create<AppStore>((set, get) => ({
  initialized: false,
  runtime: null,
  settings: initialSettings,
  history: [],
  activeRun: null,
  pendingDraft: null,

  initialize: async () => {
    if (get().initialized) return;
    const [runtime, persisted] = await Promise.all([
      getRuntimeInfo(),
      loadPersistedState(),
    ]);
    const interruptedAt = new Date().toISOString();
    let historyChanged = false;
    const history = (persisted?.history ?? []).map((entry) => {
      if (entry.status !== "Running" && entry.status !== "Pending") return entry;
      historyChanged = true;
      return {
        ...entry,
        status: "Interrupted" as const,
        phase: "Interrupted" as const,
        finishedAt: entry.finishedAt ?? interruptedAt,
        error: "The application closed before this run completed.",
      };
    });
    const fallbackThreads = persisted?.settings.defaultThreads || runtime.logicalCpus || 1;
    const settings = {
      defaultThreads: fallbackThreads,
      defaultOutputDir: persisted?.settings.defaultOutputDir ?? "",
    };
    set({ runtime, settings, history, initialized: true });
    if (
      persisted == null ||
      historyChanged ||
      persisted.settings.defaultThreads !== settings.defaultThreads ||
      persisted.settings.defaultOutputDir !== settings.defaultOutputDir
    ) {
      void savePersistedState({ settings, history });
    }
  },

  saveSettings: (settings) => {
    set({ settings });
    persistFromStore(get, true);
  },

  addHistory: (entry) => {
    set((state) => ({
      history: [entry, ...state.history].slice(0, MAX_HISTORY),
    }));
    persistFromStore(get, true);
  },

  updateHistory: (id, patch) => {
    set((state) => ({
      history: state.history.map((entry) =>
        entry.id === id ? { ...entry, ...patch } : entry,
      ),
    }));
    persistFromStore(get);
  },

  clearHistory: () => {
    set({ history: [] });
    persistFromStore(get, true);
  },

  setActiveRun: (activeRun) => set({ activeRun }),
  setPendingDraft: (pendingDraft) => set({ pendingDraft }),
}));
