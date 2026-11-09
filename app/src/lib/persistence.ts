/**
 * Persistent desktop settings and run history storage.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { LazyStore } from "@tauri-apps/plugin-store";
import type { PersistedAppState } from "../types";

const STORE_FILE = "fastnlr-state.json";
const STORE_KEY = "state";
const LOCAL_KEY = "fastnlr-state";

let store: LazyStore | null = null;

function isTauriRuntime() {
  return typeof window !== "undefined" && "__TAURI_INTERNALS__" in window;
}

function getStore() {
  store ??= new LazyStore(STORE_FILE);
  return store;
}

export async function loadPersistedState(): Promise<PersistedAppState | null> {
  try {
    if (isTauriRuntime()) {
      return (await getStore().get<PersistedAppState>(STORE_KEY)) ?? null;
    }
    const value = window.localStorage.getItem(LOCAL_KEY);
    return value ? (JSON.parse(value) as PersistedAppState) : null;
  } catch {
    return null;
  }
}

export async function savePersistedState(state: PersistedAppState): Promise<void> {
  try {
    if (isTauriRuntime()) {
      await getStore().set(STORE_KEY, state);
      await getStore().save();
      return;
    }
    window.localStorage.setItem(LOCAL_KEY, JSON.stringify(state));
  } catch {
    // 持久化失败不应阻断正在进行的注释任务。
  }
}
