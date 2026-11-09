/**
 * Desktop route shell and lazy page loading.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import { lazy, Suspense, useEffect } from "react";
import { Navigate, Route, Routes } from "react-router-dom";
import { AppShell } from "./components/AppShell";
import { useAppStore } from "./store/appStore";

const NewRunPage = lazy(() =>
  import("./pages/NewRunPage").then((module) => ({ default: module.NewRunPage })),
);
const RunMonitorPage = lazy(() =>
  import("./pages/RunMonitorPage").then((module) => ({
    default: module.RunMonitorPage,
  })),
);
const ResultsPage = lazy(() =>
  import("./pages/ResultsPage").then((module) => ({ default: module.ResultsPage })),
);
const HistoryPage = lazy(() =>
  import("./pages/HistoryPage").then((module) => ({ default: module.HistoryPage })),
);
const SettingsPage = lazy(() =>
  import("./pages/SettingsPage").then((module) => ({ default: module.SettingsPage })),
);
const AboutPage = lazy(() =>
  import("./pages/AboutPage").then((module) => ({ default: module.AboutPage })),
);

export default function App() {
  const initialized = useAppStore((state) => state.initialized);
  const initialize = useAppStore((state) => state.initialize);

  useEffect(() => {
    void initialize();
  }, [initialize]);

  if (!initialized) {
    return (
      <div className="boot-screen">
        <div className="brand-mark large">
          <span>F</span>
        </div>
        <div>
          <strong>FastNLR</strong>
          <p>Preparing the annotation workspace...</p>
        </div>
      </div>
    );
  }

  return (
    <Suspense
      fallback={
        <div className="boot-screen compact">
          <div className="brand-mark">
            <span>F</span>
          </div>
          <p>Loading workspace...</p>
        </div>
      }
    >
      <Routes>
        <Route element={<AppShell />}>
          <Route index element={<Navigate to="/new" replace />} />
          <Route path="/new" element={<NewRunPage />} />
          <Route path="/run" element={<RunMonitorPage />} />
          <Route path="/results/:historyId" element={<ResultsPage />} />
          <Route path="/history" element={<HistoryPage />} />
          <Route path="/settings" element={<SettingsPage />} />
          <Route path="/about" element={<AboutPage />} />
          <Route path="*" element={<Navigate to="/new" replace />} />
        </Route>
      </Routes>
    </Suspense>
  );
}
