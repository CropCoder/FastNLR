/**
 * Application version, citation, and project information.
 * Author: Jiwen Zhao (https://github.com/CropCoder)
 */
import {
  BookOpen,
  Dna,
  ExternalLink,
  Github,
  Scale,
  ShieldCheck,
} from "lucide-react";
import { PageHeader } from "../components/PageHeader";
import { useAppStore } from "../store/appStore";

export function AboutPage() {
  const runtime = useAppStore((state) => state.runtime);

  return (
    <div className="page">
      <PageHeader
        eyebrow="About FastNLR"
        title="Fast, self-contained NLR annotation"
        description="A native desktop interface for the FastNLR Rust pipeline, built for plant genome and protein analyses."
      />

      <section className="about-hero panel">
        <div className="about-mark">
          <Dna size={42} />
        </div>
        <div>
          <span className="eyebrow">Desktop application</span>
          <h2>FastNLR v{runtime?.version}</h2>
          <p>
            Six-frame genome scanning, motif assembly, HMM-assisted domain evidence,
            and publication-ready output formats in one local application.
          </p>
        </div>
      </section>

      <div className="about-grid">
        <article className="panel about-card">
          <BookOpen size={22} />
          <h3>Citation</h3>
          <p>
            A dedicated manuscript is in preparation. Until then, cite the FastNLR
            repository and the software version used for analysis.
          </p>
        </article>
        <article className="panel about-card">
          <GitHubLink />
          <h3>Source and releases</h3>
          <p>Download the latest builds, report issues, and inspect the implementation.</p>
          <a className="text-link" href="https://github.com/CropCoder/FastNLR" target="_blank" rel="noreferrer">
            github.com/CropCoder/FastNLR
            <ExternalLink size={14} />
          </a>
        </article>
        <article className="panel about-card">
          <Scale size={22} />
          <h3>License</h3>
          <p>FastNLR and this desktop application are distributed under GPL-3.0-only.</p>
        </article>
        <article className="panel about-card">
          <ShieldCheck size={22} />
          <h3>Local processing</h3>
          <p>
            Input sequences and results remain on this computer. The application does
            not transmit analysis data or usage telemetry.
          </p>
        </article>
      </div>

      <div className="about-footer">
        Created by Jiwen Zhao · CropCoder
      </div>
    </div>
  );
}

function GitHubLink() {
  return <Github size={22} />;
}
