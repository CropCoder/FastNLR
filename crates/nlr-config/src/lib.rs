//! nlr-config — parse mot.txt (PWM) and store.txt (CDF).
//!
//! Faithful reimplementation of Java `MotifDefinition`:
//! - mot.txt: `motif@pos@aa score`, pos 0-based, aa is a single A-Z char, integer score;
//! - store.txt: `motif@score pvalue`, integer index + double p-value;
//! - PWM uses `[aa - 'A']` (ASCII 65) as the first dimension index; non A-Z chars yield 0;
//! - `motif_names` records the order motifs first appear in mot.txt (scan iteration order).

/// Motif id (1..=20), matching nlr-core's MotifId.
pub type MotifId = u8;

/// Motif definition (PWM + CDF + length + order of appearance).
pub struct MotifDefinition {
    /// `pwm[id]` = flattened PWM, `[aa(0..25) * width + pos]`, length `26 * width`.
    pwm: Vec<Vec<i32>>,
    /// `cdf[id]` = CDF array indexed by integer score.
    cdf: Vec<Vec<f64>>,
    /// `lengths[id]` = motif length (amino acids).
    lengths: Vec<u16>,
    /// Order of first appearance (mot.txt order).
    motif_names: Vec<MotifId>,
    /// Max motif length.
    max_length: u16,
}

impl MotifDefinition {
    /// Load from mot.txt and store.txt file paths (equivalent to Java constructor).
    pub fn load(pwm_file: &std::path::Path, cdf_file: &std::path::Path) -> std::io::Result<Self> {
        let mot_text = read_mmap_str(pwm_file)?;
        let store_text = read_mmap_str(cdf_file)?;
        Self::load_from_str(&mot_text, &store_text)
    }

    /// Load from in-memory mot.txt / store.txt text (enables `include_str!` embedding).
    pub fn load_from_str(mot: &str, store: &str) -> std::io::Result<Self> {
        let (lengths, motif_names) = Self::lengths_from(mot);
        let pwm = Self::pwm_from(mot, &lengths);
        let cdf = Self::cdf_from(store);
        let max_length = lengths.iter().copied().max().unwrap_or(0);
        Ok(MotifDefinition {
            pwm,
            cdf,
            lengths,
            motif_names,
            max_length,
        })
    }

    /// Parse "motif_N" -> id N.
    fn parse_motif_id(s: &str) -> Option<MotifId> {
        s.strip_prefix("motif_")?.parse::<MotifId>().ok()
    }

    /// First pass: compute each motif's length (equivalent to `loadMotifLengths`).
    fn lengths_from(text: &str) -> (Vec<u16>, Vec<MotifId>) {
        let mut lengths = vec![0u16; 21];
        let mut motif_names: Vec<MotifId> = Vec::new();
        let mut seen = [false; 21];

        for line in text.lines() {
            let first = match line.split_whitespace().next() {
                Some(f) => f,
                None => continue,
            };
            let parts: Vec<&str> = first.split('@').collect();
            if parts.len() < 2 {
                continue;
            }
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => continue,
            };
            // mot.txt position is 0-based; length = max pos + 1.
            let pos = parts[1].parse::<usize>().unwrap_or(0) + 1;
            if !seen[id as usize] {
                seen[id as usize] = true;
                motif_names.push(id);
            }
            if pos as u16 > lengths[id as usize] {
                lengths[id as usize] = pos as u16;
            }
        }
        (lengths, motif_names)
    }

    /// Second pass: fill PWM (equivalent to `loadPWM`).
    fn pwm_from(text: &str, lengths: &[u16]) -> Vec<Vec<i32>> {
        let mut pwm: Vec<Vec<i32>> = (0..21)
            .map(|i| vec![0; 26 * lengths[i] as usize])
            .collect();

        for line in text.lines() {
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            let parts: Vec<&str> = cols[0].split('@').collect();
            if parts.len() < 3 {
                continue;
            }
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => continue,
            };
            let position = parts[1].parse::<usize>().unwrap_or(0);
            let aa = parts[2].as_bytes().first().copied().unwrap_or(b'A');
            let score = cols[1].parse::<i32>().unwrap_or(0);

            let width = lengths[id as usize] as usize;
            let aa_idx = (aa as usize).saturating_sub(b'A' as usize);
            if aa_idx < 26 && position < width {
                pwm[id as usize][aa_idx * width + position] = score;
            }
        }
        pwm
    }

    /// Parse store.txt (CDF), two passes (equivalent to `loadCDF`).
    fn cdf_from(text: &str) -> Vec<Vec<f64>> {
        // First pass: find max index per motif.
        let mut max_idx = vec![0usize; 21];
        for line in text.lines() {
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.is_empty() {
                continue;
            }
            let parts: Vec<&str> = cols[0].split('@').collect();
            if parts.len() < 2 {
                continue;
            }
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => continue,
            };
            let idx = parts[1].parse::<usize>().unwrap_or(0);
            if idx > max_idx[id as usize] {
                max_idx[id as usize] = idx;
            }
        }

        // Second pass: fill.
        let mut cdf: Vec<Vec<f64>> = (0..21).map(|i| vec![0.0; max_idx[i] + 1]).collect();
        for line in text.lines() {
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < 2 {
                continue;
            }
            let parts: Vec<&str> = cols[0].split('@').collect();
            if parts.len() < 2 {
                continue;
            }
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => continue,
            };
            let idx = parts[1].parse::<usize>().unwrap_or(0);
            let val = cols[1].parse::<f64>().unwrap_or(0.0);
            if idx < cdf[id as usize].len() {
                cdf[id as usize][idx] = val;
            }
        }
        cdf
    }

    /// Order of first appearance (scan iteration order).
    #[inline]
    pub fn motif_names(&self) -> &[MotifId] {
        &self.motif_names
    }

    /// Motif length (amino acids).
    #[inline]
    pub fn length(&self, id: MotifId) -> u16 {
        self.lengths[id as usize]
    }

    /// Max motif length (sequences shorter than this skip scanning).
    #[inline]
    pub fn max_length(&self) -> u16 {
        self.max_length
    }

    /// Query PWM score (equivalent to `getScore`): aa < 'A' returns 0.
    #[inline]
    pub fn score(&self, id: MotifId, position: usize, aa: u8) -> i32 {
        // Java: `if ((int) aa - 65 < 0) return 0;` i.e. aa < 'A' returns 0.
        if aa < b'A' {
            return 0;
        }
        let aa_idx = (aa - b'A') as usize;
        if aa_idx >= 26 {
            return 0;
        }
        let width = self.lengths[id as usize] as usize;
        if position >= width {
            return 0;
        }
        self.pwm[id as usize][aa_idx * width + position]
    }

    /// Query CDF (equivalent to `getCDF`). CDF is right-tail cumulative (decreasing).
    #[inline]
    pub fn cdf(&self, id: MotifId, score: i32) -> f64 {
        let idx = score.max(0) as usize;
        let arr = &self.cdf[id as usize];
        if idx < arr.len() {
            arr[idx]
        } else {
            // Score beyond table upper bound -> highly significant (right-tail p ~= 0).
            // Java throws out-of-bounds here; returning 0 is equivalent to highly significant.
            0.0
        }
    }

    /// Precompute integer thresholds where `p < thresh <=> score >= T(id)`.
    ///
    /// Note: CDF is **right-tail cumulative** (`p = P(random score >= score)`),
    /// monotonically decreasing — higher score is more significant.
    /// Returns `thresholds[id]`: when score >= thresholds[id], pvalue < thresh.
    pub fn score_thresholds(&self, thresh: f64) -> Vec<i32> {
        let mut thresholds = vec![0i32; 21];
        for id in 1..=20u8 {
            let arr = &self.cdf[id as usize];
            // Find first score where cdf[score] < thresh (CDF is decreasing).
            let mut t = arr.len() as i32;
            for (score, &v) in arr.iter().enumerate() {
                if v < thresh {
                    t = score as i32;
                    break;
                }
            }
            thresholds[id as usize] = t;
        }
        thresholds
    }
}

/// Load from standard filenames (mot.txt / store.txt) in a directory.
pub fn load_default(dir: &std::path::Path) -> std::io::Result<MotifDefinition> {
    MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt"))
}

/// Memory-map a file read-only and return its contents as a string.
fn read_mmap_str(path: &std::path::Path) -> std::io::Result<String> {
    let file = std::fs::File::open(path)?;
    // SAFETY: read-only mapping; file lifetime held by Mmap.
    let mmap = unsafe { memmap2::Mmap::map(&file) }?;
    let bytes: &[u8] = &mmap;
    std::str::from_utf8(bytes)
        .map(|s| s.to_string())
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// Parse an exported TSV motif line into (motif_id, protein sequence id, position, pvalue).
/// Reused by `-c` import.
pub fn parse_motif_line(line: &str) -> Option<(MotifId, String, u64, f64)> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 6 {
        return None;
    }
    let id = cols[0]
        .trim_start_matches("motif_")
        .parse::<MotifId>()
        .ok()?;
    let position = cols[2].parse::<u64>().ok()?;
    let pvalue = cols[4].parse::<f64>().ok()?;
    Some((id, cols[1].to_string(), position, pvalue))
}
