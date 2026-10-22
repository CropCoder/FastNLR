//! nlr-config — parse mot.txt (PWM) and store.txt (CDF).
//!
//! Motif definition parsing:
//! - mot.txt: `motif@pos@aa score`, pos 0-based, aa is a single A-Z char, integer score;
//! - store.txt: `motif@score pvalue`, integer index + double p-value;
//! - PWM uses `[aa - 'A']` (ASCII 65) as the first dimension index; non A-Z chars yield 0;
//! - `motif_names` records the order motifs first appear in mot.txt (scan iteration order).

/// Motif id (1-based), matching nlr-core's MotifId.
/// 上限由 mot.txt 里出现的最大 id 决定（不再硬编码 20）——这样外部 mot.txt/store.txt
/// 可以带更多 motif（例如从特定基因家族重挖出来的 CC/TIR motif）。
pub type MotifId = u8;

/// Motif definition (PWM + CDF + length + order of appearance).
#[derive(Debug)]
pub struct MotifDefinition {
    /// mot.txt 里出现的最大 motif id（数组尺寸 = 该值 + 1）。
    max_motif_id: MotifId,
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
    /// mot.txt 里最大的 motif id（外部库可能 > 20，数组尺寸 = 该值 + 1）。
    #[inline]
    pub fn max_motif_id(&self) -> MotifId {
        self.max_motif_id
    }

    /// Whether an id has a parsed PWM in this profile.
    #[inline]
    pub fn has_motif(&self, id: MotifId) -> bool {
        self.lengths.get(id as usize).copied().unwrap_or(0) > 0
    }

    /// Load from `mot.txt` and `store.txt` file paths.
    pub fn load(pwm_file: &std::path::Path, cdf_file: &std::path::Path) -> std::io::Result<Self> {
        let mot_text = read_mmap_str(pwm_file)?;
        let store_text = read_mmap_str(cdf_file)?;
        Self::load_from_str(&mot_text, &store_text)
    }

    /// Load from in-memory mot.txt / store.txt text (enables `include_str!` embedding).
    pub fn load_from_str(mot: &str, store: &str) -> std::io::Result<Self> {
        if mot.trim().is_empty() {
            return Err(invalid_data("mot.txt is empty"));
        }
        if store.trim().is_empty() {
            return Err(invalid_data("store.txt is empty"));
        }
        let (lengths, motif_names) = Self::lengths_from(mot)?;
        if motif_names.is_empty() {
            return Err(invalid_data("mot.txt contains no motif entries"));
        }
        let slots = lengths.len();
        let max_motif_id = (slots.saturating_sub(1)) as MotifId;
        let pwm = Self::pwm_from(mot, &lengths)?;
        let cdf = Self::cdf_from(store)?;
        let max_length = lengths.iter().copied().max().unwrap_or(0);

        // 一致性校验：mot.txt 与 store.txt 的最大 id 各自独立推导，二者不一致时
        // `score_thresholds()`/`cdf()` 会按 mot.txt 的 id 索引 store.txt 的表而越界 panic。
        // 这里提前拦截并给出明确的 motif 列表，而不是等到查表时崩溃。
        let missing: Vec<String> = (1..=max_motif_id)
            .filter(|&id| lengths[id as usize] > 0)
            .filter(|&id| match cdf.get(id as usize) {
                Some(table) => table.is_empty(),
                None => true,
            })
            .map(|id| format!("motif_{}", id))
            .collect();
        if !missing.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "store.txt has no CDF table for {}: mot.txt declares up to motif_{} but store.txt only covers up to motif_{}",
                    missing.join(", "),
                    max_motif_id,
                    cdf.len().saturating_sub(1),
                ),
            ));
        }

        Ok(MotifDefinition {
            pwm,
            cdf,
            lengths,
            motif_names,
            max_length,
            max_motif_id,
        })
    }

    /// Parse "motif_N" -> id N.
    fn parse_motif_id(s: &str) -> Option<MotifId> {
        s.strip_prefix("motif_")?.parse::<MotifId>().ok()
    }

    /// First pass: compute each motif's length.
    fn lengths_from(text: &str) -> std::io::Result<(Vec<u16>, Vec<MotifId>)> {
        // 第一遍：扫出最大 motif id（决定数组尺寸，支持 >20 的外部库）
        let mut max_id = 0usize;
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() != 2 {
                return Err(invalid_data(format!(
                    "invalid mot.txt line, expected 2 columns: {line}"
                )));
            }
            let parts: Vec<&str> = cols[0].split('@').collect();
            if parts.len() != 3 {
                return Err(invalid_data(format!(
                    "invalid mot.txt key, expected motif_ID@position@AA: {line}"
                )));
            }
            let id = Self::parse_motif_id(parts[0])
                .filter(|&id| id > 0)
                .ok_or_else(|| invalid_data(format!("invalid motif id in line: {line}")))?;
            parts[1]
                .parse::<usize>()
                .map_err(|_| invalid_data(format!("invalid motif position in line: {line}")))?;
            let aa = parts[2].as_bytes();
            if aa.len() != 1 || !aa[0].is_ascii_uppercase() {
                return Err(invalid_data(format!(
                    "invalid amino acid in mot.txt line: {line}"
                )));
            }
            cols[1]
                .parse::<i32>()
                .map_err(|_| invalid_data(format!("invalid PWM score in line: {line}")))?;
            max_id = max_id.max(id as usize);
        }
        let slots = max_id + 1;
        let mut lengths = vec![0u16; slots];
        let mut motif_names: Vec<MotifId> = Vec::new();
        let mut seen = vec![false; slots];

        let mut seen_keys = std::collections::HashSet::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let first = match line.split_whitespace().next() {
                Some(f) => f,
                None => continue,
            };
            let parts: Vec<&str> = first.split('@').collect();
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => return Err(invalid_data(format!("invalid motif id in line: {line}"))),
            };
            // mot.txt position is 0-based; length = max pos + 1.
            let pos = parts[1]
                .parse::<usize>()
                .map_err(|_| invalid_data(format!("invalid motif position in line: {line}")))?
                + 1;
            if !seen_keys.insert((id, pos - 1, parts[2].to_string())) {
                return Err(invalid_data(format!(
                    "duplicate PWM entry in mot.txt: {line}"
                )));
            }
            if !seen[id as usize] {
                seen[id as usize] = true;
                motif_names.push(id);
            }
            if pos as u16 > lengths[id as usize] {
                lengths[id as usize] = pos as u16;
            }
        }
        Ok((lengths, motif_names))
    }

    /// Second pass: fill the position-weight matrix.
    fn pwm_from(text: &str, lengths: &[u16]) -> std::io::Result<Vec<Vec<i32>>> {
        let mut pwm: Vec<Vec<i32>> = (0..lengths.len())
            .map(|i| vec![0; 26 * lengths[i] as usize])
            .collect();

        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            let parts: Vec<&str> = cols[0].split('@').collect();
            let id = match Self::parse_motif_id(parts[0]) {
                Some(i) => i,
                None => return Err(invalid_data(format!("invalid motif id in line: {line}"))),
            };
            let position = parts[1]
                .parse::<usize>()
                .map_err(|_| invalid_data(format!("invalid motif position in line: {line}")))?;
            let aa = parts[2]
                .as_bytes()
                .first()
                .copied()
                .ok_or_else(|| invalid_data(format!("empty amino acid in line: {line}")))?;
            let score = cols[1]
                .parse::<i32>()
                .map_err(|_| invalid_data(format!("invalid PWM score in line: {line}")))?;

            let width = lengths[id as usize] as usize;
            let aa_idx = (aa as usize).saturating_sub(b'A' as usize);
            if aa_idx < 26 && position < width {
                pwm[id as usize][aa_idx * width + position] = score;
            }
        }
        Ok(pwm)
    }

    /// Parse and validate `store.txt` (right-tail CDF).
    fn cdf_from(text: &str) -> std::io::Result<Vec<Vec<f64>>> {
        let mut entries: Vec<std::collections::HashMap<usize, f64>> = vec![Default::default()];
        let mut seen = std::collections::HashSet::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() != 2 {
                return Err(invalid_data(format!(
                    "invalid store.txt line, expected 2 columns: {line}"
                )));
            }
            let parts: Vec<&str> = cols[0].split('@').collect();
            if parts.len() != 2 {
                return Err(invalid_data(format!(
                    "invalid store.txt key, expected motif_ID@score: {line}"
                )));
            }
            let id = Self::parse_motif_id(parts[0])
                .filter(|&id| id > 0)
                .ok_or_else(|| invalid_data(format!("invalid motif id in line: {line}")))?;
            let idx = parts[1]
                .parse::<usize>()
                .map_err(|_| invalid_data(format!("invalid CDF score index in line: {line}")))?;
            let val = cols[1]
                .parse::<f64>()
                .map_err(|_| invalid_data(format!("invalid CDF p-value in line: {line}")))?;
            if !val.is_finite() || !(0.0..=1.0).contains(&val) {
                return Err(invalid_data(format!(
                    "CDF p-value must be finite and within [0,1]: {line}"
                )));
            }
            if !seen.insert((id, idx)) {
                return Err(invalid_data(format!(
                    "duplicate CDF entry in store.txt: {line}"
                )));
            }
            if entries.len() <= id as usize {
                entries.resize_with(id as usize + 1, Default::default);
            }
            entries[id as usize].insert(idx, val);
        }
        if seen.is_empty() {
            return Err(invalid_data("store.txt contains no CDF entries"));
        }

        let mut cdf = Vec::with_capacity(entries.len());
        for table in entries {
            let max_idx = table.keys().copied().max().unwrap_or(0);
            let mut values = vec![1.0; max_idx + 1];
            for (&idx, &val) in &table {
                values[idx] = val;
            }
            if values.windows(2).any(|w| w[0] < w[1]) {
                return Err(invalid_data(
                    "CDF must be monotonically non-increasing with score",
                ));
            }
            cdf.push(values);
        }
        Ok(cdf)
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

    /// Maximum motif length across the profile.
    #[inline]
    pub fn max_length(&self) -> u16 {
        self.max_length
    }

    /// Query the PWM score; characters below `A` return 0.
    #[inline]
    pub fn score(&self, id: MotifId, position: usize, aa: u8) -> i32 {
        // Characters below 'A' return 0.
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

    /// Query the CDF. The CDF is right-tail cumulative (decreasing).
    #[inline]
    pub fn cdf(&self, id: MotifId, score: i32) -> f64 {
        let idx = score.max(0) as usize;
        let arr = &self.cdf[id as usize];
        // 分数超出 CDF 表覆盖范围时，真实 p 值只会比表内最后一个更小，因此用表内最后一个
        // （最小）p 值作为保守上界。旧实现直接返回 0.0，等于把"未知"当成"极显著"，
        // 一旦自定义 store.txt 的分数区间不足就会静默产生大量假阳性。
        arr.get(idx)
            .copied()
            .unwrap_or_else(|| arr.last().copied().unwrap_or(1.0))
    }

    /// Precompute integer thresholds where `p < thresh <=> score >= T(id)`.
    ///
    /// Note: CDF is **right-tail cumulative** (`p = P(random score >= score)`),
    /// monotonically decreasing — higher score is more significant.
    /// Returns `thresholds[id]`: when score >= thresholds[id], pvalue < thresh.
    pub fn score_thresholds(&self, thresh: f64) -> Vec<i32> {
        let mut thresholds = vec![0i32; self.cdf.len()];
        for id in 1..=self.max_motif_id {
            let arr = &self.cdf[id as usize];
            // 找首个满足 cdf[score] < thresh 的分数（CDF 递减）。若表内没有任何分数满足，
            // 说明该 motif 不可能达到该显著性，用 i32::MAX 表示"没有窗口能通过预筛"；
            // 这与 cdf() 对越界分数返回表内最小 p 值的语义一致（也 >= thresh，会被拒）。
            let mut t = i32::MAX;
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

fn invalid_data(message: impl Into<String>) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message.into())
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
