//! nlr-core — FastNLR domain model.
//!
//! This crate has no IO/external dependencies and only contains:
//! - Basic enums such as strand and domain category;
//! - Hard-coded static rule tables (rank, category, color, seed combinations, signatures, consensus sequences);
//! - Core types and business rules for Motif / MotifList (coordinate mapping, sorting, merging, completeness checks).

pub mod motif;
pub mod motif_list;
pub mod signature_def;
pub mod strand;

pub use motif::Motif;
pub use motif_list::MotifList;
pub use signature_def::AnnotatorSignatureDefinition;
pub use strand::Strand;
