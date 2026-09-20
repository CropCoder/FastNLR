//! nlr-seq -- sequence processing: six-frame translation, reverse complement, codon table, FASTA reading, sequence chopping.
//!
//! Sequence-processing primitives: translation, reverse complement, codon table, FASTA reading, and chopping.

pub mod codon;
pub mod chopper;
pub mod fasta;
pub mod translate;

pub use chopper::SequenceChopper;
pub use fasta::FastaReader;
pub use translate::BioSequence;
