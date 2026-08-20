//! nlr-seq -- sequence processing: six-frame translation, reverse complement, codon table, FASTA reading, sequence chopping.
//!
//! Faithful reimplementation of Java `BioSequence` / `FastaReader` / `SequenceChopper`.

pub mod codon;
pub mod chopper;
pub mod fasta;
pub mod translate;

pub use chopper::SequenceChopper;
pub use fasta::FastaReader;
pub use translate::BioSequence;
