# Cigar String Parser

A Rust library to parse CIGAR strings and, in particular, identify junction positions. CIGAR strings are a format used to summarize alignments in bioinformatics, present in BAM and SAM formats. They are used by NGS aligners to summarize how a read aligns to a genome.

## How to Use

You have two options to create a Cigar object:

```rust
// This returns a Result allowing fine-grained control
let cig = Cigar::from_str("35M110N45M3I45M10N11M");

// This will panic if the CIGAR string is not properly formatted
let cig = Cigar::from("35M110N45M3I45M10N11M");
```
To get the junction positions, if any, use:.
```rust
let cig = Cigar::from("35M110N45M3I45M10N");
let results = cig.get_skipped_pos_on_ref(500);
assert_eq!(results, Some(vec![535, 645, 738, 748]));
```
I wrote this as a standalone library so you can integrate it with any tools that read BAM files, such as rust-htslib.
Suggestions and comments are welcome!
