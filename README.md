# Cigar String Parser

A Rust library to parse CIGAR strings and, in particular, identify junction positions. CIGAR strings are a format used to summarize alignments in bioinformatics, present in BAM and SAM formats. They are used by NGS aligners to summarize how a read aligns to a genome.

## How to Use

in your toml add this dependancy
```rust
CigarParser = { git = "https://github.com/rLannes/CigarParser" }
```

This crate exposes the Cigar structure that you can load using the following:

```rust
extern crate CigarParser;
use CigarParser::cigar::Cigar;
```

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
let results = cig.get_junction_position(&500);// 500 is the starting position of the alignment
assert_eq!(results, Some(vec![535, 645, 738, 748]));
```

To Compute the coverage one can parse each read instead of using the pileup method.
```rust
let cig = Cigar::from("5M15N5M");
let results = cig.get_reference_cover(&500); // 500 is the starting position of the alignment
assert_eq!(results, Some(vec![500, 501, 502, 503, 504, 520, 521, 522, 523, 524]));
```

To check if the read overlap fully an interval: 
```rust
let cig = Cigar::from("5M15N5M");
let results = cig.does_it_match_an_intervall(&500, 501, 503);// 500 is the starting position of the alignment
assert_eq!(results, true);
```

To get the end of the alignment : 
```rust
let cig = Cigar::from("5M15N5M");
let results = cig.get_end_of_aln(&500);// 500 is the starting position of the alignment
assert_eq!(results, 519);
```

get_end_of_aln


 wrote this as a standalone library so you can integrate it with any tools that read BAM files, such as rust-htslib.
Suggestions and comments are welcome!
