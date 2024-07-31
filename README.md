# CigarParser
A Rust library to parse Cigar string and, in particular, identify junction position. 
Cigar string are a format to summarize alignment use in bioformatic present in BAM and SAM format.
They are used by NGS aligner to summarize how a read align to a genome.

        how to use
        you have two options to create a Cigar object:
        # This return a Results allowing fine grained control
        let cig = Cigar::from_str("35M110N45M3I45M10N11M");
        # Thiw will panics if the Cigar string is not properly formated
        let cig = Cigar::from("35M110N45M3I45M10N11M");
        
        # to get the juction if any use:
        let cig = Cigar::from("35M110N45M3I45M10N");
        let results = cig.get_skipped_pos_on_ref(&500);
        assert_eq!(results, Some(vec![535, 645, 738, 748]))

I wrote it at a standalaone so you can integrate it with any tools that reads bam file such as rust-htslib.

Suggestion and comment are welcome!
