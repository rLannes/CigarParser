

use std::fmt;

pub mod cigar{
    

    //! # CIGAR String Parsing and Manipulation
    //! 
    //! This module provides functionality for working with CIGAR strings,
    //! which are used in bioinformatics to represent sequence alignments.
    //! 
    //! A CIGAR string consists of a series of length + operation pairs that describe
    //! how a sequence aligns to a reference. For example, "35M110N45M" means:
    //! - 35 bases match the reference
    //! - 110 bases are skipped in the reference (e.g., an intron)
    //! - 45 more bases match the reference
    //! 
    //! - Parsing CIGAR strings
    //! - Analyzing alignment properties
    //! - Detecting splicing junctions
    //! - Working with soft/hard clipped sequences
    //!
    //! ## Example
    //! 
    //! ```rust
    //! use CigarParser::cigar::Cigar;
    //! 
    //! let cigar = Cigar::from_str("35M110N45M").unwrap();
    //! 
    //! // Check if the alignment spans an intron
    //! assert_eq!(cigar.has_skipped(), true);
    //! 
    //! // Get the positions of splicing junctions (given alignment start at position 100)
    //! let junctions = cigar.get_junction_position(100);
    //! assert_eq!(junctions, Some(vec![135, 245]));
    //! ```


    #![allow(dead_code)]
    use strand_specifier_lib::{Strand};
    use std::path::Display;
    use std::str::FromStr;
    use std::fmt;
    #[derive(Debug, PartialEq)]
    /// Basic Cigar Operation, does not accept "X" or "=".
    pub enum CigarOperation{
        Nskipped(i64),
        Match(i64),
        Insertion(i64),
        Deletion(i64),
        Soft(i64),
        Hard(i64),
        Padded(i64),
        Unaligned,
        Invalid
    }

    impl CigarOperation{
        /// Returns whether this CIGAR operation consumes positions in the reference sequence.
        ///
        /// Operations that consume reference: M (match), D (deletion), N (skipped)
        /// Operations that don't: I (insertion), S (soft clip), H (hard clip), P (padding)
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::{CigarOperation};
        /// 
        /// let match_op = CigarOperation::Match(50);
        /// assert_eq!(match_op.consume_ref(), true);
        /// 
        /// let insertion = CigarOperation::Insertion(10);
        /// assert_eq!(insertion.consume_ref(), false);
        /// ```
        pub fn consume_ref(&self) -> bool{
            match self{
                CigarOperation::Nskipped(_) => true,
                CigarOperation::Match(_)  => true,
                CigarOperation::Insertion(_)  => false,
                CigarOperation::Deletion(_)  => true,
                CigarOperation::Soft(_)  => false,
                CigarOperation::Hard(_)  => false,
                CigarOperation::Padded(_)  => false,
                _ => false
            }
        }
        
        /// Returns whether this CIGAR operation consumes positions in the query sequence.
        ///
        /// Operations that consume query: M (match), I (insertion), S (soft clip)
        /// Operations that don't: D (deletion), N (skipped), H (hard clip), P (padding)
        ///
        /// # Examples
        /// ``
        /// use CigarParser::cigar::{CigarOperation};
        /// 
        /// let match_op = CigarOperation::Match(50);
        /// assert_eq!(match_op.consume_que(), true);
        /// 
        /// let deletion = CigarOperation::Deletion(5);
        /// assert_eq!(deletion.consume_que(), false);
        /// ```
        pub fn consume_que(&self) -> bool{


            match self{
                CigarOperation::Nskipped(_)  => false,
                CigarOperation::Match(_)  => true,
                CigarOperation::Insertion(_)  => true,
                CigarOperation::Deletion(_)  => false,
                CigarOperation::Soft(_)  => true,
                CigarOperation::Hard(_)  => false,
                CigarOperation::Padded(_)  => false,
                _ => false
            }
        }
    }


    #[derive(Debug, thiserror::Error)]
    pub enum CigarError{
        #[error("Error While parsing Cigar String")]
        ParseCigarError,
    }
    

    #[derive(Debug, PartialEq)]
    /// Representation of Cigar Operation 
    /// This is the main structure users interact with.
    ///
    /// Note: This does not check the logic of operation, for example a cigar string starting or ending by N should not be possible,
    /// at the moment this code does not check for it.
    ///
    /// Warning: from(&str) can panic! use from_str(&str) for a Result<> 
    /// 
    /// let cig = Cigar::from_str("35M110N45M3I45M10N");
    /// assert_eq!(cig.has_skipped(), true);
    /// let cig = Cigar::from_str("35M45M3I45M");
    /// assert_ne!(cig.has_skipped(), true);
    /// 
    /// let cig = Cigar::from_str("35M110N45M3I45M10N");
    /// assert_eq!(cig.get_junction_position(&500), Some([535, 645, 738, 748]));
    pub struct Cigar{
        cigar: Vec<CigarOperation>,
    }


    /// Create a new Cigar struct from a &str. the &str must be a valid cigar string without "X" or "=" operation
    /// Will return an error if the cigar string is not valid.
    impl FromStr for Cigar {

        type Err = CigarError;
        fn from_str(str: &str) -> Result<Self, CigarError> {
            let mut operations = Vec::new();
            let mut length = 0 as i64;

            for c in str.chars() {
                if c.is_ascii_digit() {
                    length = length * 10 + c.to_digit(10).unwrap() as i64;
                } else {
                    let op = match c {
                        'M' => CigarOperation::Match(length),
                        'I' => CigarOperation::Insertion(length),
                        'D' => CigarOperation::Deletion(length),
                        'N' => CigarOperation::Nskipped(length),
                        'S' => CigarOperation::Soft(length),
                        'H' => CigarOperation::Hard(length),
                        'P' => CigarOperation::Padded(length),
                        '*' => CigarOperation::Unaligned,
                        _ => CigarOperation::Invalid,
                    };
                    if op == CigarOperation::Invalid{
                        return Err(CigarError::ParseCigarError);
                    }
                    if op == CigarOperation::Unaligned{
                        return Ok(Cigar{cigar:Vec::new()});
                    }
                    operations.push(op);
                    length = 0;
                }
            }
                Ok(Cigar{
                cigar: operations
                }
            )
        }
        
    }
    
    



    

    //#[deprecated(note = "use `Cigar::from_str` which returns a Result instead")]
    /// Create a new Cigar struct from a &str. the &str must be a valid cigar string without "X" or "=" operation
    /// Will Panic if the cigar string is not valid.
    impl From<&str> for Cigar {

        fn from(str: &str) -> Self {
            let mut operations = Vec::new();
            let mut length = 0 as i64;

            for c in str.chars() {
                if c.is_ascii_digit() {
                    length = length * 10 + c.to_digit(10).unwrap() as i64;
                } else {
                    let op = match c {
                        'M' => CigarOperation::Match(length),
                        'I' => CigarOperation::Insertion(length),
                        'D' => CigarOperation::Deletion(length),
                        'N' => CigarOperation::Nskipped(length),
                        'S' => CigarOperation::Soft(length),
                        'H' => CigarOperation::Hard(length),
                        'P' => CigarOperation::Padded(length),
                        '*' => CigarOperation::Unaligned,
                        _ => panic!("Invalid CIGAR operation"),
                    };
                    operations.push(op);
                    length = 0;
                }
            }
            Cigar{
                cigar: operations
            }
        }
    }


    impl Cigar{
        /// Checks if the CIGAR string contains any skipped regions (N operations).
        ///
        /// Skipped regions typically represent introns in RNA-seq alignments where
        /// the read spans across splice junctions.
        ///
        /// # Returns
        /// `true` if at least one N operation is present, `false` otherwise.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let spliced = Cigar::from_str("35M110N45M").unwrap();
        /// assert_eq!(spliced.has_skipped(), true);
        /// 
        /// let continuous = Cigar::from_str("80M").unwrap();
        /// assert_eq!(continuous.has_skipped(), false);
        /// ```        
        pub fn has_skipped(&self) -> bool{
            self.cigar.iter()
            .any(|e| match e{
                CigarOperation::Nskipped(_) => true,
                _ => false
            })
        }

        pub fn get_read_length_from_cigar(&self) -> i64 {
            let mut res: i64 = 0;
            for cigar_op in self.cigar.iter(){
                    match cigar_op{
                        CigarOperation::Match(n) |  CigarOperation::Insertion(n) | CigarOperation::Soft(n)  => {
                            res += n ;
                        }
                        _  => ()
                    }
                }
            res
        }
        
        /// Returns the positions of all junction boundaries in the reference sequence.
        ///
        /// For each N (skipped) operation, this returns a pair of coordinates:
        /// the end position of the preceding aligned region and the start position of the following aligned region.
        /// This is useful for identifying gaps or junctions in the alignment.
        ///
        /// # Arguments
        /// * `pos` - The starting position of the alignment on the reference sequence
        ///
        /// # Returns
        /// `Some(Vec<i64>)` containing junction coordinates in pairs [region_end, region_start],
        /// or `None` if no N operations are present.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let cigar = Cigar::from_str("35M110N45M").unwrap();
        /// let junctions = cigar.get_junction_position(100);
        /// assert_eq!(junctions, Some(vec![135, 245])); // first region ends at 135, second starts at 245
        /// 
        /// let cigar = Cigar::from_str("35M110N45M10N30M");
        /// let junctions = cigar.get_junction_position(100).unwrap();
        /// assert_eq!(junctions, Some(vec![135, 245, 290, 300])); // two junctions
        /// ```
        pub fn get_skipped_pos_on_ref(&self, pos: i64) -> Option<Vec<i64>>{
            // test skipped so we avoid allocation if we don't need it
            if self.has_skipped(){
                let mut ref_pos = pos; // copy
                let mut results = Vec::new();

                for cigar_op in self.cigar.iter(){
                    // By definition it is impossible to have to consecutive same (N) operation.
                    match cigar_op{
                        CigarOperation::Nskipped(n) => {
                            results.push(ref_pos);
                            ref_pos += n;
                            results.push(ref_pos);
                            },
                        CigarOperation::Match(n) | CigarOperation::Deletion(n) => { ref_pos += n; }, 
                        _  => ()
                    }
                }
                if results.is_empty(){ // should never happend, but I do like fail safe
                    None
                }
                else{
                    Some(results)
                }
            }
            else{
                None
            }
        }

            
        /// Returns the positions of all skipped regions on the reference sequence.
        ///
        /// This is an alias for `get_junction_position` that uses a reference parameter.
        /// Each N operation produces two coordinates: start and end of the skipped region.
        ///
        /// # Arguments
        /// * `pos` - Reference to the starting position of the alignment
        ///
        /// # Returns
        /// `Some(Vec<i64>)` with skipped region boundaries, or `None` if no N operations exist.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let cigar = Cigar::from_str("35M110N45M").unwrap();
        /// let skipped = cigar.get_skipped_pos_on_ref(&100);
        /// assert_eq!(skipped, Some(vec![135, 245]));
        /// ```
        pub fn get_junction_position(&self, pos: i64) -> Option<Vec<i64>>{
            // test skipped so we avoid allocation if we don't need it
            if self.has_skipped(){
                let mut ref_pos = pos; // copy
                let mut results = Vec::new();

                for cigar_op in self.cigar.iter(){
                    /// By definition it is impossible to have to consecutive same (N) operation.
                    match cigar_op{
                        CigarOperation::Nskipped(n) => {
                            results.push(ref_pos);
                            ref_pos += n;
                            results.push(ref_pos);
                            },
                        CigarOperation::Match(n) | CigarOperation::Deletion(n) => { ref_pos += n; }, 
                        _  => ()
                    }
                }
                if results.is_empty(){ // should never happend, but I do like fail safe
                    None
                }
                else{
                    Some(results)
                }
            }
            else{
                None
            }
        }
        /// Returns the number of soft-clipped bases at the appropriate end based on strand.
        ///
        /// For plus strand reads, checks the 3' end (last operation).
        /// For minus strand reads, checks the 5' end (first operation).
        /// This is important for strand-specific analysis of read trimming.
        ///
        /// # Arguments
        /// * `strand` - The strand orientation of the read
        ///
        /// # Returns
        /// `Some(i64)` with the number of soft-clipped bases, or `None` if no soft clipping
        /// at the relevant end.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// use strand_specifier_lib::Strand;
        /// 
        /// let cigar = Cigar::from_str("10S70M20S").unwrap();
        /// assert_eq!(cigar.get_soft_clipped_n(&Strand::Plus), Some(20));   // 3' end
        /// assert_eq!(cigar.get_soft_clipped_n(&Strand::Minus), Some(10));  // 5' end
        /// ```
        pub fn get_soft_clipped_n(&self, strand: &Strand) -> Option<i64>{
            let mut soft_n = None;
            if *strand == Strand::Minus{
                match self.cigar[0]{
                    CigarOperation::Soft (n) => {soft_n = Some(n)},
                    _ => ()
                }

            }
            if *strand == Strand::Plus{
                match self.cigar.last(){
                    Some(CigarOperation::Soft(n)) => {soft_n = Some(*n)},
                    _ => ()
                }

            }
            soft_n
        }

        /// Checks if the read has significant soft clipping at the strand-appropriate end.
        ///
        /// This function is strand-aware: for plus strand reads it checks the 3' end,
        /// for minus strand reads it checks the 5' end. Returns `false` for unstranded reads.
        ///
        /// # Arguments
        /// * `strand` - The strand orientation of the read
        /// * `delta` - Minimum number of soft-clipped bases to consider significant
        ///
        /// # Returns
        /// `true` if soft clipping exceeds the delta threshold, `false` otherwise.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// use strand_specifier_lib::Strand;
        /// 
        /// let cigar = Cigar::from_str("5S70M25S").unwrap();
        /// assert_eq!(cigar.soft_clipped_end(&Strand::Plus, 10), true);   // 25 > 10
        /// assert_eq!(cigar.soft_clipped_end(&Strand::Minus, 10), false); // 5 < 10
        /// assert_eq!(cigar.soft_clipped_end(&Strand::NA, 10), false);    // unstranded
        /// ```
        pub fn soft_clipped_end(&self, strand: &Strand, delta: i64) -> bool{
            match *strand {
                Strand::Minus => {
                    match self.cigar[0]{
                        CigarOperation::Soft (n) => {if n > delta{return true;}},
                        _ => {return false;}
                    };
                },
                Strand::Plus => {
                    match self.cigar.last(){
                        Some(CigarOperation::Soft(n)) => {if *n > delta{return true;}},
                        _ => {return false;}
                    }
                },
                Strand::NA => {
                    return false
                }
            }

            false
        }
        

        /// Checks if a genomic interval is completely covered by a single match operation.
        ///
        /// This function determines whether the specified interval [st, end] (inclusive)
        /// falls entirely within one of the Match (M) operations of the alignment.
        /// Useful for validating that a region of interest has continuous coverage.
        ///
        /// # Arguments
        /// * `pos` - Reference to the starting position of the alignment
        /// * `st` - Start coordinate of the interval to check (inclusive)
        /// * `end` - End coordinate of the interval to check (inclusive)
        ///
        /// # Returns
        /// `true` if the interval is completely covered by a single match operation.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let cigar = Cigar::from_str("100M").unwrap();
        /// assert_eq!(cigar.does_it_match_an_intervall(&500, 520, 580), true);  // within match
        /// assert_eq!(cigar.does_it_match_an_intervall(&500, 520, 620), false); // extends beyond
        /// ```
        pub fn does_it_match_an_intervall(&self, pos: i64, st:i64, end:i64) -> bool{
            let mut ref_pos = pos;
            let mut flag: bool = false;
            for cigar_op in self.cigar.iter(){
                match cigar_op{
                    CigarOperation::Nskipped(n) | CigarOperation::Deletion(n) => {ref_pos += n;},
                    CigarOperation::Match(n) => {
                        if (st >= ref_pos) & (end <= ref_pos + n) {
                           flag = true;
                        }
                        ref_pos += n;
                    },
                    _ => (), // does not consme the reference
                }
            }
            flag
        } 
        /// Calculates the end position of the alignment on the reference sequence.
        ///
        /// This sums all operations that consume reference sequence positions
        /// (M, D, N operations) to determine where the alignment ends.
        ///
        /// # Arguments
        /// * `pos` - Reference to the starting position of the alignment
        ///
        /// # Returns
        /// The final reference coordinate covered by this alignment.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let cigar = Cigar::from_str("50M100N75M").unwrap();
        /// let end = cigar.get_end_of_aln(&1000);
        /// assert_eq!(end, 1225); // 1000 + 50 + 100 + 75
        /// ```
        pub fn get_end_of_aln(&self, pos: i64) -> i64{
            let mut ref_pos = pos;
            for cigar_op in self.cigar.iter(){
                match cigar_op{
                    CigarOperation::Nskipped(n) | CigarOperation::Deletion(n) | CigarOperation::Match(n)=> {ref_pos += n;},
                    _ => (), // does not consme the reference
                }
            }
            ref_pos
        } 


        /// Returns all reference coordinate ranges covered by match operations.
        ///
        /// This function identifies all genomic intervals that are covered by
        /// match (M) operations, skipping over deletions, insertions, and skipped regions.
        /// Returns pairs of coordinates [start, end] for each contiguous matched region.
        ///
        /// # Arguments
        /// * `st` - Starting position of the alignment on the reference
        ///
        /// # Returns
        /// `Vec<i64>` containing alternating start and end coordinates of covered regions.
        ///
        /// # Examples
        /// ```
        /// use CigarParser::cigar::Cigar;
        /// 
        /// let cigar = Cigar::from_str("50M100N75M").unwrap();
        /// let coverage = cigar.get_reference_cover(1000);
        /// assert_eq!(coverage, vec![1000, 1050, 1150, 1225]); // two covered regions
        /// ```
        pub fn get_reference_cover(&self, st:i64)-> Vec<i64>{
            let mut ref_pos = st;
            let mut result : Vec<i64> = Vec::new();
            for cigar_op in self.cigar.iter(){
                match cigar_op{
                CigarOperation::Nskipped(n) | CigarOperation::Deletion(n) => {
                        ref_pos += n;
                    },
                    CigarOperation::Match(n) =>{
                        result.push(ref_pos);
                        result.push(ref_pos + n);
                        ref_pos += n;
                    },
                    _ => ()
                }
            }
            
        result
    }

        
    }

    impl fmt::Display for Cigar {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let mut tobeformatted= Vec::new();
            for op in &self.cigar{
                let x = match op{
                    CigarOperation::Match(length) => format!("{}M", length),
                    CigarOperation::Insertion(length) => format!("{}I", length),
                    CigarOperation::Deletion(length) => format!("{}D", length),
                    CigarOperation::Nskipped(length) => format!("{}N", length),
                    CigarOperation::Soft(length) => format!("{}S", length),
                    CigarOperation::Hard(length) => format!("{}H", length),
                    CigarOperation::Padded(length) => format!("{}P", length),
                    CigarOperation::Unaligned => format!("*"),
                    _ => panic!("Invalid CIGAR operation"),
                };
                tobeformatted.push(x);
            }
            write!(f, "{}", tobeformatted.join(""))
        }
    }


    #[cfg(test)]
    mod tests {
        use crate::cigar::{Cigar, CigarOperation};
        use super::*;
        #[test]
        fn test_from() {
            let cig = Cigar::from_str("35M110N45M3I45M10N").unwrap();
            assert_eq!(cig, Cigar{ cigar: vec![CigarOperation::Match(35), CigarOperation::Nskipped(110), CigarOperation::Match(45), CigarOperation::Insertion(3),
            CigarOperation::Match(45), CigarOperation::Nskipped(10)]});
        }
        #[test]
        #[should_panic]
        fn test_from_panic() {
            let cig = Cigar::from("35M110N45M3I45M10N50K");

        }
        #[test]
        fn test_from_str(){
            let cig = Cigar::from_str("35M110N45M3I45M10N");
            assert_eq!(cig.unwrap(), Cigar{ cigar: vec![CigarOperation::Match(35), CigarOperation::Nskipped(110), CigarOperation::Match(45), CigarOperation::Insertion(3),
                CigarOperation::Match(45), CigarOperation::Nskipped(10)]});
        }
        #[test]
        fn test_pos(){
            let cig = Cigar::from("35M110N45M3I45M10N");
            let results = cig.get_skipped_pos_on_ref(500);
            assert_eq!(results, Some(vec![535, 645, 735, 745]))
        }
        #[test]
        fn test_pos_none(){
            let cig = Cigar::from("35M45M3I45M");
            let results = cig.get_skipped_pos_on_ref(500);
            assert_eq!(results, None)
        }   
        #[test]
        fn test_pos_2(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.get_skipped_pos_on_ref(16946);
            //assert_eq!(results, None)
        }   
        #[test]
        fn test_match_1(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(500, 550, 560);
            //println!("{:?}", results);
            assert_eq!(results, true)
        }   
        #[test]
        fn test_match_outofbound(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(500, 550, 585);
            assert_eq!(results, false)
            //assert_eq!(results, None)
        } 
        #[test]  
        fn test_match_s_eq_e(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(500, 550, 550);
            assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
        #[test]  
        fn soft(){
            let cig = Cigar::from("100M45S");
            let results = cig.soft_clipped_end(&Strand::Plus, 10);
            assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
        #[test]  
        fn test_match_justinbound(){
            let cig = Cigar::from("2S80M53373N169M45S");
            let results = cig.does_it_match_an_intervall(500, 575, 579);
            assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
        #[test]  
        fn softR(){
            let cig = Cigar::from("2S80M53373N169M45S");
            let results = cig.soft_clipped_end(&Strand::Minus, 10);
            assert_eq!(results, false);
            let results = cig.soft_clipped_end(&Strand::Plus, 10);
            assert_eq!(results, true);
            //assert_eq!(results, None)
        }   

        #[test]  
        fn map_region(){
            let cig = Cigar::from("11M214030N240M");
            let res = cig.get_reference_cover(20672897);
            println!("{:?}", res);
            
            //assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
    }
}

