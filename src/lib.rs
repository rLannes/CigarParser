
pub mod cigar{
    
    use std::str::FromStr;

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
        /// utility function return if the operation consume the reference
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
        
        pub fn consume_que(&self) -> bool{
            /// utility function return if the operation consume the Query
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

    #[derive(Debug, PartialEq)]
    /// Representation of Cigar Operation 
    /// This is the main structure users interact with.
    /// Right now it does only allow to interact with from(&str) and get_skipped_pos_on_ref().
    ///
    /// Note: This does not check the logic of operation, for example a cigar string starting or ending by N should not be possible,
    /// at the moment this code does not check for it.
    ///
    /// Warning: from(&str) can panic! use from_str(&str) for a Result<> 
    /// 
    /// let cig = Cigar::from("35M110N45M3I45M10N");
    /// assert_eq!(cig.has_skipped(), true);
    /// let cig = Cigar::from("35M45M3I45M");
    /// assert_ne!(cig.has_skipped(), true);
    /// 
    /// let cig = Cigar::from("35M110N45M3I45M10N");
    /// assert_eq!(cig.get_skipped_pos_on_ref(&500), Some([535, 645, 738, 748]));
    pub struct Cigar{
        cigar: Vec<CigarOperation>,
    }

    #[derive(Debug, PartialEq, Eq)]
    pub struct ParseCigarError;
    /// Create a new Cigar struct from a &str. the &str must be a valid cigar string without "X" or "=" operation
    /// Will return an error if the cigar string is not valid.
    impl FromStr for Cigar {

        type Err = ParseCigarError;
        fn from_str(str: &str) -> Result<Self, ParseCigarError> {
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
                        return Err(ParseCigarError);
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
        
        fn has_skipped(&self) -> bool{
            self.cigar.iter()
            .any(|e| match e{
                CigarOperation::Nskipped(_) => true,
                _ => false
            })
        }
        /// given a Cigar string and the start of the alignment of a read
        /// return all the skipped position junction.
        /// usefull to identify putative splicing junction.
        /// Will not allocate any memory if the Cigar does not contain at least one Skipped operation (N)
        /// This function does not check fo integer overflow, but with i64 there a no genome that come close to this size.
        /// larger genome are in the range of 10**9 and split in chr in the range up to 10**8
        /// i64 is in the range of 10**15.
        pub fn get_skipped_pos_on_ref(&self, pos: &i64) -> Option<Vec<i64>>{
            // test skipped so we avoid allocation if we don't need it
            if self.has_skipped(){
                let mut ref_pos = *pos; // copy
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

        /// this function return true if the reads fully match region defined by st(art) and end.
        /// inclusive of both end
        /// st >= interval, en <= intervall // TODO make end exclusive
        /// st <= end,  st == end should work as expected. 
        pub fn does_it_match_an_intervall(&self, pos: &i64, st:i64, end:i64) -> bool{
            let mut ref_pos = *pos;
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
        pub fn get_end_of_aln(&self, pos: &i64) -> i64{
            let mut ref_pos = *pos;
            for cigar_op in self.cigar.iter(){
                match cigar_op{
                    CigarOperation::Nskipped(n) | CigarOperation::Deletion(n) | CigarOperation::Match(n)=> {ref_pos += n;},
                    _ => (), // does not consme the reference
                }
            }
            ref_pos
        } 

        
    }

    #[cfg(test)]
    mod tests {
        use crate::cigar::{Cigar, CigarOperation};
        use super::*;
        #[test]
        fn test_from() {
            let cig = Cigar::from("35M110N45M3I45M10N");
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
            let results = cig.get_skipped_pos_on_ref(&500);
            assert_eq!(results, Some(vec![535, 645, 735, 745]))
        }
        #[test]
        fn test_pos_none(){
            let cig = Cigar::from("35M45M3I45M");
            let results = cig.get_skipped_pos_on_ref(&500);
            assert_eq!(results, None)
        }   
        #[test]
        fn test_pos_2(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.get_skipped_pos_on_ref(&16946);
            //assert_eq!(results, None)
        }   
        #[test]
        fn test_match_1(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(&500, 550, 560);
            //println!("{:?}", results);
            assert_eq!(results, true)
        }   
        #[test]
        fn test_match_outofbound(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(&500, 550, 580);
            assert_eq!(results, false)
            //assert_eq!(results, None)
        } 
        #[test]  
        fn test_match_s_eq_e(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(&500, 550, 550);
            assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
        #[test]  
        fn test_match_justinbound(){
            let cig = Cigar::from("2S80M53373N169M");
            let results = cig.does_it_match_an_intervall(&500, 575, 579);
            assert_eq!(results, true)
            //assert_eq!(results, None)
        }   
    }
}

