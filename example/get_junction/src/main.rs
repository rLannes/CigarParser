


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of Input file
    #[arg(short, long, required = true)]
    input: String,
    /// Prefix name  to be used for Output file
    #[arg(short, long, required = true)]
    output_file_prefix: String
}




fn main() {

    //let mut bam = IndexedReader::from_path(bam_file).unwrap();
    let mut bam = bam::Reader::from_path(&bam_file).unwrap();
    for r in bam.records() {
        //counter += 1;
        //if counter % 1_000_000 == 0 {
        //    println!("Contig: {}; {} reads done", contig, counter);
        //}
        record = r.unwrap();
        pos_s = record.pos();
        cig = Cigar::from_str(&record.cigar().to_string()).unwrap();
        //pos_e = cig.get_end_of_aln(&pos_s);
        flag = record.flags();
        // QC from bam flag value and mapq
        //if (!check_flag(flag, flag to assert in, flag to exclue)) || (record.mapq() < mapq) {
        //    continue;
        //}
        // compute the junction
        cig.get_junction_position(pos_s);


    }
    println!("Hello, world!");

}
