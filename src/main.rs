use std::env;
use clam::prelude::*;



fn main() {
    let args: Vec<String> = env::args().collect();
    let algorithm = &args[1];
    let seq1 = &args[2];
    let seq2= &args[3];
    if algorithm == "Needleman-Wunsch" {
        algorithms_lib::needleman_wunsch::align(seq1.to_string(), seq2.to_string());
    }
    else if algorithm == "Smith-Waterman" {
        algorithms_lib::smith_waterman::align(seq1.to_string(), seq2.to_string());
    }
    else {
        println!("Algorithm is unknown. Please run program again!");
    }
}

#[cfg(test)]
mod test {
use clam::metric_from_name;
    fn test1(){
        let seq1 = "TGATGTG";
        let seq2 = "CAGGTGG";
        let metric = metric_from_name("NeedlemanWunsch", false).unwrap();
        let score:u8 = metric.one_to_one(seq1.as_bytes(), seq2.as_bytes());
        assert_eq!(score, 0);
        //approx_eq!(i64, metric.one_to_one(&seq1, &seq2))
    }
}