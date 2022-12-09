
use std::env;

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