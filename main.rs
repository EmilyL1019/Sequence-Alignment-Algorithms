mod Needleman_Wunsch;
mod Smith_Waterman;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let algorithm = &args[1];
    let seq1 = &args[2];
    let seq2= &args[3];
    if algorithm == "Needleman-Wunsch" {
        Needleman_Wunsch::align(seq1.to_string(), seq2.to_string());
    }
    else if algorithm == "Smith-Waterman" {
        Smith_Waterman::align(seq1.to_string(), seq2.to_string());
    }
    else {
        println!("Algorithm is unknown. Please run program again!");
    }
}