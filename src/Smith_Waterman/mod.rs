// Adapted from Needleman_Wunsch program on GitHub
mod grid;
mod alignment;

use clam::prelude::*;

pub fn align(mut seq1: String, mut seq2: String) -> (Vec<String>, Vec<String>, i32){
    // Get the length
    let len1 = seq1.len() as i32;
    let len2 = seq2.len() as i32;
    // Create the grid
    let (grid, directions) = grid::create_grid(&mut seq1, &mut seq2, len1, len2);
    let high_cell = alignment::highest_cell(&grid);
    // Build and print alignment
    let (aligned_seq1, aligned_seq2) = alignment::build_best_alignment(&grid, &mut directions.clone(), high_cell, seq1.clone(), seq2.clone());
    let score = alignment::score(&aligned_seq1, &aligned_seq2);
    alignment::print_alignments(&aligned_seq1, &aligned_seq2, score);
    return (aligned_seq1, aligned_seq2, score); 
}

// Duplicate of align without print
pub fn align_no_print(mut seq1: String, mut seq2: String) -> i32{
    // Get the length
    let len1 = seq1.len() as i32;
    let len2 = seq2.len() as i32;
    // Create the grid
    let (grid, directions) = grid::create_grid(&mut seq1, &mut seq2, len1, len2);
    let high_cell = alignment::highest_cell(&grid);
    // Build and print alignment
    let (aligned_seq1, aligned_seq2) = alignment::build_best_alignment(&grid, &mut directions.clone(), high_cell, seq1, seq2);
    let score = alignment::score(&aligned_seq1, &aligned_seq2);
    return score; 
}

pub fn clam_align(seq1: String, seq2: String) -> i8{
    return align_no_print(seq1, seq2) as i8;
}

#[derive(Debug)]
/// Implements Smith-Waterman sequence alignments
pub struct SmithWaterman {
    pub seq1: String,
    pub seq2: String
}

impl <T: Number, U: Number> Metric<T, U> for SmithWaterman {
    /// Returns name of distance function
    fn name(&self) -> String {
        "SmithWaterman".to_string()
    }

    /// Returns score of best Smith-Waterman alignment
    /// 
    /// # Arguments
    /// 
    /// * 'x' - An array of a unknown type representing the first sequence
    /// * 'y' - An array of a unknown type representing the second sequence
    fn one_to_one(&self, x: &[T], y: &[T]) -> U {
        let xu = x.iter().map(|v| v.to_be_bytes()[0]).collect();
        let yu = y.iter().map(|v| v.to_be_bytes()[0]).collect();
        let x_str:String = String::from_utf8(xu).unwrap();
        let y_str:String = String::from_utf8(yu).unwrap();
        let score = clam_align(x_str, y_str);
        U::from(score).unwrap()
    }

    /// Returns boolean indictating if this distance function is expensive
    fn is_expensive(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use clam::Metric;
    use crate::Smith_Waterman::grid::Direction;
    use crate::Smith_Waterman::grid::create_grid;
    use crate::Smith_Waterman::alignment::build_best_alignment;
    use crate::Smith_Waterman::alignment::highest_cell;
    use crate::Smith_Waterman::alignment::print_alignments;
    use crate::Smith_Waterman::alignment::score;
    use crate::Smith_Waterman::align;
    use super::SmithWaterman;
    use super::grid;

    struct ImportantExcerpt<'a> {
        part: &'a str,
    }

    #[test]
    fn test1() {
        let grid: Vec<i32> = vec![
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0,
        0, 1, 0, 1, 3, 2, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 2, 2, 3, 2, 1, 0, 1,
        0, 0, 0, 0, 1, 1, 3, 4, 3, 2, 1,
        0, 1, 0, 0, 1, 2, 2, 3, 3, 2, 3,
        0, 0, 2, 1, 0, 1, 1, 2, 2, 2, 2,
        0, 0, 1, 3, 2, 1, 0, 1, 3, 3, 2,
        0, 0, 1, 2, 2, 1, 0, 0, 2, 2, 2
        ];
        let directions:Vec<Direction> = vec![Direction::None, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Left, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalLeft, Direction::Left, Direction::DiagonalUpLeft, Direction::Up,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::DiagonalUp, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Up, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::Up, Direction::DiagonalUp, Direction::Up, Direction::DiagonalUp, Direction::Diagonal, Direction::Up,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Left, Direction::Left, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::DiagonalUp, Direction::Diagonal
        ];
        let mut seq1 : String = "GTCAGGATCT".to_string();
        let mut seq2 : String = "ATCAAGGCCA".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 10, 10);
        // Check values
        for i in 0..120 {
            assert_eq!(grid[i], ftn_grid[i]);
        }
        // Check directions
         for i in 0..119 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![73]);
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn clam_ftn_test1() {
        let seq1 : String = "GTCAGGATCT".to_string();
        let seq2 : String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn test2() {
        let grid = vec![0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 2, 1, 0,
        0, 1, 0, 1, 1, 0,
        0, 0, 0, 0, 2, 2,
        0, 0, 0, 1, 1, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 2, 1,
        ];
        let directions:Vec<Direction> = vec![Direction::None, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, Direction::DiagonalUp,  
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalLeft  
        ];
        let mut seq1 : String = "ATGCAGGA".to_string();
        let mut seq2 : String = "CTGAA".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 8, 5);
        // Check values
        for i in 0..53 {
            assert_eq!(grid[i], ftn_grid[i]);
        }
        // Check directions
        for i in 0..52 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![21, 34, 35, 52]);
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1,&aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["TG", "TGCA", "TGCA", "GA"]);
        assert_eq!(aligned_seq2, vec!["TG", "TG-A", "TGAA", "GA"]);
        assert_eq!(score, 2);
    }
    
    #[test]
    fn clam_ftn_test2() {
        let seq1: String = "ATGCAGGA".to_string();
        let seq2: String = "CTGAA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TG", "TGCA", "TGCA", "GA"]);
        assert_eq!(aligned_seq2, vec!["TG", "TG-A", "TGAA", "GA"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn test3() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 2, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 3, 2, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 3, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 1, 1, 1, 2, 2, 1,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 3, 2, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 4, 3, 2, 2, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 1, 3, 3, 2, 3, 2,
        0, 0, 0, 1, 1, 0, 0, 0, 2, 2, 1, 2, 2, 4, 3, 4,
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 1, 3, 5, 4,
        0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 3, 2, 4, 4,
        0, 0, 2, 1, 0, 0, 2, 1, 0, 1, 1, 0, 2, 2, 3, 3,
        0, 0, 1, 1, 0, 0, 1, 3, 2, 1, 0, 2, 1, 1, 3, 2,
        0, 0, 1, 0, 0, 0, 1, 2, 2, 3, 2, 1, 1, 0, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 1, 1, 3, 4, 3, 2, 1, 1, 1,
        0, 0, 0, 2, 1, 0, 0, 0, 0, 2, 3, 3, 2, 3, 2, 2,
        0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 2, 4, 3, 2, 4, 3,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 2, 3, 3, 2, 3, 3,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 3, 2, 2, 2, 2, 2,
        0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 2, 2, 1, 1, 1, 1,
        ];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Left, Direction::Left, Direction::Diagonal, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalLeft, Direction::Up, Direction::DiagonalUp, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::Left, Direction::Left, Direction::Up, Direction::DiagonalUp, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Left, Direction::Left, Direction::Diagonal, Direction::Left, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::Left, Direction::UpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::DiagonalLeft, Direction::Diagonal, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Up, Direction::DiagonalUp, Direction::Diagonal, Direction::Left, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUp, Direction::DiagonalUp, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Left,
        Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::UpLeft, Direction::Up, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::Up, Direction::DiagonalUp, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Left, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::UpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalLeft, Direction::UpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, Direction::DiagonalUp, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::Left, Direction::Left, Direction::Up, Direction::DiagonalUp, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Diagonal, Direction::Left, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Up, Direction::Up, Direction::Diagonal, Direction::Left, Direction::UpLeft, Direction::Diagonal, Direction::Left, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Up, Direction::Diagonal, Direction::DiagonalLeft, Direction::Up, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::UpLeft, Direction::DiagonalUp, Direction::Diagonal, Direction::Up, Direction::DiagonalUp, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUp, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::DiagonalUp, Direction::DiagonalUp
        ];
        let mut seq1 : String = "AAGTAAGGTGCAGAATGAAA".to_string();
        let mut seq2 : String = "CATTCAGGAAGCTGT".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 20, 15);
        //  Check values
        for i in 0..335 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..334 {
            assert_eq!(directions[i], ftn_directions[i]);
            
        }   
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![174]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AGTAAGGTG"]);
        assert_eq!(aligned_seq2, vec!["AGGAAGCTG"]);
        assert_eq!(score, 5);
    }

    #[test]
    fn clam_ftn_test3() {
        let seq1: String = "GTCAGGATCT".to_string();
        let seq2: String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn test4() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 1, 0, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 0,
        0, 0, 0, 0, 0, 1, 0, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0
        ];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft
        ];
        let mut seq1 : String = "TGACTG".to_string();
        let mut seq2 : String = "AAGGTACAA".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 6, 9);
        //  Check values
        for i in 0..69 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..68 {
            assert_eq!(directions[i], ftn_directions[i]);
        }   
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![47]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AC"]);
        assert_eq!(aligned_seq2, vec!["AC"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn clam_ftn_test4() {
        let seq1 : String = "GTCAGGATCT".to_string();
        let seq2 : String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn test5() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 0, 1,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 2, 1,
        0, 0, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 1, 1,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 2, 1,
        ];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, Direction::Diagonal, 
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, 
        Direction::Up, Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left
        ];
        let mut seq1 : String = "CTAGATGAG".to_string();
        let mut seq2 : String = "TTCAGT".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 9, 6);
        //  Check values
        for i in 0..69 {
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..68 {
            assert_eq!(directions[i], ftn_directions[i]);
        }   
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![33, 48, 68]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AG", "AGAT", "AG"]);
        assert_eq!(aligned_seq2, vec!["AG", "AG-T", "AG"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn clam_ftn_test5() {
        let seq1 : String = "GTCAGGATCT".to_string();
        let seq2 : String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

     #[test]
    fn test6() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0
        ];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft
        ];
        let mut seq1 : String = "TTGATGT".to_string();
        let mut seq2 : String = "AAACTACA".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 7, 8);
        //  Check values
        for i in 0..63{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..62 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![14, 23, 37, 38, 39, 42, 44, 50, 68]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["T", "T", "A", "A", "A", "A", "A", "T", "T"]);
        assert_eq!(aligned_seq2, vec!["T", "T", "A", "A", "A", "A", "A", "T", "T"]);
        assert_eq!(score, 1);
    }

    #[test]
    fn clam_ftn_test6() {
        let seq1 : String = "GTCAGGATCT".to_string();
        let seq2 : String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn test7() {
        let grid = vec![0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0
        ];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft
        ];
        let mut seq1: String = "AAAAA".to_string();
        let mut seq2: String = "TTTTT".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 5, 5);
        //  Check values
        for i in 0..36{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..35 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![-1]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(score, 0);
    }

    #[test]
    fn clam_ftn_test7() {
        let seq1 : String = "GTCAGGATCT".to_string();
        let seq2 : String = "ATCAAGGCCA".to_string();
        let (aligned_seq1, aligned_seq2, score) = align(seq1, seq2);
        assert_eq!(aligned_seq1, vec!["TC-AGG", "TCA-GG"]);
        assert_eq!(aligned_seq2, vec!["TCAAGG", "TCAAGG"]);
        assert_eq!(score, 4);
    }

    #[test]
    fn test8() {
        let grid = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 2,
        0, 1, 2, 1, 0, 1, 2, 2, 1, 0, 1,
        0, 0, 1, 1, 0, 0, 1, 1, 1, 2, 1];
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left,Direction::Left, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::Diagonal,
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Up, 
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal,
        Direction::Up, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Diagonal, Direction::Diagonal, Direction::Left, Direction::DiagonalUpLeft, Direction::DiagonalUp,
        Direction::Up, Direction::DiagonalUpLeft, Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Up, Direction::DiagonalUp, Direction::Diagonal, Direction::Diagonal, Direction::Left
        ];
        let mut seq1: String = "CAGAATATTA".to_string();
        let mut seq2: String = "TTGCTTTGAT".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 10, 10);
        //  Check values
        for i in 0..121{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..120 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![53, 76, 98, 101, 105, 106, 119]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, high_cell, seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["GA", "AT", "GAAT", "AT", "TT", "TT", "TT", "TT-A"]);
        assert_eq!(aligned_seq2, vec!["GA", "AT", "GA-T", "AT", "TT", "TT", "TT", "TTGA"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn test9() {
        let grid = vec![0, 0, 0, 0, 
        0, 1, 0, 0,
        0, 0, 2, 1,
        0, 0, 1, 1,
        0, 0, 0, 2];
        let directions:Vec<Direction> = vec![Direction::None, Direction::Left, Direction::Left, Direction::Left, 
        Direction::Up, Direction::Diagonal, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, 
        Direction::Up, Direction::DiagonalUpLeft, Direction::Diagonal, Direction::Left,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUp, Direction::Diagonal,
        Direction::Up, Direction::DiagonalUpLeft, Direction::DiagonalUpLeft, Direction::Diagonal
        ];
        let mut seq1: String = "ATTG".to_string();
        let mut seq2: String = "ATG".to_string();
        let (ftn_grid, mut ftn_directions) = create_grid(&mut seq1, &mut seq2, 4, 3);
        //  Check values
        for i in 0..20{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..20 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![10, 19]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &mut ftn_directions, vec![10, 19], seq1, seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AT", "TG", "ATTG"]);
        assert_eq!(aligned_seq2, vec!["AT", "TG", "AT-G"]);
        assert_eq!(score, 2);
    }

    #[test]
    fn struct_test1(){
        let seq1_c: String = "GTCAGGATCT".to_string();
        let seq2_c: String = "ATCAAGGCCA".to_string();
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 4 as i8);
    }

    #[test]
    fn struct_test2(){
        let seq1_c: String = "ATGCAGGA".to_string();
        let seq2_c: String = "CTGAA".to_string();
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 2 as i8);
    }

    #[test]
    fn struct_test3(){
        let seq1_c: String = "AAGTAAGGTGCAGAATGAAA".to_string();
        let seq2_c: String = "CATTCAGGAAGCTGT".to_string();        
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 5 as i8);
    }

    #[test]
    fn struct_test4(){
        let seq1_c: String = "TGACTG".to_string();
        let seq2_c: String = "AAGGTACAA".to_string();
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 2 as i8);
    }

    #[test]
    fn struct_test5(){
        let seq1_c: String = "CTAGATGAG".to_string();
        let seq2_c: String = "TTCAGT".to_string();
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 2 as i8);
    }

    #[test]
    fn struct_test6(){
        let seq1_c: String = "TTGATGT".to_string();
        let seq2_c: String = "AAACTACA".to_string();        
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 1 as i8);
    }

    #[test]
    fn struct_test7(){
        let seq1_c: String = "AAAAA".to_string();
        let seq2_c: String = "TTTTT".to_string();
        let seq1: String = seq1_c.clone();
        let seq2: String = seq2_c.clone();
        let metric_sw: SmithWaterman = SmithWaterman{seq1, seq2};
        let score:i8 = metric_sw.one_to_one(seq1_c.as_bytes(), seq2_c.as_bytes());
        assert_eq!(score, 0 as i8);
    }
}