// Adapted from Needleman_Wunsch program on GitHub
mod grid;
mod alignment;
pub fn align(mut seq1: String, mut seq2: String) {
    // Get the length
    let len1 = seq1.len() as i32;
    let len2 = seq2.len() as i32;
    // Create the grid
    let (grid, directions) = grid::create_grid(&mut seq1, &mut seq2, len1, len2);
    let high_cell = alignment::highest_cell(&grid);
    // Build and print alignment
    let (aligned_seq1, aligned_seq2) = alignment::build_best_alignment(&grid, &directions, high_cell, &mut seq1, &mut seq2);
    let score = alignment::score(&aligned_seq1, &aligned_seq2);
    alignment::print_alignments(&aligned_seq1, &aligned_seq2, score);
        
}

#[cfg(test)]
mod tests {
    use crate::Smith_Waterman::grid::Direction;
    use crate::Smith_Waterman::grid::create_grid;
    use crate::Smith_Waterman::alignment::build_best_alignment;
    use crate::Smith_Waterman::alignment::highest_cell;
    use crate::Smith_Waterman::alignment::print_alignments;
    use crate::Smith_Waterman::alignment::score;
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
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, 
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 10, 10);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
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
        let directions:Vec<Direction> = vec![Direction::Left, Direction::Left, Direction::Left, Direction::Left, Direction::Left, 
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 8, 5);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1,&aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 20, 15);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AGTAAGGTG"]);
        assert_eq!(aligned_seq2, vec!["AGGAAGCTG"]);
        assert_eq!(score, 5);
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 6, 9);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AC"]);
        assert_eq!(aligned_seq2, vec!["AC"]);
        assert_eq!(score, 2);
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 9, 6);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["AG", "AGAT", "AG"]);
        assert_eq!(aligned_seq2, vec!["AG", "AG-T", "AG"]);
        assert_eq!(score, 2);
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
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 7, 8);
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
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(aligned_seq1, vec!["T", "T", "A", "A", "A", "A", "A", "T", "T"]);
        assert_eq!(aligned_seq2, vec!["T", "T", "A", "A", "A", "A", "A", "T", "T"]);
        assert_eq!(score, 1);
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
        let mut seq1 : String = "TTTTT".to_string();
        let mut seq2 : String = "AAAAA".to_string();
        let (ftn_grid, ftn_directions) = create_grid(&mut seq1, &mut seq2, 5, 5);
        //  Check values
        for i in 0..35{
            assert_eq!(grid[i], ftn_grid[i]); 
        }
        // Check directions
        for i in 0..34 {
            assert_eq!(directions[i], ftn_directions[i]);
        }
        let high_cell = highest_cell(&ftn_grid);
        assert_eq!(high_cell, vec![-1]);
        // Create and check alignments
        let (aligned_seq1, aligned_seq2) = build_best_alignment(&ftn_grid, &ftn_directions, high_cell, &mut seq1, &mut seq2);
        let score = score(&aligned_seq1, &aligned_seq2);
        print_alignments(&aligned_seq1, &aligned_seq2, score);
        assert_eq!(score, 0);
    }
}
