// Adapted from Needleman_Wunsch program on GitHub
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Direction {
    Diagonal,
    Up,
    Left,
    DiagonalUp,
    DiagonalLeft,
    UpLeft,
    DiagonalUpLeft,
    None
}
struct ImportantExcerpt<'a> {
    part: &'a str,
}
// Creates grid from user entered sequences
pub fn create_grid<'a>(mut seq1: &'a mut String, mut seq2: &'a mut String, len1: i32, len2: i32) -> (Vec<i32>, Vec<Direction>) {
    // Create new vec and set up row and column labels
    let mut score_grid:Vec<i32>= Vec::new();
    let total_cells = (len1 + 1) * (len2 + 1);
    let mut directions:Vec<Direction> = Vec::new();
    for _ in 0..total_cells {
        score_grid.push(0);
    }
    // Middle sections
    let mut i = 0;
    while i < total_cells{
        directions = max_cell_score(&mut seq1, &mut seq2, &mut score_grid, &i, &mut directions).to_vec();
        i = i + 1;
    }
    return (score_grid, directions);
}

// Returns direction of the best score for a given cell
fn max_cell_score<'a>(seq1: &'a mut String, seq2: &'a mut String, score_grid: &'a mut Vec<i32>, cell: &'a i32, directions: &'a mut Vec<Direction>) -> &'a mut Vec<Direction> {
    // Assign every cell but the first in the first row with Left
    if (cell == &0) {
        directions.push(Direction::None);
    }
    else if cell <= &(seq2.len() as i32) {
        directions.push(Direction::Left);
    }
    // Assign every cell in the first column with Up
    else if cell % &(seq2.len() as i32 + 1) == 0 {
        directions.push(Direction::Up);
    }
    else {
        // determine if location is a match
        let seq_1_char_index = (((cell - (cell % (seq2.len() + 1) as i32)) / (seq2.len() + 1) as i32) - 1) as usize;
        let seq_2_char_index = ((cell - 1) % (seq2.len() + 1) as i32) as usize;
        let seq_1_char = seq1.chars().nth(seq_1_char_index).unwrap();
        let seq_2_char = seq2.chars().nth(seq_2_char_index).unwrap();
        let mut match_point = 0;
        // Get surrounding scores
        if seq_1_char == seq_2_char {
            // Match scoore
            match_point = 1;
        }
        else {
            match_point = -1            
        }
        // Prevent negative scores
        // Gap score
        let mut from_above:i32 = score_grid[(cell - (seq2.len() + 1) as i32) as usize] - 1;
        if from_above < 0 {
            from_above = 0;
        }
        let mut from_left:i32 = score_grid[(cell - 1) as usize] - 1;
        if from_left < 0 {
            from_left = 0;
        }
        // Base score
        let mut from_diagonal:i32 = score_grid[(cell - seq2.len() as i32 - 2) as usize] + match_point;
        if from_diagonal < 0 {
            from_diagonal = 0;
        }
        // Find best score
        // Save best directions for each cell
        if (from_diagonal > from_left) && (from_diagonal > from_above) {
            directions.push(Direction::Diagonal);
            score_grid[*cell as usize] = from_diagonal;
        }
        else if (from_above > from_left) && (from_above > from_diagonal) {
            directions.push(Direction::Up);
            score_grid[*cell as usize] = from_above;
        }
        else if (from_left > from_above) && (from_left > from_diagonal) {
            directions.push(Direction::Left);
            score_grid[*cell as usize] = from_left;
        }
        else if (from_above > from_left) && (from_above == from_diagonal) {
            directions.push(Direction::DiagonalUp);
            score_grid[*cell as usize] = from_diagonal;
        }
        else if (from_left > from_above) && (from_left == from_diagonal) {
            directions.push(Direction::DiagonalLeft);
            score_grid[*cell as usize] = from_diagonal;
        }
        else if (from_above == from_left) && (from_above > from_diagonal) {
            directions.push(Direction::UpLeft);
            score_grid[*cell as usize] = from_above;
        }
        else {
            directions.push(Direction::DiagonalUpLeft);
            score_grid[*cell as usize] = from_diagonal;  
        }
    }
    return directions;
}