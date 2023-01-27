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
    let mut score_grid: Vec<i32> = Vec::new();
    let total_cells = (len1 + 1) * (len2 + 1);
    let mut directions:Vec<Direction> = Vec::new();
    for _ in 0..total_cells {
        score_grid.push(0);
    }
    // Row label
    let mut i = 0i32;
    while i <= len2 {
        score_grid[i as usize] = 0 - i;
        i = i + 1;
    }
    // Column label
    i = 0;
    while i <= len1 {
        score_grid[(i * (len2 + 1)) as usize] = 0 - i;
        i = i + 1;
    }
    // Middle sections
    i = len2 + 2;
    while i < total_cells {
        directions = max_cell_score(&mut seq1, &mut seq2, &mut score_grid, &i, &mut directions).to_vec();
        if i % (len2 + 1) == len2 as i32 {
            i = i + 1;
        }
        i = i + 1;
    }
    return (score_grid, directions);
}

// Returns direction of the best score for a given cell
fn max_cell_score<'a>(seq1: &'a mut String, seq2: &'a mut String, score_grid: &'a mut Vec<i32>, cell: &'a i32, directions: &'a mut Vec<Direction>) -> &'a mut Vec<Direction> {
    // determine if location is a match
    let seq1_char_index = ((cell - (cell % ((seq2.len() + 1) as i32))) / (seq2.len() + 1) as i32 - 1) as usize;
    let seq2_char_index = (cell % (seq2.len() + 1) as i32 - 1) as usize;
    let seq1_char = seq1.chars().nth(seq1_char_index).unwrap();
    let seq2_char = seq2.chars().nth(seq2_char_index).unwrap();
    let mut match_point = -1;
    // Get surrounding scores
    if seq1_char == seq2_char {
        // Match scoore
        match_point = 1;
    }
    // Gap score
    let from_above:i32 = score_grid[(cell - (seq2.len() + 1) as i32) as usize] - 1;
    let from_left:i32 = score_grid[(cell - 1) as usize] - 1;
    // Base score
    let from_diagonal:i32 = score_grid[(cell - (seq2.len() + 1) as i32 - 1) as usize] + match_point;
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
    return directions;
}
