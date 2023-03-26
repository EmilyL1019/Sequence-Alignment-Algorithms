// Adapted from Needleman_Wunsch program on GitHub
use super::grid::Direction;

// Finds all of the cells with the highest score
pub fn highest_cell<'a>(score_grid: &'a Vec<i32>) -> Vec<i32> {
    let mut high_score: i32= 0;
    let mut highest_cells: Vec<i32> = vec![];
    for i in 0..score_grid.len() {
        // If higher score is found, empty highest cells vector and put current cell inside
        if score_grid[i] > high_score {
            high_score = score_grid[i];
            highest_cells = vec![i as i32];
        }
        // If another cell with the highest score is found, add to the vector
        else if score_grid[i] == high_score {
            highest_cells.append(&mut vec![i as i32]);            
        }
    }
    // If highest cell score = 0, return -1
    if high_score == 0 {
        highest_cells = vec![-1];
    }
    return highest_cells;
}

// Finds the starting cell of alignment
fn lowest_cell<'a>(mut cell: usize, seq1: Vec<char>, seq2: Vec<char>, seq2len: usize) -> usize{
    for i in 0..(seq1.len() - 1) {
        if seq1[seq1.len() - i - 1] == '-' {
            cell = cell - 1;
        }
        else {
            if seq2[seq2.len() - i - 1] == '-' {
                cell = cell - seq2len - 1;
            }
            else {
                cell = cell - seq2len - 2;
            }            
        }
    }
    return cell;
}

fn update_directions<'a>(directions: &'a mut Vec<Direction>, mut cell: usize, seq2len: usize, aligned_seq1: &'a mut Vec<char>, aligned_seq2: &'a mut Vec<char>) -> &'a mut Vec<Direction> {
    //println!("aligned seq1 len = {}", aligned_seq1.len());
    let mut i: usize = 0;
    while i < aligned_seq1.len() {
        // if cell does not have multiple parts move to next cell if possible
        if directions[cell] == Direction::Diagonal || directions[cell] == Direction::Left || directions[cell] == Direction::Up {
            // If at starting cell update direction to none
            if i == aligned_seq1.len() - 1 {
                directions[cell] = Direction::None
            }
            else {
                if aligned_seq1[i + 1] == '-' {
                    cell = cell + 1
                }
                else {
                    if aligned_seq2[i + 1] == '-' {
                        cell = cell + seq2len + 1;
                    }
                    else {
                        cell = cell + seq2len + 2;
                    }
                }
            }
            i = i + 1;
        }
        // If cell with multiple directions is hit, remove one direction and stop editting directions vec
        else {
            if directions[cell] == Direction::DiagonalUpLeft {
                directions[cell] = Direction::UpLeft;
            }
            else if directions[cell] == Direction::DiagonalUp{
                directions[cell] = Direction::Up;
            }
            else if directions[cell] == Direction::DiagonalLeft || directions[cell] == Direction::UpLeft {
                directions[cell] = Direction::Left;
            }
            i = aligned_seq1.len();
        }
    }
    return directions;
}

// get sequence characters
fn get_seq_char<'a>(cell: i32, seq1: &'a String, seq2: &'a String) -> (char, usize, char, usize) {
    let seq1_char_index = (((cell - (cell % (seq2.len() + 1) as i32)) / (seq2.len() + 1) as i32) - 1) as usize;
    let seq2_char_index = ((cell - 1) % (seq2.len() + 1) as i32) as usize;
    let seq1_char:char = seq1.chars().nth(seq1_char_index).unwrap();
    let seq2_char:char = seq2.chars().nth(seq2_char_index).unwrap();
    return (seq1_char, seq1_char_index, seq2_char, seq2_char_index);
}

// Find best alignment
fn priv_best_alignment<'a>(score_grid: &'a Vec<i32>, directions: &'a mut Vec<Direction>, mut cell: i32, seq1: String, seq2: String, aligned_seq1: &'a mut Vec<char>, aligned_seq2: &'a mut Vec<char>) -> (String, String){
    // Look for directions if cell is not 0
    while score_grid[cell as usize] != 0 {
        let (seq1_char, seq1_char_index, seq2_char, seq2_char_index) = get_seq_char(cell, &seq1, &seq2) ;
        // If the current direction index includes Diagonal, add the two corresponding characters to the sequence strings    
        if (directions[cell as usize] == Direction::Diagonal) || (directions[cell as usize] == Direction::DiagonalLeft) || (directions[cell as usize] == Direction::DiagonalUp) || (directions[cell as usize] == Direction::DiagonalUpLeft) {
            aligned_seq1.insert(0, seq1_char);
            aligned_seq2.insert(0, seq2_char);
            // Move to the cell to the diagonally left of the current cell
            if (cell as i32 - seq2.len() as i32 - 2) as i32 >= 0 {
                cell = cell - seq2.len() as i32 - 2;
            }
        }
        // If the current direction index includes Up, add a gap to the first sequence and the corresponding character to the second sequence 
        else if (directions[cell as usize] == Direction::Up) || (directions[cell as usize] == Direction::UpLeft) {
            aligned_seq1.insert(0, seq1_char);
            aligned_seq2.insert(0, '-');
            // Move to the cell to the diagonally left of the current cell
            if (cell as i32 - seq2.len() as i32 - 1) as i32 >= 0 && (score_grid[cell as usize] > 0) && seq1_char_index > 0 && seq2_char_index > 0{
                cell = cell - seq2.len() as i32 - 1;
            }
        }
        // If the current direction index includes Left, add the corresponing character to the second sequence and a gap to the second sequence
        else if directions[cell as usize] == Direction::Left {
            aligned_seq1.insert(0, '-');
            aligned_seq2.insert(0, seq2_char);
            // Move to the cell to the diagonally left of the current cell
            if (cell - 1) as i32 >= 0 && (score_grid[cell as usize] > 0) && seq1_char_index > 0 && seq2_char_index > 0{
                cell = cell - 1;
            }
        }
    }
    // Base case: Turn finished aligned sequences into strings and return them
    let str_aligned_seq1: String = aligned_seq1.to_vec().into_iter().collect();
    let str_aligned_seq2: String = aligned_seq2.to_vec().into_iter().collect();
    return (str_aligned_seq1, str_aligned_seq2);
}

fn best_alignment<'a>(score_grid: &'a Vec<i32>,directions: &'a mut Vec<Direction>, cell: i32, seq1: String, seq2: String) -> (String, String){
    let mut aligned_sequence1:Vec<char> = vec![];
    let mut aligned_sequence2:Vec<char> = vec![];
    return priv_best_alignment(score_grid, directions, cell, seq1, seq2, &mut aligned_sequence1, &mut aligned_sequence2);
}

pub fn build_best_alignment <'a> (score_grid: &'a Vec<i32>, directions:&'a mut Vec<Direction>, cell: Vec<i32>, seq1: String, seq2: String) -> (Vec<String>, Vec<String>){
    // if highest cell is -1, return null
    if cell == vec![-1] {
        return (vec![], vec![]);
    } 
    let (mut sequences1, mut sequences2): (Vec<String>, Vec<String>) = (vec![], vec![]);
    // Start alignment with each highest cell
    for i in cell {
        // Each starting cell gets the original starting grid
        let mut copy_directions = directions.clone();
        // Find every alignment starting from this cell
        while copy_directions[i as usize] != (Direction::None){
            let (new_sequence1, new_sequence2) = best_alignment(score_grid, &mut copy_directions, i, seq1.clone(), seq2.clone());
            sequences1.append(&mut vec![new_sequence1.clone()]);
            sequences2.append(&mut vec![new_sequence2.clone()]);
            let low_cell = lowest_cell(i as usize,new_sequence1.chars().collect(), new_sequence2.chars().collect(), seq2.len());
            // Directions are updated to remove the recently discovered path
            copy_directions = update_directions(&mut copy_directions, low_cell, seq2.len(), &mut new_sequence1.chars().collect(), &mut new_sequence2.chars().collect()).to_vec();
        }
    }
    return (sequences1, sequences2);
}

// Calculates the score
pub fn score<'a>(str_aligned_seq1: &'a Vec<String>, str_aligned_seq2: &'a Vec<String>) -> i32{
    let mut score: i32 = 0;
    if str_aligned_seq1.len() > 0 {
        for i in 0..str_aligned_seq1[0].len(){
            // Check for match and add 1 to score
            if str_aligned_seq1[0].chars().nth(i) == str_aligned_seq2[0].chars().nth(i) {
                score = score + 1;
            }
            else {
                score = score - 1
            }
        }
    }
    return score;
}

// Prints all alignments
pub fn print_alignments<'a>(str_aligned_seq1: &'a Vec<String>, str_aligned_seq2: &'a Vec<String>, score: i32) {
    if str_aligned_seq1.len() == 0 {
        println!("There are no optimal alignments!");
    }
    else {
        println!("Best Alignments:");
        for i in 0..str_aligned_seq1.len() {
            println!("{}      {}", str_aligned_seq1[i], str_aligned_seq2[i]);
        }
    }
    println!("Score: {}", score);
}
