use std::clone;

// Adapted from Needleman_Wunsch program on GitHub// Transforms the vector of characters into one sequence string
use super::grid::Direction;
fn char_to_string(characters: Vec<char>) -> String {
    let sequence: String = characters.into_iter().collect();
    return sequence;
}

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

fn updateDirections<'a>(directions: &'a mut Vec<Direction>, grid: &'a Vec<i32>, cell: usize, previousDirection: Direction, seq2len: usize) -> &'a mut Vec<Direction>{
    if (cell >= directions.len()){
        return directions;
    }
    if (previousDirection ==  Direction::Diagonal) && (cell >= seq2len + 2) {
        //println!("Diagonal");
        if directions[cell - seq2len - 2] == Direction::None || grid[cell - seq2len - 2] == 0{
            if directions[cell as usize] == Direction::DiagonalUpLeft {
                println!("{}: D->UL ", cell);
                directions[cell as usize] = Direction::UpLeft;
                 
            }
            else if directions[cell as usize] == Direction::DiagonalUp {
                println!("{}: D->U ", cell);
                directions[cell as usize] = Direction::Up;
                 
            }
            else if directions[cell as usize] == Direction::DiagonalLeft {
                println!("{}: D->L ", cell);
                directions[cell as usize] = Direction::Left;
                 
            }
            else {
                println!("{}: D->None ", cell);
                directions[cell as usize] = Direction::None;
                if cell + 1 < grid.len() {
                    let mut directions = updateDirections(directions, grid, cell + 1, Direction::Left, seq2len);
                    if cell + seq2len + 1 < grid.len() {
                        directions = updateDirections(directions, grid, cell + seq2len + 1, Direction::Up, seq2len);
                        if cell + seq2len + 1 < grid.len() {
                            directions = updateDirections(directions, grid, cell + seq2len + 2, Direction::Diagonal, seq2len);
                        }
                    }
                }
                // 
            }
        }
    }
    else if (previousDirection == Direction::Left) && (directions[cell - 1] == Direction::None || grid[cell - 1] == 0) {
        if directions[cell - 1] == Direction::None || grid[cell - 1] == 0 {
            if directions[cell as usize] == Direction::DiagonalUpLeft {
                println!("{}: L->DU ", cell);
                directions[cell as usize] = Direction::DiagonalUp;
                 
            }
            else if directions[cell as usize] == Direction::DiagonalLeft {
                println!("{}: L->D ", cell);
                directions[cell as usize] = Direction::Diagonal;
                 
            }
            else if directions[cell as usize] == Direction::UpLeft {
                println!("{}: L->U ", cell);
                directions[cell as usize] = Direction::Up;
                 
            }
            else {
                println!("{}: L->None ", cell);
                directions[cell as usize] = Direction::None;
                if ((cell + 1) % (seq2len + 1) > 0){
                    let mut directions = updateDirections(directions, grid, cell + 1, previousDirection, seq2len);
                    if cell + seq2len + 1 < grid.len() {
                        directions = updateDirections(directions, grid, cell + seq2len + 1, previousDirection, seq2len);
                        if cell + seq2len + 2 < grid.len() {
                            directions = updateDirections(directions, grid, cell + seq2len + 2, previousDirection, seq2len);
                        }
                    }
                }
                 
            }
        }
    }
    else if (previousDirection == Direction::Up) && (directions[cell as usize] == Direction::None || grid[cell - seq2len] == 0) {
        println!("Up");
        if directions[cell - seq2len - 1] == Direction::None || grid[cell - seq2len - 1] == 0 {
            if directions[cell as usize] == Direction::DiagonalUpLeft {
                println!("{}: U->DL ", cell);
                directions[cell as usize] = Direction::DiagonalLeft;
                 
            }
            else if directions[cell as usize] == Direction::DiagonalUp{
                println!("{}: U->D ", cell);
                directions[cell as usize] = Direction::Diagonal;
                 
            }
            else if directions[cell as usize] == Direction::UpLeft{
                println!("{}: U->L ", cell);
                directions[cell as usize] = Direction::Left;
            }
            else {
                println!("{}: U->None ", cell);
                directions[cell as usize] = Direction::None;
                if cell + 1 < grid.len() {
                    println!("gridlen: {}", grid.len());
                    let mut directions = updateDirections(directions, grid, cell + 1, previousDirection, seq2len);
                    if cell + seq2len + 1 < grid.len() {
                        directions = updateDirections(directions, grid, cell + seq2len + 1, previousDirection, seq2len);
                        if cell + seq2len + 1 < grid.len() {
                            directions = updateDirections(directions, grid, cell + seq2len + 2, previousDirection, seq2len);
                        }
                    }
                }
                 
            }
        }
    }
    //println!("Finished update");
    return directions;
}

// get sequence characters
fn get_seq_char<'a>(cell: i32, seq1: &'a String, seq2: &'a String) -> (char, usize, char, usize) {
    let seq1_char_index = (((cell - (cell % (seq2.len() + 1) as i32)) / (seq2.len() + 1) as i32) - 1) as usize;
    let seq2_char_index = ((cell - 1) % (seq2.len() + 1) as i32) as usize;
    let seq1_char:char = seq1.chars().nth(seq1_char_index).unwrap();
    let seq2_char:char = seq2.chars().nth(seq2_char_index).unwrap();
    println!("{}, {}", seq1_char, seq2_char);
    return (seq1_char, seq1_char_index, seq2_char, seq2_char_index);
}

// Build alignments
fn get_aligned_seq<'a>(new_aligned_seq1: &'a mut Vec<Vec<char>>, new_aligned_seq2: &'a mut Vec<Vec<char>>) -> (Vec<String>, Vec<String>){
    let mut str_aligned_seq1: Vec<String> = Vec::new();
    let mut str_aligned_seq2: Vec<String> = Vec::new();
    for i in 0..(new_aligned_seq1.len() as i32) {
        let char_aligned_seq1:Vec<char> = new_aligned_seq1[i as usize].to_vec().into_iter().collect();
        let char_aligned_seq2:Vec<char> = new_aligned_seq2[i as usize].to_vec().into_iter().collect();
        str_aligned_seq1.append(&mut vec![char_to_string(char_aligned_seq1)]);
        str_aligned_seq2.append(&mut vec![char_to_string(char_aligned_seq2)]);
    }
    return (str_aligned_seq1, str_aligned_seq2);
}

// Find best alignment
fn priv_best_alignment<'a>(score_grid: &'a Vec<i32>, directions: &'a mut Vec<Direction>, previousDirection: Direction, cell: i32, seq1: String, seq2: String, aligned_seq1: &'a mut Vec<Vec<char>>, aligned_seq2: &'a mut Vec<Vec<char>>, aligned_seq_index: &'a mut usize, highcell: i32) -> (Vec<String>, Vec<String>){
    let (seq1_char, seq1_char_index, seq2_char, seq2_char_index) = get_seq_char(cell, &seq1, &seq2) ;
    let mut new_aligned_seq1:&'a mut Vec<Vec<char>> = aligned_seq1;
    let mut new_aligned_seq2: &'a mut Vec<Vec<char>> = aligned_seq2;

    println!("This is cell {}!", cell);
    // If the current direction index includes Diagonal, add the two corresponding characters to the sequence strings    
    if (directions[cell as usize] == Direction::Diagonal) || (directions[cell as usize] == Direction::DiagonalLeft) || (directions[cell as usize] == Direction::DiagonalUp) || (directions[cell as usize] == Direction::DiagonalUpLeft) {
        println!("Diagonal");
        let copy_aligned_seq1:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
        let copy_aligned_seq2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
        copy_aligned_seq1[*aligned_seq_index].insert(0, seq1_char);
        copy_aligned_seq2[*aligned_seq_index].insert(0, seq2_char);
        //println!("192");
        let new_directions = updateDirections(directions, score_grid, cell as usize, Direction::Diagonal, seq2.len());
        // Move to the cell to the diagonally left of the current cell
        if ((cell as i32 - seq2.len() as i32 - 2) as i32 == 0 || (score_grid[(cell as i32 - seq2.len() as i32 - 2) as usize] > 0)) && seq1_char_index > 0 && seq2_char_index > 0{
            println!("Repeat");
            priv_best_alignment(score_grid, new_directions, Direction::Diagonal, cell - seq2.len() as i32 - 2, seq1, seq2, copy_aligned_seq1, copy_aligned_seq2, aligned_seq_index, highcell);
        }
    }
    // If the current direction index includes Up, add a gap to the first sequence and the corresponding character to the second sequence 
    else if (directions[cell as usize] == Direction::Up) || (directions[cell as usize] == Direction::UpLeft) {
        println!("Up");
        let copy_aligned_seq1:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
        let copy_aligned_seq2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
        copy_aligned_seq1[*aligned_seq_index].insert(0, seq1_char);
        copy_aligned_seq2[*aligned_seq_index].insert(0, '-');
        let new_directions = updateDirections(directions, score_grid, cell as usize, Direction::Up, seq2.len());
        // Move to the cell to the diagonally left of the current cell
        if (cell as i32 - seq2.len() as i32 - 1) as i32 >= 0 && (score_grid[cell as usize] > 0) && seq1_char_index > 0 && seq2_char_index > 0{
            priv_best_alignment(score_grid, new_directions, Direction::Up, cell - seq2.len() as i32 - 1, seq1, seq2, copy_aligned_seq1, copy_aligned_seq2, aligned_seq_index, highcell);
        }
    }
    // If the current direction index includes Left, add the corresponing character to the second sequence and a gap to the second sequence
    else if directions[cell as usize] == Direction::Left {
        println!("Left");
        let copy_aligned_seq2:&mut Vec<Vec<char>> = &mut new_aligned_seq2;
        let copy_aligned_seq1:&mut Vec<Vec<char>> = &mut new_aligned_seq1;
        copy_aligned_seq1[*aligned_seq_index].insert(0, '-');
        copy_aligned_seq2[*aligned_seq_index].insert(0, seq2_char);
        // Move to the cell to the diagonally left of the current cell
        let new_directions = updateDirections(directions, score_grid, cell as usize, Direction::Left, seq2.len());
        if (cell - 1) as i32 >= 0 && (score_grid[cell as usize] > 0) && seq1_char_index > 0 && seq2_char_index > 0{
            priv_best_alignment(score_grid, new_directions, Direction::Left, cell - 1, seq1, seq2, copy_aligned_seq1, copy_aligned_seq2, aligned_seq_index, highcell);
        }
    }
    println!("Got alignments");
    // Base case: Turn finished aligned sequences into strings and return them
        return get_aligned_seq(new_aligned_seq1, new_aligned_seq2);
}

    
fn best_alignment<'a>(score_grid: &'a Vec<i32>,directions: &'a mut Vec<Direction>, cell: i32, seq1: String, seq2: String) -> (Vec<String>, Vec<String>){
    let mut index:usize = 0;
    let mut str_aligned_seq1:Vec<String> = vec![];
    let mut str_aligned_seq2:Vec<String> = vec![];
    let mut aligned_sequence1:Vec<Vec<char>> = vec![[].to_vec()];
    let mut aligned_sequence2:Vec<Vec<char>> = vec![[].to_vec()];
    let (mut group_aligned_seq1, mut group_aligned_seq2) = priv_best_alignment(score_grid, directions, Direction::None, cell, seq1, seq2, &mut aligned_sequence1, &mut aligned_sequence2, &mut index, cell);
    str_aligned_seq1.append(&mut group_aligned_seq1);
    str_aligned_seq2.append(&mut group_aligned_seq2);
    for i in 0..str_aligned_seq1.len() {
        let mut j = str_aligned_seq1[i].len();
        while j < str_aligned_seq1[0].len() {
            let char1 = str_aligned_seq1[i - 1].chars().nth(j).unwrap();
            let char2 = str_aligned_seq2[i - 1].chars().nth(j).unwrap();
            str_aligned_seq1[i].insert(j, char1);
            str_aligned_seq2[i].insert(j, char2);
            j = j + 1;
        }
    }
    return (str_aligned_seq1, str_aligned_seq2);
}

pub fn build_best_alignment <'a> (score_grid: &'a Vec<i32>, directions:&'a mut Vec<Direction>, cell: Vec<i32>, seq1: String, seq2: String) -> (Vec<String>, Vec<String>){
    // if highest cell is -1, return null
    if cell == vec![-1] {
        return (vec![], vec![]);
    } 
    println!("# cells: {}", directions.len());
    let (mut sequences1, mut sequences2): (Vec<String>, Vec<String>) = (vec![], vec![]);
    for i in cell {
        println!("{}", i);
        let mut copy_directions = directions.clone();
        while copy_directions[i as usize] != (Direction::None){
            println!("New sequence");
            let (mut new_sequences1, mut new_sequences2) = best_alignment(score_grid, &mut copy_directions, i, seq1.clone(), seq2.clone());
            sequences1.append(&mut new_sequences1);
            sequences2.append(&mut new_sequences2);
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
