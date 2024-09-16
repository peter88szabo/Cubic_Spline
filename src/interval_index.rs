#![allow(dead_code)]

//The function loops over the array starting from index 1 up to n-1. 
//It checks where xval lies between two consecutive values
//and returns the indices for the bracket.

pub fn interval_bracket(x: &[f64], xval: f64) -> (usize, usize) {
    let n = x.len();

    if n < 2 {
        panic!("Array length must be at least 2.");
    }

    for i in 1..n - 1 {
        if xval < x[i] {
            return (i - 1, i); // Return the indices for left and right
        }
    }

    // If xval is greater than all elements, return the last bracket
    (n - 2, n - 1)
}

//In terms of speed, the implementation using binary search implmentation
//should faster than a linear search approach owing to the binary search algorithm,
//which has a time complexity of O(log n). In contrast, a linear search
//has a time complexity of O(n).

pub fn interval_bracket_binarysearch(x: &[f64], xval: f64) -> (usize, usize) {
    let n = x.len();

    if xval <= x[0] {
        // If xval is less than or equal to the first element, return the first interval
        return (0, 1);
    } else if xval >= x[n - 1] {
        // If xval is greater than or equal to the last element, return the last interval
        return (n - 2, n - 1);
    }

    // Perform binary search to find the bracketing interval
    match x.binary_search_by(|probe| probe.partial_cmp(&xval).unwrap()) {
        Ok(index) => {
            // Exact match found, adjust to ensure we return a valid interval
            if index == 0 {
                (0, 1) // If the exact match is the first element
            } else if index == n - 1 {
                (n - 2, n - 1) // If the exact match is the last element
            } else {
                (index - 1, index) // Use the previous and next elements as the interval
            }
        }
        Err(index) => {
            // No exact match found, return the bracketing indices
            (index - 1, index)
        }
    }
}



/*
fn main() {
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let xval = 2.5;
    let (left, right) = r8vec_bracket_binarysearch(x.len(), &x, xval);

    println!("Left: {}, Right: {}", left, right);
}
*/
