pub fn penta(
    a1: &mut [f64],       // A(I,I-2)
    a2: &mut [f64],       // A(I,I-1)
    a3: &mut [f64],       // A(I,I)
    a4: &mut [f64],       // A(I,I+1)
    a5: &mut [f64],       // A(I,I+2)
    b: &mut [f64],        // right-hand side vector
) -> Vec<f64> {           // returns the solution vector
    //==============================================================================
    // PENTA solves a pentadiagonal system of linear equations.
    //
    // Ax = b
    //
    // The diagonal and upper and lower first and second diagonals are:
    //
    //   A(I,I-2) -> A1(I)
    //   A(I,I-1) -> A2(I)
    //   A(I,I)   -> A3(I)
    //   A(I,I+1) -> A4(I)
    //   A(I,I+2) -> A5(I)
    //
    //==============================================================================

    let n = b.len();
    
    let mut x = vec![0.0; n];
    
    let mut xmult: f64;

    // Forward elimination
    for i in 1..n - 1 {
        xmult = a2[i] / a3[i - 1];
        a3[i] -= xmult * a4[i - 1];
        a4[i] -= xmult * a5[i - 1];
        b[i] -= xmult * b[i - 1];

        xmult = a1[i + 1] / a3[i - 1];
        a2[i + 1] -= xmult * a4[i - 1];
        a3[i + 1] -= xmult * a5[i - 1];
        b[i + 1] -= xmult * b[i - 1];
    }

    // Handle the last equation for row n
    xmult = a2[n - 1] / a3[n - 2];
    a3[n - 1] -= xmult * a4[n - 2];
    x[n - 1] = (b[n - 1] - xmult * b[n - 2]) / a3[n - 1];

    // Handle the second-to-last equation for row n-1
    x[n - 2] = (b[n - 2] - a4[n - 2] * x[n - 1]) / a3[n - 2];

    // Back substitution to calculate the rest of the solution
    for i in (0..n - 2).rev() {
        x[i] = (b[i] - a4[i] * x[i + 1] - a5[i] * x[i + 2]) / a3[i];
    }

    x
}

/*

fn main() {

    // Coeffs of pentadiagonal matrix
    let mut a1 = vec![0.0, 0.0, 1.0, 1.0, 1.0];  // A(I,I-2)
    let mut a2 = vec![0.0, 1.0, 1.0, 1.0, 1.0];  // A(I,I-1)
    let mut a3 = vec![2.0, 2.0, 2.0, 2.0, 2.0];  // A(I,I)
    let mut a4 = vec![1.0, 1.0, 1.0, 1.0, 0.0];  // A(I,I+1)
    let mut a5 = vec![1.0, 1.0, 1.0, 0.0, 0.0];  // A(I,I+2)

    // Right-hand side
    let mut b = vec![5.0, 5.0, 5.0, 5.0, 5.0];

    // Solve the pentadiagonal system
    let x = penta(&mut a1, &mut a2, &mut a3, &mut a4, &mut a5, &mut b);

    println!("Solution: {:?}", x);
}
*/
