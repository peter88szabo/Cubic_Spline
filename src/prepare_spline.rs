//mod solve_pentadiagonal;
use crate::solve_pentadiagonal::penta;

pub fn set_spline_cubic(
    t: &[f64],           // the points where data is specified
    y: &[f64],           // the data values to be interpolated
    left_bc_ind: i32,    // the left boundary condition flag
    left_bc_value: f64,  // the left boundary value, if needed
    right_bc_ind: i32,   // the right boundary condition flag
    right_bc_value: f64, // the right boundary value, if needed
) -> Vec<f64> {
    //==============================================================================
    // SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
    //
    // For data interpolation, the user must call SPLINE_CUBIC_SET to
    // determine the second derivative data, passing in the data to be
    // interpolated, and the desired boundary conditions.
    //
    // The data to be interpolated, plus the SPLINE_CUBIC_SET output,
    // defines the spline. The user may then call SPLINE_CUBIC_VAL to
    // evaluate the spline at any point.
    //
    // The cubic spline is a piecewise cubic polynomial. The intervals
    // are determined by the "knots" or abscissas of the data to be
    // interpolated. The cubic spline has continuous first and second
    // derivatives over the entire interval of interpolation.
    //
    //==============================================================================

    let n = t.len();
    
    if n < 2 {
        panic!("Error in set_spline_cubic()! The number of knots must be at least 2.");
    }

    for i in 0..n - 1 {
        if t[i + 1] <= t[i] {
            panic!("Error in set_spline_cubic() The knots must be strictly increasing.");
        }
    }

    // Initialize arrays for solving the tridiagonal system
    let mut a1 = vec![0.0; n];
    let mut a2 = vec![0.0; n];
    let mut a3 = vec![0.0; n];
    let mut a4 = vec![0.0; n];
    let mut a5 = vec![0.0; n];
    let mut b = vec![0.0; n];

    // Set the first equation
    if left_bc_ind == 0 {
        b[0] = 0.0;
        a3[0] = 1.0;
        a4[0] = -1.0;
    } else if left_bc_ind == 1 {
        b[0] = (y[1] - y[0]) / (t[1] - t[0]) - left_bc_value;
        a3[0] = (t[1] - t[0]) / 3.0;
        a4[0] = (t[1] - t[0]) / 6.0;
    } else if left_bc_ind == 2 {
        b[0] = left_bc_value;
        a3[0] = 1.0;
        a4[0] = 0.0;
    } else if left_bc_ind == 3 {
        b[0] = 0.0;
        a3[0] = -(t[2] - t[1]);
        a4[0] = t[2] - t[0];
        a5[0] = -(t[1] - t[0]);
    } else {
        panic!("ERROR in spline_cubic_set()! The boundary flag left_bc_ind must be 0, 1, 2 or 3.");
    }

    // Set the intermediate equations
    for i in 1..n - 1 {
        b[i] = (y[i + 1] - y[i]) / (t[i + 1] - t[i]) - (y[i] - y[i - 1]) / (t[i] - t[i - 1]);
        a2[i] = (t[i + 1] - t[i]) / 6.0;
        a3[i] = (t[i + 1] - t[i - 1]) / 3.0;
        a4[i] = (t[i] - t[i - 1]) / 6.0;
    }

    // Set the last equation
    if right_bc_ind == 0 {
        b[n - 1] = 0.0;
        a2[n - 1] = -1.0;
        a3[n - 1] = 1.0;
    } else if right_bc_ind == 1 {
        b[n - 1] = right_bc_value - (y[n - 1] - y[n - 2]) / (t[n - 1] - t[n - 2]);
        a2[n - 1] = (t[n - 1] - t[n - 2]) / 6.0;
        a3[n - 1] = (t[n - 1] - t[n - 2]) / 3.0;
    } else if right_bc_ind == 2 {
        b[n - 1] = right_bc_value;
        a2[n - 1] = 0.0;
        a3[n - 1] = 1.0;
    } else if right_bc_ind == 3 {
        b[n - 1] = 0.0;
        a1[n - 1] = -(t[n - 1] - t[n - 2]);
        a2[n - 1] = t[n - 1] - t[n - 3];
        a3[n - 1] = -(t[n - 2] - t[n - 3]);
    } else {
        panic!("ERROR in spline_cubic_set()! The boundary index right_bc_ind must be 0, 1, 2 or 3.");
    }

    // Solve the linear system using the external penta function from solve_pentadiagonal.rs
    let ypp = penta(&mut a1, &mut a2, &mut a3, &mut a4, &mut a5, &mut b);

    // Return the second derivatives
    ypp
}

/*
fn main() {
    let t = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = vec![0.0, 1.0, 4.0, 9.0, 16.0];
    let left_bc_ind = 0;
    let left_bc_value = 0.0;
    let right_bc_ind = 2;
    let right_bc_value = 0.0;

    // Get the second derivatives
    let ypp = spline_cubic_set(&t, &y, left_bc_ind, left_bc_value, right_bc_ind, right_bc_value);

    println!("Second derivatives: {:?}", ypp);
}
*/
