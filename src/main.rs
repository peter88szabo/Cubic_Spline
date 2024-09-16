mod interval_index;
mod prepare_spline;
mod solve_pentadiagonal;
mod spline_value;

use crate::prepare_spline::set_spline_cubic;
use crate::spline_value::spline_cubic_value;

fn main() {
    // Define 10 points for fitting a sine curve between 0 and 2π
    let t: Vec<f64> = (0..36).map(|i| i as f64 * (2.0 * std::f64::consts::PI / 35.0)).collect();
    let y: Vec<f64> = t.iter().map(|&x| x.sin()).collect();

    // Boundary conditions
    let left_bc_ind = 0;
    let right_bc_ind = 0;
    let left_bc_value = 0.0;
    let right_bc_value = 0.0;


    //set the spline by calculating the second derivatives at the knots
    let ypp = set_spline_cubic(&t, &y, left_bc_ind, left_bc_value, right_bc_ind, right_bc_value);

    // Evaluate the spline at an intermediate point (for example, tval = π)
    let mut tval = std::f64::consts::PI;
    let spline_val = spline_cubic_value(&t, &y, &ypp, tval, 0);  // 0 for spline value

    println!("Spline value at t = π: {}", spline_val);

    // Optionally, you can also compute the first and second derivatives at t = π
    let spline_derivative_1 = spline_cubic_value(&t, &y, &ypp, tval, 1);  // 1 for first derivative
    let spline_derivative_2 = spline_cubic_value(&t, &y, &ypp, tval, 2);  // 2 for second derivative

    println!("Spline first derivative at t = π: {}", spline_derivative_1);
    println!("Spline second derivative at t = π: {}", spline_derivative_2);
    let spline_val = spline_cubic_value(&t, &y, &ypp, tval/2.0, 0);
    println!("Spline value at t = π/2: {}\n", spline_val);

    let n = t.len();

    println!("       t         yknot        yval         y'val        y''val ");
    for i in 0..n{
        tval = t[i];
        let spline_val = spline_cubic_value(&t, &y, &ypp, tval, 0);
        let spline_der_1 = spline_cubic_value(&t, &y, &ypp, tval, 1);
        let spline_der_2 = spline_cubic_value(&t, &y, &ypp, tval, 2);

        println!("{:10.5} {:12.5} {:12.5} {:12.5} {:12.5}", t[i], y[i], spline_val, spline_der_1, spline_der_2);
    }
}
