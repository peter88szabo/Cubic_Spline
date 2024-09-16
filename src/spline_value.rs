//use crate::interval_index::interval_bracket;
use crate::interval_index::interval_bracket_binarysearch;

pub fn spline_cubic_value(
    t: &[f64],       // knot values
    y: &[f64],       // data values at the knots
    ypp: &[f64],     // second derivatives of the spline at the knots
    tval: f64,       // the point at which the spline is evaluated
    derivative: i32  // decides which value to return: 0 for value, 1 for first derivative, 2 for second derivative
) -> f64 {
    //==============================================================================
    // SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
    //
    // SPLINE_CUBIC_SET must have already been called to define the values of YPP.
    //
    // For any point TVAL in the interval T(LEFT), T(RIGHT), the form of
    // the spline is:
    //
    //   SPL(T) = A + B * (T - T(LEFT))
    //              + C * (T - T(LEFT))^2
    //              + D * (T - T(LEFT))^3
    //
    // Where:
    //   A = Y(LEFT)
    //   B = ( Y(RIGHT) - Y(LEFT) ) / ( T(RIGHT) - T(LEFT) )
    //       - ( YPP(RIGHT) + 2 * YPP(LEFT) ) * ( T(RIGHT) - T(LEFT) ) / 6
    //   C = YPP(LEFT) / 2
    //   D = ( YPP(RIGHT) - YPP(LEFT) ) / ( 6 * ( T(RIGHT) - T(LEFT) ) )
    //
    // Output:
    //   Returns the value of the spline, first derivative, or second derivative based on `derivative`.
    //
    //==============================================================================
    
    // Determine the interval [left, right] that contains tval
    let (left, right) = interval_bracket_binarysearch(t, tval);

    // Compute differences
    let dt = tval - t[left];
    let h = t[right] - t[left];

    match derivative {
        // Case 0: Calculate 0th derivative (the spline value)
        0 => {
            let yval = y[left]
                + dt * ((y[right] - y[left]) / h
                - (ypp[right] / 6.0 + ypp[left] / 3.0) * h
                + dt * (0.5 * ypp[left]
                + dt * ((ypp[right] - ypp[left]) / (6.0 * h))));
            yval
        }
        // Case 1: Return the first derivative
        1 => {
            let ypval = (y[right] - y[left]) / h
                - (ypp[right] / 6.0 + ypp[left] / 3.0) * h
                + dt * (ypp[left] + dt * (0.5 * (ypp[right] - ypp[left]) / h));
            ypval
        }
        // Case 2: Return the second derivative
        2 => {
            let yppval = ypp[left] + dt * (ypp[right] - ypp[left]) / h;
            yppval
        }
        // Default case: return the spline value if derivative is invalid
        _ => {
            panic!("Invalid derivative choice: {}. It must be either 0 (spline value), 1 (first derivative), or 2 (second derivative).", derivative);
        }
    }
}


/*
fn main() {
    let t = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = vec![0.0, 1.0, 4.0, 9.0, 16.0];
    let ypp = vec![0.0, 0.0, 0.0, 0.0, 0.0]; // Example second derivatives (these would normally come from a spline calculation)

    let tval = 2.5;

    // Return the spline value (derivative = 0)
    let spline_value = spline_cubic_value(t.len(), &t, &y, &ypp, tval, 0);
    println!("At t = {}, spline value = {}", tval, spline_value);

    // Return the first derivative (derivative = 1)
    let first_derivative = spline_cubic_value(t.len(), &t, &y, &ypp, tval, 1);
    println!("At t = {}, first derivative = {}", tval, first_derivative);

    // Return the second derivative (derivative = 2)
    let second_derivative = spline_cubic_value(t.len(), &t, &y, &ypp, tval, 2);
    println!("At t = {}, second derivative = {}", tval, second_derivative);
}
*/
