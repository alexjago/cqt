//! Calculate the values of a window function.
//!
//! `left_bound`, `right_bound` and `n` in samples.
//!
//! `n` is sample position within the window,
//! values outside the window are to be zero,
//! `left_bound` is inclusive, `right_bound` not.
//! If you only know the width, set that as the `right_bound`
//! and then set `left_bound = 0`.
use std::f64;

/// Rectangular window
///
/// left_bound`, `right_bound and `n` in samples.
pub fn rect(left_bound: usize, right_bound: usize, n: usize) -> f64 {
    if (n < left_bound) || (n >= right_bound){
        return 0.0;
    } else if n >= right_bound {
        return 0.0;
    } else {
        return 1.0;
    }
}

/// Sine window
///
/// `left_bound`, `right_bound` and `n` in samples
///
/// see [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Sine_window).
pub fn sine(left_bound: usize, right_bound: usize, n: usize) -> f64 {
    let width = right_bound - left_bound;
    if (n < left_bound) || (n >= right_bound){
        return 0.0;
    } else {
        return (f64::consts::PI * n as f64 / ((width - 1) as f64)).sin();
    }
}

/// Hamming 'equiripple' window
///
/// `left_bound`, `right_bound` and `n` in samples
///
/// see [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Hann_and_Hamming_windows).
pub fn hamming(left_bound: usize, right_bound: usize, n: usize) -> f64 {
    let width = right_bound - left_bound;
    if (n < left_bound) || (n >= right_bound){
        return 0.0;
    } else {
        return 0.53836 + 0.46164 * ( (2.0 * f64::consts::PI * n as f64) / ((width - 1) as f64)).cos();
    }
}

/* Unit tests will go here */
