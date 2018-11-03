//! Constant Q Transform
//!
//! This is a frequency transform on time-series data.
//! It's like a Fourier transform, except that the bins are
//! logarithmically spaced rather than linearly spaced.
//! Humans perceive pitch in a similar way.
//! Therefore, the CQT is useful for musical analysis.

pub mod window;
use std::f64;

/// Directly calculates the CQT for a single frequency bin, from some samples.
///
/// Currently, samples at the centre of data.
/// Panics if too few samples.
/// Doesn't care about phase.
///
/// * `root_freq`: root frequency (e.g. 220hz)
/// * `bin`: which bin this is (e.g. 2) -- might be 1-indexed in the source literature?
/// * `rate`: sampling frequency (e.g. 44100hz)
/// * `resolution`: bins per octave (e.g. 12)
/// * `data`: handle to array of samples
/// * `window_fn`: a windowing function matching `fn_name(left_bound: usize, right_bound: usize, n: usize) -> f64`
pub fn naive_cqt<W>(root_freq: f64, bin: u8, rate: f64, resolution: u8, data: &[i16], window_fn: &W) -> u16
    where W: Fn(usize, usize, usize) -> f64{
    // first we need to calculate the centre frequency for bin `k`

    let f_k: f64 = calc_f_k(root_freq, bin, resolution);

    let N_k: f64 = calc_N_k(rate, resolution, f_k);

    //println!("data length: {}\t, N_k: {}", data.len(), N_k);

    let left_bound: usize = (data.len() - (N_k.floor() as usize)) / 2;
    let right_bound: usize = (data.len() + (N_k.floor() as usize)) / 2;

    // X_cq[n] = sum{j=left..right}(x[j] * conj(a_k(j - n + N_k/2))
    // and a_k(t) = (1/N_k) * window(t) * exp(-2*pi*i*t*f_k/f_s)
    // and window(t) is some window fn that's zero for samples more than N_k/2 away from n
    // so the complex conjugate of that is the reversal of the
    //  sign of the imaginary part. The conjugate of the sum equals
    //  the sum of the conjugates, so we should be able to safely
    //  ignore it given a real-only everything else.
    // Better question: can we formulate the exponential in such a way
    //  as to only get the amplitude and have that answer work?
    // No. sum of amplitudes doesn't equal amplitude of sum

    let mut sum_real: f64 = 0.0;
    let mut sum_imag: f64 = 0.0;

    for j in left_bound..right_bound {
        let real_common =
            (data[j] as f64) * N_k.recip() * window_fn(left_bound, right_bound, j);
        let complex_common = exp_complex(
            (j as i16 + (N_k.floor() as i16) / 2) - (data.len() / 2) as i16,
            f_k,
            rate,
        );

        sum_real += real_common * complex_common.1;
        sum_imag += real_common * complex_common.0;
    }

    return sum_real.hypot(sum_imag).round() as u16;
}

/// Calculate the value of centre frequency *f_k* in Hz.
///
/// `root_freq * (2.0 ^ (bin / 12.0))` if 12/octave
pub fn calc_f_k(root_freq: f64, bin: u8, resolution: u8) -> f64 {
    // Schoerkhuber and Klapuri have a 1-indexed `bin` here. We should not.
    return root_freq * 2.0_f64.powf( bin as f64 / resolution as f64 );
}

/// Width (N) of the *kth* bin in Hz.
///
/// `(sample rate / ()((2 ^ (1 / resolution)) - 1)) * f_k)`
pub fn calc_N_k(rate: f64, resolution: u8, f_k: f64) -> f64 {
    return rate / ((2.0_f64.powf(1.0 / resolution as f64) - 1.0) * f_k);
}

/// Calculate the complex and real (respectively) values of the exponent in the atom function.
///
/// The atom function: `a_k(t) = (1/N_k) * window(t) * exp(-2*pi*i*t*f_k/f_s)`.
/// This function considers just the exponential part of the atom -
/// the other parts of the atom are purely real.
/// We can
pub fn exp_complex(j: i16, f_k: f64, rate: f64) -> (f64, f64) {
    return (-2.0 * (f64::consts::PI) * (j as f64) * f_k / rate).sin_cos();
}

/// Calculate the CQT over some number of octaves.
///
/// Generally the same parameters as `naive_cqt`, except now we want the number of `octaves`,
/// rather than which `bin` within an octave.
///
/// Loops over `naive_cqt` so you don't have to.

pub fn naive_cqt_octaves<W>(root_freq: f64, octaves: u8, rate: f64, resolution: u8, data: &[i16], window_fn: &W) -> Vec<u16>
    where W: Fn(usize, usize, usize) -> f64{

    // We can't just return an u16 any more.
    // Instead, return a Vec of them.
    let mut output: Vec<u16> = Vec::with_capacity((resolution * octaves) as usize);

    for i in 0..(resolution * octaves){
        output.push(naive_cqt(root_freq, i as u8, rate, resolution, data, window_fn));
    }

    return output;
}

/* *** Tests! *** */

#[cfg(test)]
mod tests {
    use super::*;
    use window;
    use std::f64;

    #[test]
    fn special_values() {
        let root_freq: f64 = 440.0;
        let bin: u8 = 0;
        let rate: f64 = 44100.0;
        let resolution: u8 = 12;

        // 440 * 2 ^ (2/12) = 493.88330125612

        let f_k: f64 = calc_f_k(root_freq, 2, resolution);

        assert_eq!(494.0, f_k.round());

        // Width of bin `k` in Hz
        // precalc'd the 2^(1/12) - 1 = 0.05946309436
        // and since 44100 is preset too
        // 44100 / 0.05946309436 = 741636.480150375
        let N_k_precalc: f64 = 741636.480150375 / f_k;

        assert_eq!(
            N_k_precalc.round(),
            calc_N_k(rate, resolution, f_k).round()
        );
    }

    #[test]
    // show the result of this one with `cargo test -- --nocapture`
    fn waves() {
        use std::f64;

        let root_freq: f64 = 110.0;
        let bin: u8 = 0;
        let rate: f64 = 44100.0;
        let resolution: u8 = 12;

        let mut circular_buffer = [0_i16; 22050]; // half a second's worth of samples

        /*
		// sine combination of ~A5 and ~E6
        println!("test type: sine, A5 and E6");
		for i in 0..22050 {
			circular_buffer[i] = ( (32768.0 / (3.0 + 1.0)) * (
									(i as f64 * 2.0 * f64::consts::PI * 882.0 / rate).sin() * 3.0 +
									(i as f64 * 2.0 * f64::consts::PI * 1320.0 / rate).sin() * 1.0
								) ) as i16;
		}
		*/


        // saw wave combination of all 12 tones in octave
        println!("saw wave combination: (do, mi, le) ");
        for i in 0..22050 {
            let mut addend:f64;
            for j in 0..3 {
                let freq = root_freq * (2.0_f64).powf((j as f64) / 3.0);
                addend = ((32768.0 / 3.0)
                        * (2.0
                        * ((i as f64 * freq / rate) - (0.5 + (i as f64 * freq / rate)).floor())));
                circular_buffer[i] += addend as i16;
            }
        }

        for j in 0..3 {
            let freq = root_freq * (2.0_f64).powf((j as f64) / 3.0);
            println!("freq {}: {}", j, freq);
        }


        /*
		// square waves
        println!("test type: square waves");
		for i in 0..22050 {
			let mut addend = 0.0_f64;
			for j in 0..12 {
				let freq = root_freq * (2.0_f64).powf((j as f64)/12.0);
				addend = ((32768.0 / 12.0) *  (i as f64 * 2.0 * f64::consts::PI * freq / rate).sin().signum() );
				circular_buffer[i] += addend as i16;
			}
		}
		*/

        let mut cb_start = [
            0_i16, 0_i16, 0_i16, 0_i16, 0_i16, 0_i16, 0_i16, 0_i16, 0_i16, 0_i16,
        ];
        for i in 0..10 {
            cb_start[i] = circular_buffer[i];
        }
        println!("Circular buffer, 1st 10 values: {:?}", cb_start);

        println!("analysis root_freq: {}; rate: {}; resolution: {}", root_freq, rate, resolution);
        println!("f_k\tc");

        let vals = naive_cqt_octaves(root_freq, 4, rate, resolution, &circular_buffer, &window::rect);

        for i in 0..48 {
            let c = vals[i];
            println!(
                "{}\t{}",
                calc_f_k(root_freq, i as u8, resolution).round(),
                c
            );
        }
    }
}
