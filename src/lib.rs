//! Constant Q Transform

pub mod cqt {
    use std::f64;

    /// Directly calculates the CQT for a single frequency from some samples.
    ///
    /// Currently, samples at the centre of data.
    /// Panics if too few samples.
    /// Doesn't care about phase.
    ///
    /// * root_freq: root frequency (e.g. 220hz)
    /// * bin: which bin this is (e.g. 2) -- might be 1-indexed in the source literature?
    /// * rate: sampling frequency (e.g. 44100hz)
    /// * resolution: bins per octave (e.g. 12)
    /// * data: handle to array of samples
    pub fn naive_cqt(root_freq: f64, bin: u8, rate: f64, resolution: u8, data: &[i16]) -> u16 {
        // first we need to calculate the centre frequency for bin `k`

        let f_k: f64 = calc_f_k(root_freq, bin, resolution);

        let N_k: f64 = calc_N_k(bin, rate, resolution, f_k);

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
                (data[j] as f64) * (1.0 / N_k) * window((right_bound - left_bound), j);
            let complex_common = exp_complex(
                (j as i16 + (N_k.floor() as i16) / 2) - (data.len() / 2) as i16,
                f_k,
                rate,
            );

            sum_real += real_common * complex_common.1;
            sum_imag += real_common * complex_common.0;
        }

        return (sum_real.powi(2) + sum_imag.powi(2)).sqrt().round() as u16;
    }

    /// Value of centre frequency in Hz.
    /// `root_freq * (2.0 ^ (bin / 12.0))` if 12/octave
    pub fn calc_f_k(root_freq: f64, bin: u8, resolution: u8) -> f64 {
        return root_freq * (2.0_f64).powf(((bin as f64) - 1.0) / (resolution as f64));
    }

    /// Width of bin `k` in Hz
    /// `(sample rate / ((2 ^ (1 / resolution)) - 1)) / f_k`
    pub fn calc_N_k(bin: u8, rate: f64, resolution: u8, f_k: f64) -> f64 {
        return (rate / ((2.0_f64).powf(1.0 / resolution as f64) - 1.0)) / f_k;
    }

    /// Calculate the complex and real (respectively) values of the exponent in the atom function.
    /// Euler's formula, essentially: `exp(i*x) = cos(x) + i*sin(x)`
    /// and `a_k(t) = (1/N_k) * window(t) * exp(-2*pi*i*t*f_k/f_s)`.
    /// Factored out the exponential for the below
    pub fn exp_complex(j: i16, f_k: f64, rate: f64) -> (f64, f64) {
        return (-2.0 * (f64::consts::PI) * (j as f64) * f_k / rate).sin_cos();
    }

    /// Calculate the values of a window function.
    /// `width` and `n` in samples,
    /// `n` is sample position within the window,
    /// values outside the window should be zero,
    /// rectangular for now.
    /// New DIY window: the 'quarter-taper'.
    /// We iterate within bounds, so nothing from here should need to be zero.
    /// A variety of options; uncomment whichever is to be tried.
    pub fn window(width: usize, n: usize) -> f64 {
        return 1.0; // rect
        /*	return (f64::consts::PI * n as f64 / ((width - 1) as f64)).sin(); // sin */
        /*	return 0.54 + 0.46 * ( (2.0 * f64::consts::PI * n as f64) / ((width - 1) as f64)).cos() */
        /*
        // derp
        //	if n < width / 4 {
        //		return 4.0 * n as f64 / width as f64;
        //	} else if n > (3 * width / 4){
        //		return 1.0 - (4.0 * (n - ((3 * width) / 4)) as f64 / width as f64);
        //	} else {
        //		return 1.0;
        //	}
        	*/
    }
}

#[cfg(test)]
mod tests {
    use cqt;
    use std::f64;

    #[test]
    fn special_values() {
        let root_freq: f64 = 440.0;
        let bin: u8 = 0;
        let rate: f64 = 44100.0;
        let resolution: u8 = 12;

        // 440 * 2 ^ ((0 - 1)/12) = 415.3046975799451

        let f_k: f64 = cqt::calc_f_k(root_freq, bin, resolution);

        assert_eq!(415.0, f_k.round());

        // Width of bin `k` in Hz
        // precalc'd the 2^(1/12) - 1 = 0.05946309436
        // and since 44100 is preset too
        // 44100 / 0.05946309436 = 741636.480150375
        let N_k_precalc: f64 = 741636.480150375 / f_k;

        assert_eq!(
            N_k_precalc.round(),
            cqt::calc_N_k(bin, rate, resolution, f_k).round()
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
		for i in 0..22050 {
			circular_buffer[i] = ( (32768.0 / (3.0 + 1.0)) * (
									(i as f64 * 2.0 * f64::consts::PI * 882.0 / rate).sin() * 3.0 +
									(i as f64 * 2.0 * f64::consts::PI * 1320.0 / rate).sin() * 1.0
								) ) as i16;
		}
		*/

        // saw wave combination of all 12 tones in octave
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

        /*
		// square waves
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
        println!("{:?}", cb_start);

        println!("root_freq: {}", root_freq);
        println!("f_k\tc");
        for i in 0..48 {
            let c = cqt::naive_cqt(root_freq, i as u8, rate, resolution, &circular_buffer);
            println!(
                "{}\t{}",
                cqt::calc_f_k(root_freq, i as u8, resolution).round(),
                c
            );
        }
    }
}
