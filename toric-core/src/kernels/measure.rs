use num_complex::Complex64;

/// Calculate probabilities of all qubits in `mask` being 0 (p0) or all of them
/// being 1 (p1).
pub(crate) fn calc_probabilities(amps: &[Complex64], m_mask: usize) -> (f64, f64) {
    let (mut p0, mut p1) = (0.0, 0.0);
    for (i, a) in amps.iter().enumerate() {
        if i & m_mask == 0 {
            p0 += a.norm_sqr()
        } else {
            p1 += a.norm_sqr()
        }
    }
    (p0, p1)
}

/// "Collapses" the `StateVector` according to the measurement `outcome`. A
/// measurement `outcome` of `true` matches with a measurement of 1 and `false`
/// with a measurement of 0.
pub(crate) fn collapse(amps: &mut [Complex64], m_mask: usize, outcome: bool, renorm: f64) {
    if outcome {
        for (i, amp) in amps.iter_mut().enumerate() {
            if i & m_mask == 0 {
                *amp = Complex64::new(0.0, 0.0);
            } else {
                *amp /= renorm;
            }
        }
    } else {
        for (i, amp) in amps.iter_mut().enumerate() {
            if i & m_mask == 0 {
                *amp /= renorm;
            } else {
                *amp = Complex64::new(0.0, 0.0);
            }
        }
    }
}
