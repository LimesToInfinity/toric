//! # Toric
//! Toric is a quantum circuit simulation software for classical computers.
//! The state vector structure convention used within Toric is little Endian.

use num_complex::Complex64;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;
use rand_core::Rng;
use std::f64::consts::PI;

pub mod gates;
use crate::gates::Gate;

#[derive(Clone, Debug, PartialEq)]
/// A quantum state vector with `n` qubits and `2^n` amplitudes.
pub struct StateVector {
    n: usize,
    amps: Vec<Complex64>,
}

/// A quantum state vector in little Endian form
impl StateVector {
    pub const MAX_QUBITS: usize = 30;

    /// Creates a new `StateVector` with `n` qubits and sets amplitude at 
    /// `index` to `1`.
    pub fn basis(n: usize, index: usize) -> Self {
        assert!(n <= Self::MAX_QUBITS, "Too many qubits: {n}");
        assert!(index < 1 << n);
        let mut amps = vec![Complex64::new(0.0, 0.0); 1 << n];
        amps[index] = Complex64::new(1.0, 0.0);
        Self { n, amps }
    }

    /// Getter function: Returns number of qubits in `StateVector`.
    pub fn size(&self) -> usize {
        self.n
    }

    /// Getter function: Returns reference to the amplitudes (`amps`) of
    /// `StateVector`.
    pub fn amplitudes(&self) -> &Vec<Complex64> {
        &self.amps
    }

    /// Creates a new `StateVector` with `n` qubits and probability of `1` for 
    /// the state `|0...0>`.
    pub fn new(n: usize) -> Self {
        Self::basis(n, 0)
    }

    /// Measures qubit `q` in the Z-Basis and uses a random `seed` to collapse
    /// the state. Renormalizes the other states afterwards.
    pub fn measure(&mut self, q: usize, seed: u64) -> bool {
        assert!(q < self.n, "Only {0} qubits in system", self.n);
        let mask = 1usize << q;
        // Probability to measure qubit q in state 1
        let p1: f64 = self.amps.iter().enumerate()
            .filter(|(i, _)| i & mask != 0)
            .map(|(_, a)| a.norm_sqr())
            .sum();
        let mut chacha_prng = ChaCha8Rng::seed_from_u64(seed);
        let mut k = (chacha_prng.next_u64() >> 11) as f64;
        let norm = (1u64 << 53) as f64;
        k /= norm;
        // true: measured 1, false: measured 0
        let outcome = p1 > k;
        if outcome {
            let renorm = p1.sqrt();
            for i in 0..self.amps.len() {
                if i & mask == 0 {
                    self.amps[i] = Complex64::new(0.0, 0.0);
                }
                else {
                    self.amps[i] /= renorm;
                }
            }
        }
        else {
            let renorm = (1.0 - p1).sqrt();
            for i in 0..self.amps.len() {
                if i & mask == 0 {
                    self.amps[i] /= renorm;
                }
                else {
                    self.amps[i] = Complex64::new(0.0, 0.0);
                }
            }
        }
        outcome
    }

    /// Apply gate `g` to qubit `t` respecting the control mask `c_mask` built
    /// from control qubits. The control mask represents which qubits have to
    /// be in state `1` for the gate to be applied.
    pub fn apply_1q_masked(&mut self, t: usize, c_mask: usize, g: Gate) {
        assert!(t < self.n, "Only {0} qubits in system", self.n);
        let t_mask = 1usize << t;
        // Fast paths and general path
        match g {
            Gate::RZ(theta) => {
                let z = Complex64::new(0.0, theta / 2.0).exp();
                let z_c = z.conj();
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps[i] *= z_c;
                        self.amps[i | t_mask] *= z;
                    }
                }
            }
            Gate::S => {
                let img = Complex64::new(0.0, 1.0);
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps[i | t_mask] *= img;
                    }
                }
            },
            Gate::T => {
                let z = Complex64::new(0.0, PI/4.0).exp();
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps[i | t_mask] *= z;
                    }
                }
            },
            Gate::X => {
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps.swap(i, i | t_mask);
                    }
                }
            },
            Gate::Y => {
                let img = Complex64::new(0.0, 1.0);
                let img_conj = img.conj();
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps.swap(i, i | t_mask);
                        self.amps[i] *= img_conj;
                        self.amps[i | t_mask] *= img;
                    }
                }
            },
            Gate::Z => {
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        self.amps[i | t_mask] = self.amps[i | t_mask].scale(-1.0);
                    }
                }
            },
            _ => {
                let u = g.create_gate();
                for i in 0..self.amps.len() {
                    if (i & t_mask == 0) && (i & c_mask == c_mask) {
                        let (a0, a1) = (self.amps[i], self.amps[i | t_mask]);
                        self.amps[i]        = u[0][0] * a0 + u[0][1] * a1;
                        self.amps[i | t_mask] = u[1][0] * a0 + u[1][1] * a1;
                    }
                }
            }
        }
    }
    /// Applies the gate `g` to target qubit `t`.
    /// # Examples
    /// ```
    /// // Applying a Hadamard gate to a single qubit.
    /// # use toric_core::gates::Gate;
    /// # use toric_core::StateVector;
    /// # use num_complex::Complex64;
    /// # use std::f64::consts::FRAC_1_SQRT_2;
    /// let mut sv = StateVector::new(1);
    /// sv.apply_1q(0, Gate::H);
    /// assert_eq!(*sv.amplitudes(), vec![Complex64::new(FRAC_1_SQRT_2, 0.0); 2]);
    /// ```
    pub fn apply_1q(&mut self, t: usize, g: Gate) {
        self.apply_1q_masked(t, 0, g);
    }

    /// Applies the 1-qubit gate `g` to target qubit `t` if control qubits `c`
    /// are all 1 else does nothing.
    /// # Examples
    /// ```
    /// // Quantum teleportation for Hadamard state 1/sqrt(2) * (|0> + |1>)
    /// # use toric_core::gates::Gate;
    /// # use toric_core::StateVector;
    /// # use num_complex::Complex64;
    /// # use std::f64::consts::FRAC_1_SQRT_2;
    ///
    /// // |000>
    /// let mut sv = StateVector::new(3);
    /// // 1/sqrt(2) (|0> + |1>)|00>
    /// sv.apply_1q(0, Gate::H);
    /// // 1/sqrt(2) (|0> - |1>)|00>
    /// sv.apply_1q(0, Gate::Z);
    /// // 1/2 (|0> - |1>)(|00> + |10>)
    /// sv.apply_1q(1, Gate::H);
    /// // 1/2 (|0> - |1>)(|00> + |11>)
    /// sv.controlled_gate(vec![1], 2, Gate::X);
    /// // 1/2 (|0> (|00> + |11>) - |1> (|10> + |01>))
    /// sv.controlled_gate(vec![0], 1, Gate::X);
    /// // 1/sqrt(2)^3 (|00> (|0> - |1>) + |10> (|0> + |1>) + |01> (|1> - |0>) + |11> (|0> + |1>))
    /// sv.apply_1q(0, Gate::H);
    ///
    /// // Apply X gate if result of measuring qubit_1 = 1 (true), else do nothing.
    /// if sv.measure(1, 349_021) {
    ///     sv.apply_1q(2, Gate::X);
    /// }
    /// // Apply Z gate if result of measuring qubit_0 = 1 (true), else do nothing.
    /// if sv.measure(0, 234_904) {
    ///     sv.apply_1q(2, Gate::Z);
    /// }
    ///
    /// let mut check_amplitudes = vec![
    ///     Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
    ///     FRAC_1_SQRT_2 * Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
    ///     Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
    ///     -FRAC_1_SQRT_2 * Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    ///
    /// let amp_diff: f64 = check_amplitudes.into_iter().enumerate()
    ///             .map(|(i, x)| (x - sv.amplitudes()[i]).norm()).sum();
    ///
    /// assert!(amp_diff < 1e-15);
    /// ```
    pub fn controlled_gate(&mut self, c: Vec<usize>, t: usize, g: Gate) {
        assert!(t < self.n, "Only {0} qubits in system", self.n);
        for i in c.iter() {
            assert!(*i != t, "Control qubit must differ from target qubit");
            assert!(*i < self.n, "Only {0} qubits in system", self.n);
        }
        let c_mask: usize = c.iter().fold(0, |m, &i| m | 1 << i);
        self.apply_1q_masked(t, c_mask, g);
    }

    pub fn swap(&mut self, q1: usize, q2: usize) {
        assert!(q1 < self.n && q2 < self.n, "Only {0} qubits in system", self.n);
        assert!(q1 != q2, "Qubit can't be swapped with itself");
        let q1_mask = 1usize << q1;
        let q2_mask = 1usize << q2;
        for i in 0..self.amps.len() {
            if (i & q1_mask == 0) && (i & q2_mask != 0) {
                self.amps.swap(i, (i | q1_mask) ^ q2_mask);
            }
        }
        
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::f64::consts::{PI, FRAC_1_SQRT_2};

    #[test]
    fn fast_paths_rz() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::RZ(PI));
        let sv_check = vec![
            Complex64::new(0.0, -1.0 / 2.0),
            Complex64::new(0.0, 1.0 / 2.0),
            Complex64::new(0.0, -1.0 / 2.0),
            Complex64::new(0.0, 1.0 / 2.0),
        ];
        let amp_diff: f64 = sv_check.into_iter().enumerate()
                    .map(|(i, x)| (x - sv.amplitudes()[i]).norm()).sum();
       
        assert!(amp_diff < 1e-15);
    }

    #[test]
    fn fast_paths_s() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::S);
        let sv_check = vec![
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(0.0, 1.0 / 2.0),
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(0.0, 1.0 / 2.0),
        ];
        assert!(sv.amps == sv_check);
    }

    #[test]
    fn fast_paths_t() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::T);
        let sv_check = vec![
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(0.0, PI / 4.0).exp().scale(1.0 / 2.0),
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(0.0, PI / 4.0).exp().scale(1.0 / 2.0),
        ];
        assert!(sv.amps == sv_check);
    }

    #[test]
    fn fast_paths_x() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(-1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(-1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::X);
        let sv_check = vec![
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(-1.0 / 2.0, 0.0),
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(-1.0 / 2.0, 0.0),
        ];
        assert!(sv.amps == sv_check);
    }

    #[test]
    fn fast_paths_y() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(-1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(-1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::Y);
        let sv_check = vec![
            Complex64::new(0.0, 1.0 / 2.0),
            Complex64::new(0.0, 1.0 / 2.0),
            Complex64::new(0.0, 1.0 / 2.0),
            Complex64::new(0.0, 1.0 / 2.0),
        ];
        assert!(sv.amps == sv_check);
    }

    #[test]
    fn fast_paths_z() {
        let mut sv = StateVector::new(2);
        sv.amps = vec![
            Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(-1.0 / 2.0, 0.0)
            , Complex64::new(1.0 / 2.0, 0.0)
            , Complex64::new(-1.0 / 2.0, 0.0)
        ];
        sv.apply_1q(0, Gate::Z);
        let sv_check = vec![
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(1.0 / 2.0, 0.0),
            Complex64::new(1.0 / 2.0, 0.0),
        ];
        assert!(sv.amps == sv_check);
    }

    #[test]
    fn new_state_vector() {
        let sv = StateVector::new(5);
        assert!(sv.amps[0] == Complex64::new(1.0, 0.0));
    }

    #[test]
    fn custom_state_vector() {
        let sv = StateVector::basis(5, 4);
        assert!(sv.amps[4] == Complex64::new(1.0, 0.0));
    }

    #[test]
    fn cnot() {
        let mut sv = StateVector::basis(5, 1);
        let controlled_not_sv = StateVector::basis(5, 3);
        sv.controlled_gate(vec![0], 1, Gate::X);
        assert!(sv.amps == controlled_not_sv.amps);
    }
    
    #[test]
    fn toffoli() {
        let mut sv = StateVector::new(3);
        sv.apply_1q(0, Gate::H);
        sv.controlled_gate(vec![0], 1, Gate::X);
        // Toffoli gate
        sv.controlled_gate(vec![0, 1], 2, Gate::X);

        let check_amps = vec![
            Complex64::new(FRAC_1_SQRT_2, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(FRAC_1_SQRT_2, 0.0)
        ];
        assert!(sv.amps == check_amps);
    }
    
    #[test]
    fn swap() {
        let mut sv = StateVector::basis(5, 1);
        let q1 = 0;
        let q2 = 2;
        let swapped_sv = StateVector::basis(5, 4);
        sv.swap(q1, q2);
        assert!(sv.amps == swapped_sv.amps);
    }

    #[test]
    fn measure() {
        let mut sv = StateVector::new(2);
        sv.apply_1q(0, Gate::H);
        sv.apply_1q(1, Gate::H);
        let result = sv.measure(1, 40_317);
        assert!(result);
    }
}
