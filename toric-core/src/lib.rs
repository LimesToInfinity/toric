//! # Toric
//! Toric is a quantum circuit simulation software for classical computers.
//! The state vector structure convention used within Toric is little Endian.

use num_complex::Complex64;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;
use rand_core::Rng;

pub mod gates;
use crate::gates::Gate;

#[derive(Debug, PartialEq)]
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
        assert!(q <= Self::MAX_QUBITS, "Only {0} qubits in system", self.n);
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

    /// Applies the gate `g` to qubit `q`.
    /// # Examples
    /// ```
    /// // Applying a Hadamard gate to a single qubit.
    /// # use toric_core::gates::Gate;
    /// # use toric_core::StateVector;
    /// # use num_complex::Complex64;
    ///
    /// let mut sv = StateVector::new(1);
    /// sv.apply_1q(0, Gate::H);
    /// assert_eq!(*sv.amplitudes(), vec![Complex64::new(1.0/2f64.sqrt(), 0.0); 2]);
    /// ```
    pub fn apply_1q(&mut self, q: usize, g: Gate) {
        assert!(q <= self.n, "Only {0} qubits in system", self.n);
        let mask = 1usize << q;
        let u = g.create_gate();
        for i in 0..self.amps.len() {
            if i & mask == 0 {
                let (a0, a1) = (self.amps[i], self.amps[i | mask]);
                self.amps[i]        = u[0][0] * a0 + u[0][1] * a1;
                self.amps[i | mask] = u[1][0] * a0 + u[1][1] * a1;
            }
        }
    }
    /// Applies the unitary operator to target qubit `t` if control qubit `c`
    /// is 1 else does nothing.
    /// # Examples
    /// ```
    /// // Quantum teleportation for Hadamard state 1/sqrt(2) * (|0> + |1>)
    /// # use toric_core::gates::Gate;
    /// # use toric_core::StateVector;
    /// # use num_complex::Complex64;
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
    /// sv.controlled_1q(1, 2, Gate::X);
    /// // 1/2 (|0> (|00> + |11>) - |1> (|10> + |01>))
    /// sv.controlled_1q(0, 1, Gate::X);
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
    /// let check_amplitudes = vec![
    ///     Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
    ///     1.0/2f64.sqrt() * Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0),
    ///     Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0),
    ///     -1.0/2f64.sqrt() * Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    /// assert_eq!(sv.amplitudes(), &check_amplitudes);
    /// ```
    pub fn controlled_1q(&mut self, c: usize, t: usize, g: Gate) {
        assert!(c <= self.n && t <= self.n, "Only {0} qubits in system", self.n);
        assert!(c != t, "Control qubit must differ from target qubit");
        let c_mask = 1usize << c;
        let t_mask = 1usize << t;
        let u = g.create_gate();
        for i in 0..self.amps.len() {
            if (i & c_mask == c_mask) && (i & t_mask == 0) {
                let (a0, a1) = (self.amps[i], self.amps[i | t_mask]);
                self.amps[i]          = u[0][0] * a0 + u[0][1] * a1;
                self.amps[i | t_mask] = u[1][0] * a0 + u[1][1] * a1;
            }
        }
    }

    pub fn swap(&mut self, q1: usize, q2: usize) {
        assert!(q1 <= self.n && q2 <= self.n, "Only {0} qubits in system", self.n);
        let q1_mask = 1usize << q1;
        let q2_mask = 1usize << q2;
        for i in 0..self.amps.len() {
            if (i & q1_mask == 0) && (i & q2_mask != 0) {
                let (a0, a1) = (self.amps[i], self.amps[(i | q1_mask) ^ q2_mask]);
                self.amps[i] = a1;
                self.amps[(i | q1_mask) ^ q2_mask] = a0;
            }
        }
        
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_new_state_vector() {
        let sv = StateVector::new(5);
        assert!(sv.amps[0] == Complex64::new(1.0, 0.0));
    }

    #[test]
    fn test_custom_state_vector() {
        let sv = StateVector::basis(5, 4);
        assert!(sv.amps[4] == Complex64::new(1.0, 0.0));
    }

    #[test]
    fn test_first_qubit_apply_not_gate() {
        let mut sv = StateVector::new(5);
        let not_sv = StateVector::basis(5, 1);
        sv.apply_1q(0, Gate::X);
        assert!(sv.amps == not_sv.amps);
    }
    
    #[test]
    fn test_controlled_cnot() {
        let mut sv = StateVector::basis(5, 1);
        let controlled_not_sv = StateVector::basis(5, 3);
        sv.controlled_1q(0, 1, Gate::X);
        assert!(sv.amps == controlled_not_sv.amps);
    }
    
    #[test]
    fn test_swap() {
        let mut sv = StateVector::basis(5, 1);
        let q1 = 0;
        let q2 = 2;
        let swapped_sv = StateVector::basis(5, 4);
        sv.swap(q1, q2);
        assert!(sv.amps == swapped_sv.amps);
    }

    #[test]
    fn test_measure() {
        let mut sv = StateVector::new(2);
        sv.apply_1q(0, Gate::H);
        sv.apply_1q(1, Gate::H);
        let result = sv.measure(1, 40_317);
        assert!(result);
    }
}
