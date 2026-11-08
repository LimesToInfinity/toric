// Quantum teleportation for Hadamard state 1/sqrt(2) * (|0> + |1>)
use num_complex::Complex64;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;
use std::f64::consts::FRAC_1_SQRT_2;
use toric_core::StateVector;
use toric_core::gates::Gate;

fn main() {
    // |000>
    let mut sv = StateVector::new(3);
    // 1/sqrt(2) (|0> + |1>)|00>
    sv.apply_1q(0, Gate::H);
    // 1/sqrt(2) (|0> - |1>)|00>
    sv.apply_1q(0, Gate::Z);
    // 1/2 (|0> - |1>)(|00> + |10>)
    sv.apply_1q(1, Gate::H);
    // 1/2 (|0> - |1>)(|00> + |11>)
    sv.controlled_gate(&[1], 2, Gate::X);
    // 1/2 (|0> (|00> + |11>) - |1> (|10> + |01>))
    sv.controlled_gate(&[0], 1, Gate::X);
    // 1/sqrt(2)^3 (|00> (|0> - |1>) + |10> (|0> + |1>) + |01> (|1> - |0>) + |11> (|0> + |1>))
    sv.apply_1q(0, Gate::H);

    let mut rng = ChaCha8Rng::seed_from_u64(349_021);

    // Apply X gate if result of measuring qubit_1 = 1 (true), else do nothing.
    if sv.measure(1, &mut rng) {
        sv.apply_1q(2, Gate::X);
    }
    // Apply Z gate if result of measuring qubit_0 = 1 (true), else do nothing.
    if sv.measure(0, &mut rng) {
        sv.apply_1q(2, Gate::Z);
    }

    let check_amplitudes = vec![
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        FRAC_1_SQRT_2 * Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        -FRAC_1_SQRT_2 * Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
    ];

    let amp_diff: f64 = check_amplitudes
        .into_iter()
        .enumerate()
        .map(|(i, x)| (x - sv.amplitudes()[i]).norm())
        .sum();

    assert!(amp_diff < 1e-15, "teleportation failed, diff = {amp_diff}");
    println!("Teleportation succeeded: {:?}", sv.amplitudes());
}

#[test]
fn teleport() {
    main();
}
