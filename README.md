<p align="center">
  <img src="assets/logo.svg" width="96" alt="toric logo">
</p>

<h1 align="center">toric</h1>

<p align="center">
<b>toric</b> is a workspace of Rust crates for efficiently simulating quantum
circuits.
</p>

<p align="center">
  <!-- <a href="https://github.com/<you>/toric/actions/workflows/ci.yml"><img src="https://github.com/<you>/toric/actions/workflows/ci.yml/badge.svg" alt="CI"></a> -->
  <img src="https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue" alt="License: MIT OR Apache-2.0">
  <img src="https://img.shields.io/badge/status-work%20in%20progress-orange" alt="Status: work in progress">
</p>

---

## Overview

**toric** aims to ease the development of quantum circuits and especially error
correction codes through an elegant API and UI.

> **Status:** early development. The API is not stable yet and the crates are
> not published on crates.io.

## Crates

| Crate | Description | Status |
|---|---|---|
| `toric-core` | State vector simulator, gates, measurement, Pauli strings | in progress |

## Progress

- [x] State vector representation (`StateVector`, little-endian qubit order)
- [x] Projective Z-basis measurement with seeds
- [ ] Gate set: Paulis, H, S, T, rotations, CNOT, CZ with exact fast paths
- [ ] 3-qubit repetition code with syndrome extraction
- [ ] Pauli noise channels and logical error rate experiments
- [ ] Pauli strings in symplectic representation
- [ ] Stabilizer (tableau) backend
- [ ] Toric code and decoder
- [ ] Toric UI

## Example
Simulating a quantum teleportation circuit

```rust
 use toric_core::gates::Gate;
 use toric_core::StateVector;
 use num_complex::Complex64;

 // |000>
 let mut sv = StateVector::new(3);
 // 1/sqrt(2) (|0> + |1>)|00>
 sv.apply_1q(0, Gate::H);
 // 1/sqrt(2) (|0> - |1>)|00>
 sv.apply_1q(0, Gate::Z);
 // 1/2 (|0> - |1>)(|00> + |10>)
 sv.apply_1q(1, Gate::H);
 // 1/2 (|0> - |1>)(|00> + |11>)
 sv.controlled_1q(1, 2, Gate::X);
 // 1/2 (|0> (|00> + |11>) - |1> (|10> + |01>))
 sv.controlled_1q(0, 1, Gate::X);
 // 1/sqrt(2)^3 (|00> (|0> - |1>) + |10> (|0> + |1>) + |01> (|1> - |0>) + |11> (|0> + |1>))
 sv.apply_1q(0, Gate::H);

 // Apply X gate if result of measuring qubit_1 = 1 (true), else do nothing.
 if sv.measure(1, 349_021) {
     sv.apply_1q(2, Gate::X);
 }
 // Apply Z gate if result of measuring qubit_0 = 1 (true), else do nothing.
 if sv.measure(0, 234_904) {
     sv.apply_1q(2, Gate::Z);
 }

 let check_amplitudes = vec![
     Complex64::new(0.0, 0.0),
     Complex64::new(0.0, 0.0),
     1.0/2f64.sqrt() * Complex64::new(1.0, 0.0),
     Complex64::new(0.0, 0.0),
     Complex64::new(0.0, 0.0),
     Complex64::new(0.0, 0.0),
     -1.0/2f64.sqrt() * Complex64::new(1.0, 0.0),
     Complex64::new(0.0, 0.0)];

 assert_eq!(sv.amplitudes(), &check_amplitudes);
```

## Design notes

- **Reproducibility:** all randomness flows through a caller-supplied,
  seedable RNG (`ChaCha8Rng`), so every experiment can be replayed bit-for-bit.
- **Symbolic gates:** gates are values of an enum, not anonymous matrices.
  Matrices are derived on demand, and backends can recognize Clifford gates.
- **Exactness where possible:** Pauli, S, CNOT and CZ are implemented as
  permutations and sign flips, so they introduce no rounding error.

## Building

```sh
cargo build
cargo test
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  <https://www.apache.org/licenses/LICENSE-2.0>)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or
  <https://opensource.org/licenses/MIT>)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms
or conditions.
