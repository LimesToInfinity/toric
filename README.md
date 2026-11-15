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
- [x] Gate set: Paulis, S, T, RZ with fast paths
- [ ] 3-qubit repetition code with syndrome extraction
- [ ] Pauli noise channels and logical error rate experiments
- [ ] Pauli strings in symplectic representation
- [ ] Stabilizer (tableau) backend
- [ ] Toric code and decoder
- [ ] Toric UI

## Example

A quantum teleportation circuit is in
[`toric-core/examples/teleport.rs`](toric-core/examples/teleport.rs):

    cargo run --example teleport

The example runs as part of `cargo test`, so it always matches the current API.

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
