use crate::gates::Gate;
use num_complex::Complex64;
use std::f64::consts::PI;

pub(crate) fn apply_1q_masked(amps: &mut [Complex64], t_mask: usize, c_mask: usize, g: Gate) {
    // Fast paths and general path
    match g {
        Gate::RZ(theta) => {
            let z = Complex64::new(0.0, theta / 2.0).exp();
            let z_c = z.conj();
            let n = amps.len();
            let apply_rz = |i0: usize, i1: usize| {
                amps[i0] *= z_c;
                amps[i1] *= z
            };
            fast_iteration(n, t_mask, c_mask, apply_rz);
        }
        Gate::S => {
            let img = Complex64::new(0.0, 1.0);
            let n = amps.len();
            let apply_s = |_i0: usize, i1: usize| {
                amps[i1] *= img;
            };
            fast_iteration(n, t_mask, c_mask, apply_s);
        }
        Gate::T => {
            let z = Complex64::new(0.0, PI / 4.0).exp();
            let n = amps.len();
            let apply_t = |_i0: usize, i1: usize| {
                amps[i1] *= z;
            };
            fast_iteration(n, t_mask, c_mask, apply_t);
        }
        Gate::X => {
            let n = amps.len();
            let apply_x = |i0: usize, i1: usize| amps.swap(i0, i1);
            fast_iteration(n, t_mask, c_mask, apply_x);
        }
        Gate::Y => {
            let img = Complex64::new(0.0, 1.0);
            let img_conj = img.conj();
            let n = amps.len();
            let apply_y = |i0: usize, i1: usize| {
                amps.swap(i0, i1);
                amps[i0] *= img_conj;
                amps[i1] *= img;
            };
            fast_iteration(n, t_mask, c_mask, apply_y);
        }
        Gate::Z => {
            let n = amps.len();
            let apply_z = |_i0: usize, i1: usize| amps[i1] *= -1.0;
            fast_iteration(n, t_mask, c_mask, apply_z);
        }
        _ => {
            let u = g.create_gate();
            let n = amps.len();
            let apply_generic = |i0: usize, i1: usize| {
                let (a0, a1) = (amps[i0], amps[i1]);
                amps[i0] = u[0][0] * a0 + u[0][1] * a1;
                amps[i1] = u[1][0] * a0 + u[1][1] * a1;
            };
            fast_iteration(n, t_mask, c_mask, apply_generic);
        }
    }
}

/// Iterates through amplitudes and applies closure `apply_gate`.
/// This iteration mechanism avoids checking amplitudes which are excluded by
/// the control qubits (c_mask). The more control qubits |C| the bigger the
/// reduction. Instead of iterating through 2^n, it only iterates through
/// 2^(n - |C|) amplitudes.
fn fast_iteration(n: usize, t_mask: usize, c_mask: usize, mut apply_gate: impl FnMut(usize, usize)) {
    let f = usize::MAX ^ c_mask ^ t_mask;
    let mut s = 0;
    loop {
        apply_gate(s | c_mask, s | c_mask | t_mask);
        s = ((s | c_mask | t_mask) + 1) & f;
        if s == n {
            break;
        }
    }
}
