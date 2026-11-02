use num_complex::Complex64;
use std::f64::consts::{PI, FRAC_1_SQRT_2};

pub enum Gate {
    H,
    RX(f64),
    RY(f64),
    RZ(f64),
    S,
    T,
    X,
    Y,
    Z,
}

impl Gate {
    pub fn create_gate(&self) -> [[Complex64; 2]; 2] {
        match self {
            Self::H => [
                [Complex64::new(FRAC_1_SQRT_2, 0.0), Complex64::new(FRAC_1_SQRT_2, 0.0)],
                [Complex64::new(FRAC_1_SQRT_2, 0.0), Complex64::new(-FRAC_1_SQRT_2, 0.0)]
            ],
            Self::RX(theta) => {
                let (s, c) = (theta/2.0).sin_cos();
                [
                    [Complex64::new(c, 0.0), Complex64::new(0.0, -s)],
                    [Complex64::new(0.0, -s), Complex64::new(c, 0.0)]
                ]
            },
            Self::RY(theta) => {
                let (s, c) = (theta/2.0).sin_cos();
                [
                    [Complex64::new(c, 0.0), Complex64::new(-s, 0.0)],
                    [Complex64::new(s, 0.0), Complex64::new(c, 0.0)]
                ]
            },
            Self::RZ(theta) => [
                [Complex64::new(0.0, -theta/2.0).exp(), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(0.0, theta/2.0).exp()]
            ],
            Self::S => [
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(0.0, 1.0)]
            ],
            Self::T => [
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(0.0, PI/4.0).exp()]
            ],
            Self::X => [
                [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]
            ],
            Self::Y => [
                [Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
                [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)]
            ],
            Self::Z => [
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)]
            ],
        }
    }
}
