use num_complex::Complex64;

pub enum Gate {
    H,
    X,
    Z,
}

impl Gate {
    pub fn create_gate(&self) -> [[Complex64; 2]; 2] {
        match self {
            Self::H => [
                [Complex64::new(1.0/2f64.sqrt(), 0.0), Complex64::new(1.0/2f64.sqrt(), 0.0)],
                [Complex64::new(1.0/2f64.sqrt(), 0.0), Complex64::new(-1.0/2f64.sqrt(), 0.0)]
            ],
            Self::X => [
                [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]
            ],
            Self::Z => [
                [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
                [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)]
            ],
        }
    }
}
