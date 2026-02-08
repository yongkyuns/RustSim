//! Anderson acceleration for fixed-point iteration

use nalgebra::DVector;
use std::collections::VecDeque;

/// Anderson acceleration for solving x = g(x)
#[derive(Debug)]
pub struct Anderson {
    m: usize,
    restart: bool,
    dx_buffer: VecDeque<DVector<f64>>,
    dr_buffer: VecDeque<DVector<f64>>,
    x_prev: Option<DVector<f64>>,
    r_prev: Option<DVector<f64>>,
}

impl Anderson {
    /// Create new Anderson accelerator
    ///
    /// # Arguments
    /// * `m` - Buffer depth (number of previous iterates to use)
    /// * `restart` - Whether to restart when buffer is full
    pub fn new(m: usize, restart: bool) -> Self {
        Self {
            m,
            restart,
            dx_buffer: VecDeque::with_capacity(m),
            dr_buffer: VecDeque::with_capacity(m),
            x_prev: None,
            r_prev: None,
        }
    }

    /// Perform one acceleration step
    ///
    /// # Arguments
    /// * `x` - Current iterate
    /// * `g` - g(x) evaluation
    ///
    /// # Returns
    /// Tuple of (accelerated iterate, residual norm)
    pub fn step(&mut self, x: &DVector<f64>, g: &DVector<f64>) -> (DVector<f64>, f64) {
        let res = g - x;
        let res_norm = res.norm();

        // Fallback to regular fixed-point iteration if m == 0
        if self.m == 0 {
            return (g.clone(), res_norm);
        }

        // No buffer yet - regular fixed-point update
        let Some(x_prev) = &self.x_prev else {
            self.x_prev = Some(x.clone());
            self.r_prev = Some(res.clone());
            return (g.clone(), res_norm);
        };
        let r_prev = self.r_prev.as_ref().unwrap();

        // Append to difference buffers
        if self.dx_buffer.len() >= self.m {
            self.dx_buffer.pop_front();
        }
        if self.dr_buffer.len() >= self.m {
            self.dr_buffer.pop_front();
        }

        self.dx_buffer.push_back(x - x_prev);
        self.dr_buffer.push_back(&res - r_prev);

        // Update history
        self.x_prev = Some(x.clone());
        self.r_prev = Some(res.clone());

        // Restart if buffer full
        if self.restart && self.dx_buffer.len() >= self.m {
            self.reset();
            return (g.clone(), res_norm);
        }

        let m_k = self.dx_buffer.len();
        let n = x.len();

        // Scalar case optimization
        if n == 1 {
            let dr_sum: f64 = self.dr_buffer.iter().map(|dr| dr[0] * dr[0]).sum();
            if dr_sum <= 1e-16 {
                return (g.clone(), res_norm);
            }
            let dx_dr: f64 = self
                .dx_buffer
                .iter()
                .zip(self.dr_buffer.iter())
                .map(|(dx, dr)| dx[0] * dr[0])
                .sum();
            let x_new = x - DVector::from_element(1, res[0] * dx_dr / dr_sum);
            return (x_new, res_norm);
        }

        // Build difference matrices for least squares
        let mut d_r = nalgebra::DMatrix::<f64>::zeros(n, m_k);
        let mut d_x = nalgebra::DMatrix::<f64>::zeros(n, m_k);

        for (i, (dx, dr)) in self.dx_buffer.iter().zip(self.dr_buffer.iter()).enumerate() {
            d_r.column_mut(i).copy_from(dr);
            d_x.column_mut(i).copy_from(dx);
        }

        // Solve least squares: dR.T @ c = res
        // Using normal equations: (dR.T @ dR) @ c = dR.T @ res
        let dr_t = d_r.transpose();
        let normal_matrix = &dr_t * &d_r;
        let rhs = &dr_t * &res;

        // Simple pseudo-inverse via Cholesky (for well-conditioned case)
        if let Some(chol) = normal_matrix.clone().cholesky() {
            let c = chol.solve(&rhs);
            let x_new = x - &d_x * &c;
            (x_new, res_norm)
        } else {
            // Fallback to regular iteration if matrix is singular
            (g.clone(), res_norm)
        }
    }

    /// Reset the accelerator state
    pub fn reset(&mut self) {
        self.dx_buffer.clear();
        self.dr_buffer.clear();
        self.x_prev = None;
        self.r_prev = None;
    }
}
