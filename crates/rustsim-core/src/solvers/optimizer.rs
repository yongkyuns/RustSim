//! Fixed-point iteration optimizer for implicit solvers
//!
//! Provides Anderson acceleration and Newton-Anderson methods for solving
//! implicit update equations in the form x = g(x).

use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

/// Anderson acceleration for fixed-point iteration
///
/// Solves nonlinear equations in fixed-point form x = g(x) by computing the
/// next iterate as a linear combination of previous iterates whose coefficients
/// minimize the least-squares residual.
///
/// # References
/// - Anderson, D. G. (1965). "Iterative Procedures for Nonlinear Integral
///   Equations". Journal of the ACM, 12(4), 547-560.
/// - Walker, H. F., & Ni, P. (2011). "Anderson Acceleration for Fixed-Point
///   Iterations". SIAM Journal on Numerical Analysis, 49(4), 1715-1735.
#[derive(Debug, Clone)]
pub struct Anderson {
    /// Buffer depth (number of stored iterates)
    m: usize,
    /// If true, clear buffer once it reaches depth m
    restart: bool,
    /// Rolling difference buffers
    dx_buffer: VecDeque<DVector<f64>>,
    dr_buffer: VecDeque<DVector<f64>>,
    /// Previous values
    x_prev: Option<DVector<f64>>,
    r_prev: Option<DVector<f64>>,
}

impl Anderson {
    /// Create a new Anderson accelerator
    ///
    /// # Arguments
    /// * `m` - Buffer depth (number of stored iterates)
    /// * `restart` - If true, clear buffer once it reaches depth m
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

    /// Reset the Anderson accelerator
    pub fn reset(&mut self) {
        self.dx_buffer.clear();
        self.dr_buffer.clear();
        self.x_prev = None;
        self.r_prev = None;
    }

    /// Perform one iteration on the fixed-point solution
    ///
    /// # Arguments
    /// * `x` - Current solution
    /// * `g` - Current evaluation of g(x)
    ///
    /// # Returns
    /// * New solution x
    /// * Residual norm (fixed-point error)
    pub fn step(&mut self, x: &DVector<f64>, g: &DVector<f64>) -> (DVector<f64>, f64) {
        // Residual (this gets minimized)
        let res = g - x;
        let res_norm = res.norm();

        // Fallback to regular FPI if m == 0
        if self.m == 0 {
            return (g.clone(), res_norm);
        }

        // If no buffer, regular fixed-point update
        if self.x_prev.is_none() {
            self.x_prev = Some(x.clone());
            self.r_prev = Some(res.clone());
            return (g.clone(), res_norm);
        }

        // Append to difference buffer
        let x_prev = self.x_prev.as_ref().unwrap();
        let r_prev = self.r_prev.as_ref().unwrap();

        self.dx_buffer.push_back(x - x_prev);
        self.dr_buffer.push_back(&res - r_prev);

        // Save values for next iteration
        self.x_prev = Some(x.clone());
        self.r_prev = Some(res.clone());

        // If buffer size m reached, restart
        if self.restart && self.dx_buffer.len() >= self.m {
            self.reset();
            return (g.clone(), res_norm);
        }

        // Get difference matrices
        let m_k = self.dx_buffer.len();
        let n = x.len();

        // Build dR matrix (m_k x n)
        let mut dr_mat = DMatrix::zeros(m_k, n);
        for (i, dr) in self.dr_buffer.iter().enumerate() {
            dr_mat.set_row(i, &dr.transpose());
        }

        // Solve least-squares problem: min ||dR^T * alpha||
        // This gives us the combination coefficients
        let dr_t = dr_mat.transpose();
        let gram = &dr_t * &dr_mat;

        // Solve (dR^T * dR) * alpha = dR^T * r
        let rhs = &dr_t * &res;

        // Use SVD for robust solution
        let alpha = match gram.clone().svd(true, true).solve(&rhs, 1e-10) {
            Ok(alpha_sol) => alpha_sol,
            Err(_) => {
                // Fallback to simple FPI if SVD fails
                return (g.clone(), res_norm);
            }
        };

        // Compute accelerated update
        let mut x_new = g.clone();
        for (i, (dx, dr)) in self.dx_buffer.iter().zip(self.dr_buffer.iter()).enumerate() {
            x_new -= alpha[i] * (dx + dr);
        }

        (x_new, res_norm)
    }
}

impl Default for Anderson {
    fn default() -> Self {
        Self::new(5, false)
    }
}

/// Newton-Anderson accelerator combining Newton's method with Anderson acceleration
///
/// Uses Newton's method when Jacobian is available, falls back to Anderson
/// acceleration otherwise.
#[derive(Debug, Clone)]
pub struct NewtonAnderson {
    anderson: Anderson,
}

impl NewtonAnderson {
    /// Create a new Newton-Anderson accelerator
    pub fn new(m: usize, restart: bool) -> Self {
        Self {
            anderson: Anderson::new(m, restart),
        }
    }

    /// Reset the optimizer
    pub fn reset(&mut self) {
        self.anderson.reset();
    }

    /// Perform one iteration step
    ///
    /// # Arguments
    /// * `x` - Current solution
    /// * `g` - Target value (right-hand side of x = g)
    /// * `jac` - Optional Jacobian matrix (dt * df/dx for implicit solvers)
    ///
    /// # Returns
    /// * New solution x
    /// * Residual norm
    pub fn step(
        &mut self,
        x: &DVector<f64>,
        g: &DVector<f64>,
        jac: Option<&DMatrix<f64>>,
    ) -> (DVector<f64>, f64) {
        if let Some(j) = jac {
            // Newton step: solve (I - J) * dx = g - x
            let n = x.len();
            let identity = DMatrix::identity(n, n);
            let system_matrix = &identity - j;
            let rhs = g - x;

            // Solve using SVD for robustness
            if let Ok(dx) = system_matrix.svd(true, true).solve(&rhs, 1e-10) {
                let x_new = x + dx;
                let residual = (g - &x_new).norm();
                return (x_new, residual);
            }
        }

        // Fallback to Anderson acceleration
        self.anderson.step(x, g)
    }
}

impl Default for NewtonAnderson {
    fn default() -> Self {
        Self::new(5, false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_anderson_simple_fixed_point() {
        // Solve x = cos(x)
        let mut anderson = Anderson::new(3, false);
        let mut x = DVector::from_vec(vec![1.0]);

        for _ in 0..20 {
            let g = DVector::from_vec(vec![x[0].cos()]);
            let (x_new, res) = anderson.step(&x, &g);
            x = x_new;
            if res < 1e-10 {
                break;
            }
        }

        // Verify solution
        assert_relative_eq!(x[0], x[0].cos(), epsilon = 1e-8);
    }

    #[test]
    fn test_newton_anderson_linear() {
        // Solve x = 0.5*x + 1 (solution is x = 2)
        let mut opt = NewtonAnderson::new(3, false);
        let mut x = DVector::from_vec(vec![0.0]);

        // Jacobian for this problem: dg/dx = 0.5
        let jac = DMatrix::from_element(1, 1, 0.5);

        for _ in 0..10 {
            let g = DVector::from_vec(vec![0.5 * x[0] + 1.0]);
            let (x_new, res) = opt.step(&x, &g, Some(&jac));
            x = x_new;
            if res < 1e-10 {
                break;
            }
        }

        assert_relative_eq!(x[0], 2.0, epsilon = 1e-8);
    }
}
