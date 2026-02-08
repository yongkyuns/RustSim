//! Operator wrappers with linearization support

use nalgebra::{DMatrix, DVector};

/// State of an operator (direct evaluation or linearized)
#[derive(Debug, Clone)]
pub enum OperatorState {
    Direct,
    Linearized {
        x0: DVector<f64>,
        f0: DVector<f64>,
        jac: DMatrix<f64>,
    },
}

/// Algebraic operator: y = f(u)
pub struct Operator<F>
where
    F: Fn(&DVector<f64>) -> DVector<f64>,
{
    func: F,
    state: OperatorState,
}

impl<F> Operator<F>
where
    F: Fn(&DVector<f64>) -> DVector<f64>,
{
    pub fn new(func: F) -> Self {
        Self {
            func,
            state: OperatorState::Direct,
        }
    }

    /// Evaluate the operator
    pub fn call(&self, x: &DVector<f64>) -> DVector<f64> {
        match &self.state {
            OperatorState::Direct => (self.func)(x),
            OperatorState::Linearized { x0, f0, jac } => {
                let dx = x - x0;
                f0 + jac * dx
            }
        }
    }

    /// Compute Jacobian numerically using central differences
    pub fn jacobian(&self, x: &DVector<f64>) -> DMatrix<f64> {
        numerical_jacobian(&self.func, x)
    }

    /// Linearize the operator at point x
    pub fn linearize(&mut self, x: &DVector<f64>) {
        let f0 = (self.func)(x);
        let jac = self.jacobian(x);
        self.state = OperatorState::Linearized {
            x0: x.clone(),
            f0,
            jac,
        };
    }

    /// Reset to direct evaluation
    pub fn reset(&mut self) {
        self.state = OperatorState::Direct;
    }
}

/// Dynamic operator: dx/dt = f(x, u, t)
pub struct DynamicOperator<F>
where
    F: Fn(&DVector<f64>, &DVector<f64>, f64) -> DVector<f64>,
{
    func: F,
    state: DynamicOperatorState,
}

#[derive(Debug, Clone, Default)]
pub struct DynamicOperatorState {
    pub f0: Option<DVector<f64>>,
    pub x0: Option<DVector<f64>>,
    pub u0: Option<DVector<f64>>,
    pub jac_x: Option<DMatrix<f64>>,
    pub jac_u: Option<DMatrix<f64>>,
}

impl<F> DynamicOperator<F>
where
    F: Fn(&DVector<f64>, &DVector<f64>, f64) -> DVector<f64>,
{
    pub fn new(func: F) -> Self {
        Self {
            func,
            state: DynamicOperatorState::default(),
        }
    }

    /// Evaluate the operator
    pub fn call(&self, x: &DVector<f64>, u: &DVector<f64>, t: f64) -> DVector<f64> {
        match &self.state.f0 {
            None => (self.func)(x, u, t),
            Some(f0) => {
                let mut result = f0.clone();

                if let (Some(x0), Some(jx)) = (&self.state.x0, &self.state.jac_x) {
                    let dx = x - x0;
                    result += jx * dx;
                }

                if let (Some(u0), Some(ju)) = (&self.state.u0, &self.state.jac_u) {
                    let du = u - u0;
                    result += ju * du;
                }

                result
            }
        }
    }

    /// Linearize the operator at point (x, u, t)
    pub fn linearize(&mut self, x: &DVector<f64>, u: &DVector<f64>, t: f64) {
        let f0 = (self.func)(x, u, t);
        let jac_x = numerical_jacobian_x(&self.func, x, u, t);
        let jac_u = numerical_jacobian_u(&self.func, x, u, t);

        self.state = DynamicOperatorState {
            f0: Some(f0),
            x0: Some(x.clone()),
            u0: Some(u.clone()),
            jac_x: Some(jac_x),
            jac_u: Some(jac_u),
        };
    }

    /// Reset to direct evaluation
    pub fn reset(&mut self) {
        self.state = DynamicOperatorState::default();
    }
}

/// Compute numerical Jacobian using central differences
pub fn numerical_jacobian<F>(func: &F, x: &DVector<f64>) -> DMatrix<f64>
where
    F: Fn(&DVector<f64>) -> DVector<f64>,
{
    let r = 1e-3;
    let tol = 1e-16;
    let n = x.len();

    let f0 = func(x);
    let m = f0.len();
    let mut jac = DMatrix::<f64>::zeros(m, n);

    for i in 0..n {
        let h = (r * x[i].abs()).max(tol);
        let mut x_plus = x.clone();
        let mut x_minus = x.clone();

        x_plus[i] += h;
        x_minus[i] -= h;

        let f_plus = func(&x_plus);
        let f_minus = func(&x_minus);

        for j in 0..m {
            jac[(j, i)] = (f_plus[j] - f_minus[j]) / (2.0 * h);
        }
    }

    jac
}

/// Compute numerical Jacobian w.r.t. x for dynamic operators
pub fn numerical_jacobian_x<F>(func: &F, x: &DVector<f64>, u: &DVector<f64>, t: f64) -> DMatrix<f64>
where
    F: Fn(&DVector<f64>, &DVector<f64>, f64) -> DVector<f64>,
{
    let r = 1e-3;
    let tol = 1e-16;
    let n = x.len();

    let f0 = func(x, u, t);
    let m = f0.len();
    let mut jac = DMatrix::<f64>::zeros(m, n);

    for i in 0..n {
        let h = (r * x[i].abs()).max(tol);
        let mut x_plus = x.clone();
        let mut x_minus = x.clone();

        x_plus[i] += h;
        x_minus[i] -= h;

        let f_plus = func(&x_plus, u, t);
        let f_minus = func(&x_minus, u, t);

        for j in 0..m {
            jac[(j, i)] = (f_plus[j] - f_minus[j]) / (2.0 * h);
        }
    }

    jac
}

/// Compute numerical Jacobian w.r.t. u for dynamic operators
pub fn numerical_jacobian_u<F>(func: &F, x: &DVector<f64>, u: &DVector<f64>, t: f64) -> DMatrix<f64>
where
    F: Fn(&DVector<f64>, &DVector<f64>, f64) -> DVector<f64>,
{
    let r = 1e-3;
    let tol = 1e-16;
    let n = u.len();

    let f0 = func(x, u, t);
    let m = f0.len();
    let mut jac = DMatrix::<f64>::zeros(m, n);

    for i in 0..n {
        let h = (r * u[i].abs()).max(tol);
        let mut u_plus = u.clone();
        let mut u_minus = u.clone();

        u_plus[i] += h;
        u_minus[i] -= h;

        let f_plus = func(x, &u_plus, t);
        let f_minus = func(x, &u_minus, t);

        for j in 0..m {
            jac[(j, i)] = (f_plus[j] - f_minus[j]) / (2.0 * h);
        }
    }

    jac
}
