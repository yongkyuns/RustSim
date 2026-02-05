//! Numerical integration solvers
//!
//! Provides various ODE solvers including:
//! - Explicit Runge-Kutta (RK4, RKDP54, etc.)
//! - Adaptive embedded pairs (RKF21, RKBS32, RKF45, RKCK54, RKV65, RKF78, RKDP87)
//! - Implicit methods (ESDIRK32, ESDIRK4, ESDIRK43, ESDIRK54, ESDIRK85, GEAR family)
//! - Strong Stability Preserving methods (SSPRK22, SSPRK33, SSPRK34)
//! - Backward Differentiation Formulas (BDF2-BDF6, GEAR)
//! - Simple methods (Euler)
//! - Adaptive timestepping with error control

mod base;
// mod bdf;  // TODO: Fix step() implementation
// mod dirk2;  // TODO: Fix borrow issues
// mod dirk3;  // TODO: Fix borrow issues
mod esdirk32;
mod esdirk4;
mod esdirk43;
mod esdirk54;
mod esdirk85;
mod euler;
// mod euler_backward;  // TODO: Fix borrow issues
mod gear;
// mod optimizer;  // TODO: Fix SVD API
mod rk4;
mod rkbs32;
mod rkck54;
mod rkdp54;
mod rkdp87;
mod rkf21;
mod rkf45;
mod rkf78;
mod rkv65;
mod ssprk;

pub use base::*;
// pub use bdf::{BDF2, BDF3, BDF4, BDF5, BDF6};
// pub use dirk2::DIRK2;  // TODO: Fix borrow issues
// pub use dirk3::DIRK3;  // TODO: Fix borrow issues
pub use esdirk32::ESDIRK32;
pub use esdirk4::ESDIRK4;
pub use esdirk43::ESDIRK43;
pub use esdirk54::ESDIRK54;
pub use esdirk85::ESDIRK85;
pub use euler::Euler;
// pub use euler_backward::EulerBackward;  // TODO: Fix borrow issues
pub use gear::{compute_bdf_coefficients, GEAR21, GEAR32, GEAR43, GEAR52A, GEAR54};
// pub use optimizer::{Anderson, NewtonAnderson};  // TODO: Fix SVD API
pub use rk4::RK4;
pub use rkbs32::RKBS32;
pub use rkck54::RKCK54;
pub use rkdp54::RKDP54;
pub use rkdp87::RKDP87;
pub use rkf21::RKF21;
pub use rkf45::RKF45;
pub use rkf78::RKF78;
pub use rkv65::RKV65;
pub use ssprk::{SSPRK22, SSPRK33, SSPRK34};
