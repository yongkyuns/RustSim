//! PID controller example using the Connection and Simulation API
//!
//! Demonstrates a PID controller regulating a first-order plant
//!
//! System:
//!   Plant: dx/dt = -x + u (first-order lag)
//!   Controller: u = PID(setpoint - x)
//!   Goal: x -> setpoint

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Adder, Constant, Integrator, PID};
use rustsim_core::{Connection, Simulation};

/// Block diagram:
/// ```text
/// setpoint ─┬─> [error] ─> [PID] ─> [plant_gain] ─> [plant] ─> x
///           │                                                    │
///           └────────────────────────────────────────────────────┘
/// ```

fn main() {
    println!("PID Controller - Connection and Simulation API Demo");
    println!("====================================================");
    println!();
    println!("System: First-order plant dx/dt = -x + u");
    println!("Controller: PID with Kp=2.0, Ki=1.0, Kd=0.5");
    println!();

    // Create PID system parameters
    let kp = 2.0;
    let ki = 1.0;
    let kd = 0.5;
    let setpoint_value = 1.0;
    let x0 = 0.0;

    // Create blocks
    let setpoint = Constant::new(setpoint_value);
    let error_calc = Adder::<2>::with_weights([1.0, -1.0]); // setpoint - x
    let controller = PID::new(kp, ki, kd);
    let plant_gain = Adder::<2>::with_weights([-1.0, 1.0]); // -x + u
    let plant = Integrator::new(x0);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(setpoint),    // 0
        Box::new(error_calc),  // 1
        Box::new(controller),  // 2
        Box::new(plant_gain),  // 3
        Box::new(plant),       // 4
    ];

    // Define connections
    // Block diagram:
    // setpoint[0] -> error_calc[0]
    // plant[0] -> error_calc[1], plant_gain[0]
    // error_calc[0] -> controller[0]
    // controller[0] -> plant_gain[1]
    // plant_gain[0] -> plant[0]
    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(),              // setpoint -> error_calc[0]
        Connection::from(4, 0).to(1, 1).to(3, 0).build(),     // plant -> error_calc[1], plant_gain[0]
        Connection::from(1, 0).to(2, 0).build(),              // error_calc -> controller
        Connection::from(2, 0).to(3, 1).build(),              // controller -> plant_gain[1]
        Connection::from(3, 0).to(4, 0).build(),              // plant_gain -> plant
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(0.01);

    println!("Initial setpoint: {}", setpoint_value);
    println!();
    println!(
        "{:>10} {:>12} {:>12} {:>12} {:>12}",
        "Time", "Setpoint", "Output", "Error", "Control"
    );
    println!(
        "{:-<10} {:-<12} {:-<12} {:-<12} {:-<12}",
        "", "", "", "", ""
    );

    // Simulate step response
    let duration = 5.0;
    let dt = sim.dt();
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        // Print every 0.1 seconds
        if i % (steps / 50) == 0 {
            let t = sim.time();
            let sp = sim.get_output(0, 0);  // setpoint
            let y = sim.get_output(4, 0);   // plant output
            let e = sim.get_output(1, 0);   // error
            let u = sim.get_output(2, 0);   // control signal
            println!("{:10.4} {:12.6} {:12.6} {:12.6} {:12.6}", t, sp, y, e, u);
        }

        sim.step();
    }

    println!();
    println!("Step response complete:");
    println!(
        "  Final output: {:.6} (setpoint: {})",
        sim.get_output(4, 0),
        setpoint_value
    );
    println!("  Final error:  {:.6}", sim.get_output(1, 0));
    println!();

    // Test setpoint change
    println!("Changing setpoint to 2.0...");
    println!();

    // Reset and create new simulation with different setpoint
    let setpoint2 = Constant::new(2.0);
    let error_calc2 = Adder::<2>::with_weights([1.0, -1.0]);
    let controller2 = PID::new(kp, ki, kd);
    let plant_gain2 = Adder::<2>::with_weights([-1.0, 1.0]);
    let plant2 = Integrator::new(sim.get_output(4, 0)); // Start from current output

    let blocks2: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(setpoint2),
        Box::new(error_calc2),
        Box::new(controller2),
        Box::new(plant_gain2),
        Box::new(plant2),
    ];

    let connections2 = vec![
        Connection::from(0, 0).to(1, 0).build(),
        Connection::from(4, 0).to(1, 1).to(3, 0).build(),
        Connection::from(1, 0).to(2, 0).build(),
        Connection::from(2, 0).to(3, 1).build(),
        Connection::from(3, 0).to(4, 0).build(),
    ];

    let mut sim2 = Simulation::new(blocks2, connections2).with_dt(0.01);

    println!(
        "{:>10} {:>12} {:>12} {:>12} {:>12}",
        "Time", "Setpoint", "Output", "Error", "Control"
    );
    println!(
        "{:-<10} {:-<12} {:-<12} {:-<12} {:-<12}",
        "", "", "", "", ""
    );

    let duration2 = 5.0;
    let steps2 = (duration2 / dt) as usize;

    for i in 0..=steps2 {
        if i % (steps2 / 50) == 0 {
            let t = sim2.time();
            let sp = sim2.get_output(0, 0);
            let y = sim2.get_output(4, 0);
            let e = sim2.get_output(1, 0);
            let u = sim2.get_output(2, 0);
            println!("{:10.4} {:12.6} {:12.6} {:12.6} {:12.6}", t, sp, y, e, u);
        }

        sim2.step();
    }

    println!();
    println!("Final state:");
    println!("  Output: {:.6} (setpoint: 2.0)", sim2.get_output(4, 0));
    println!("  Error:  {:.6}", sim2.get_output(1, 0));
    println!();
    println!("PID controller successfully regulates the plant!");
}
