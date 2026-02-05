//! Generic function blocks with const generic I/O sizes

use crate::block::{AlgebraicBlock, Block};

/// Generic function block: outputs = f(inputs)
///
/// Uses const generics for fixed I/O sizes.
///
/// # Example
///
/// ```ignore
/// // Square function: y = x^2
/// let mut sq = Function::<1, 1, _>::new(|inputs, outputs| {
///     outputs[0] = inputs[0] * inputs[0];
/// });
///
/// // 2-input, 1-output: y = x1 * x2
/// let mut mult = Function::<2, 1, _>::new(|inputs, outputs| {
///     outputs[0] = inputs[0] * inputs[1];
/// });
/// ```
#[derive(Clone)]
pub struct Function<const IN: usize, const OUT: usize, F>
where
    F: Fn(&[f64; IN], &mut [f64; OUT]),
{
    inputs: [f64; IN],
    outputs: [f64; OUT],
    func: F,
}

impl<const IN: usize, const OUT: usize, F> Function<IN, OUT, F>
where
    F: Fn(&[f64; IN], &mut [f64; OUT]),
{
    /// Create function block with given function
    pub fn new(func: F) -> Self {
        Self {
            inputs: [0.0; IN],
            outputs: [0.0; OUT],
            func,
        }
    }
}

impl<const IN: usize, const OUT: usize, F> std::fmt::Debug for Function<IN, OUT, F>
where
    F: Fn(&[f64; IN], &mut [f64; OUT]),
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Function")
            .field("inputs", &self.inputs)
            .field("outputs", &self.outputs)
            .field("func", &"<closure>")
            .finish()
    }
}

impl<const IN: usize, const OUT: usize, F> Block for Function<IN, OUT, F>
where
    F: Fn(&[f64; IN], &mut [f64; OUT]),
{
    const NUM_INPUTS: usize = IN;
    const NUM_OUTPUTS: usize = OUT;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, _t: f64) {
        (self.func)(&self.inputs, &mut self.outputs);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
    }
}

impl<const IN: usize, const OUT: usize, F> AlgebraicBlock for Function<IN, OUT, F> where
    F: Fn(&[f64; IN], &mut [f64; OUT])
{
}

/// Time-dependent function block: outputs = f(t, inputs)
#[allow(dead_code)]
#[derive(Clone)]
pub struct TimeFunction<const IN: usize, const OUT: usize, F>
where
    F: Fn(f64, &[f64; IN], &mut [f64; OUT]),
{
    inputs: [f64; IN],
    outputs: [f64; OUT],
    func: F,
}

#[allow(dead_code)]
impl<const IN: usize, const OUT: usize, F> TimeFunction<IN, OUT, F>
where
    F: Fn(f64, &[f64; IN], &mut [f64; OUT]),
{
    pub fn new(func: F) -> Self {
        Self {
            inputs: [0.0; IN],
            outputs: [0.0; OUT],
            func,
        }
    }
}

impl<const IN: usize, const OUT: usize, F> std::fmt::Debug for TimeFunction<IN, OUT, F>
where
    F: Fn(f64, &[f64; IN], &mut [f64; OUT]),
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TimeFunction")
            .field("inputs", &self.inputs)
            .field("outputs", &self.outputs)
            .field("func", &"<closure>")
            .finish()
    }
}

impl<const IN: usize, const OUT: usize, F> Block for TimeFunction<IN, OUT, F>
where
    F: Fn(f64, &[f64; IN], &mut [f64; OUT]),
{
    const NUM_INPUTS: usize = IN;
    const NUM_OUTPUTS: usize = OUT;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, t: f64) {
        (self.func)(t, &self.inputs, &mut self.outputs);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_function_square() {
        let mut sq = Function::<1, 1, _>::new(|inputs, outputs| {
            outputs[0] = inputs[0] * inputs[0];
        });

        sq.inputs_mut()[0] = 3.0;
        sq.update(0.0);
        assert_eq!(sq.get_output(0), 9.0);
    }

    #[test]
    fn test_function_multiply() {
        let mut mult = Function::<2, 1, _>::new(|inputs, outputs| {
            outputs[0] = inputs[0] * inputs[1];
        });

        mult.inputs_mut()[0] = 3.0;
        mult.inputs_mut()[1] = 4.0;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 12.0);
    }

    #[test]
    fn test_function_multi_output() {
        let mut split = Function::<1, 2, _>::new(|inputs, outputs| {
            outputs[0] = inputs[0];
            outputs[1] = -inputs[0];
        });

        split.inputs_mut()[0] = 5.0;
        split.update(0.0);
        assert_eq!(split.outputs()[0], 5.0);
        assert_eq!(split.outputs()[1], -5.0);
    }

    #[test]
    fn test_time_function() {
        let mut tf = TimeFunction::<0, 1, _>::new(|t, _inputs, outputs| {
            outputs[0] = t * 2.0;
        });

        tf.update(3.0);
        assert_eq!(tf.get_output(0), 6.0);
    }
}
