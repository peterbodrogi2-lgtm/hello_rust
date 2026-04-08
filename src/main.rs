use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::goldensectionsearch::GoldenSectionSearch;

struct MyFunc;

impl CostFunction for MyFunc {
    type Param = f64;
    type Output = f64;

    fn cost(&self, x: &f64) -> Result<f64, Error> {
        Ok((x - 2.0).powi(2) + 1.0)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let problem = MyFunc;

    let solver = GoldenSectionSearch::new(0.0_f64, 5.0_f64)?
        .with_tolerance(1e-8)?;

    let init = 2.5_f64;

    let result = Executor::new(problem, solver)
        .configure(|state| state.param(init).max_iters(100))
        .run()?;

    let state = result.state();
    let x_min = state.get_best_param().expect("no parameter returned");
    let f_min = state.get_best_cost();

    println!("x_min = {}", x_min);
    println!("f(x_min) = {}", f_min);

    Ok(())
}