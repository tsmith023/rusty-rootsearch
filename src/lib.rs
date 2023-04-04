use std;
use num_dual::*;

fn get_dual(x: f32) -> DualVec<f32, f32, 1> {
    Dual32::from_re(x)
}

fn newton<'a, F>(f: F, guess: f32, patience: i32, tolerance: f32) -> Option<f32> 
where
    F: Fn(DualVec<f32, f32, 1>) -> DualVec<f32, f32, 1> + Send + Sync + 'a
{
    let mut current: f32 = guess;
    let mut count = 0;
    loop {
        count += 1;
        let x_dual = Dual32::from_re(current).derive();
        let z_dual = f(x_dual);
        let next = x_dual.re - z_dual.re / z_dual.eps[0];
        let diff = next - current;
        if diff.abs() < tolerance {
            println!("Found route at: {}", next);
            return Some(next);
        } else {
            if count > patience {
                println!("Failed to find root with initial guess of {}", guess);
                println!("Last iteration was: {}", current);
                println!("Try updating the initial guess or increasing the tolerance or patience");
                return None;
            }
            current = next;
        }
    }
}

fn find_bisections<'a, F>(f: F, lower: f32, upper: f32, resolution: i32) -> Vec<(f32, f32)>
where
    F: Fn(DualVec<f32, f32, 1>) -> DualVec<f32, f32, 1> + Sync + Send + Copy + 'a
{
    let step = (upper - lower) / (resolution as f32 + std::f32::consts::PI);
    // Add off-set to resolution to deal with roots at middle of lower and upper range
    let mut values: Vec<(f32, f32)> = Vec::new();

    for i in 0..resolution {
        let a = lower + step * i as f32;
        let b = lower + step * (i + 1) as f32;
        let fa = f(get_dual(a));
        let fb = f(get_dual(b));
        let pos2neg = fa.re > 0.0 && fb.re < 0.0;
        let neg2pos = fa.re < 0.0 && fb.re > 0.0;
        if pos2neg || neg2pos {
            values.push((a, b));
        }
    };
    values
}

pub fn root_search<'a, F>(f: F, lower: f32, upper: f32, resolution: i32, patience: i32, tolerance: f32) -> (Vec<f32>, Vec<(f32, f32)>)
where
    F: Fn(DualVec<f32, f32, 1>) -> DualVec<f32, f32, 1> + Sync + Send + Copy + 'a
{
    if lower > upper {
        panic!("Lower bound must be greater than upper bound")
    }
    if lower == upper {
        panic!("Bounds cannot be the same")
    }
    let bisections = find_bisections(f, lower, upper, resolution);
    // TODO: Launch threading here to find the the root contained within each values pair,
    // which are the lower and upper bounds (a_dual, b_dual) of a suspected root.
    let mut roots: Vec<f32> = Vec::new();
    for bisection in &bisections {
        let res = 100;
        let step = (bisection.1 - bisection.0) / (res as f32);
        for i in 0..res {
            let guess = bisection.0 + (i as f32 * step);
            let root = newton(f, guess, patience, tolerance);
            if root.is_none() {
                break;
            }
            let root = root.unwrap();
            if bisection.0 < root && root < bisection.1 {
                roots.push(root);
                break;
            } else if bisection.0 < root && root < bisection.1 {
                roots.push(root);
                break;
            }
        }

    }
    (roots, bisections)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_sine_root_newton() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let root = newton(&sine, 2.0, 1000, 0.0001);
        assert_eq!(std::f32::consts::PI, root.unwrap())
    }

    #[test]
    fn find_sine_bisections() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let bisections = find_bisections(&sine, -5.0, 5.0, 2000);
        for bisection in &bisections {
            println!("bisection: ({},{})", bisection.0, bisection.1)
        }
        assert_eq!(bisections.len(), 3)
    }

    #[test]
    fn find_sine_roots() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let roots = root_search(&sine, -5.0, 5.0, 2000, 1000, 0.0001);
        for root in &roots.0 {
            println!("root: {}", root);
        }
        assert_eq!(roots.0.len(), 3)
    }

    #[test]
    fn find_cosine_roots() {
        fn cosine<D: DualNum<f32>>(x: D) -> D {
            x.cos()
        }
        let roots = root_search(&cosine, -5.0, 5.0, 2000, 1000, 0.0001);
        for root in &roots.0 {
            println!("root: {}", root);
        }
        assert_eq!(roots.0.len(), 4)
    }


}