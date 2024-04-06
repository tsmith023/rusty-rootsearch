use std::{env, fmt::Display, ops::{Sub, Div}};
// use std::{sync::mpsc::{Sender, Receiver, channel}, thread::{Thread,spawn, JoinHandle}};
use num_dual::{DualNumFloat,Dual32};

pub trait Derivable<T> where T: DualNumFloat {
    fn execute_derivative(&self) -> Self;
    fn zeroth_derivative(&self) -> T;
    fn first_derivative(&self) -> T;
}

pub trait Coerceable<T> where T: DualNumFloat{
    fn coerce_to(&self) -> T;
    fn coerce_from(value: T) -> Self;
}

impl Derivable<f32> for Dual32 {
    fn execute_derivative(&self) -> Self {
        return self.derivative()
    }
    fn zeroth_derivative(&self) -> f32 {
        return self.re
    }
    fn first_derivative(&self) -> f32 {
        return self.eps
    }
}

impl <T: DualNumFloat> Coerceable<T> for Dual32 {
    fn coerce_to(&self) -> T {
        return T::from(self.re).unwrap()
    }
    fn coerce_from(value: T) -> Self {
        return Dual32::from_re(value.to_f32().unwrap())
    }
}

pub struct NewtonOptions<T> where T: DualNumFloat {
    pub guess: T,
    pub patience: i32,
    pub tolerance: T
}

pub struct BisectionOptions<T> where T: DualNumFloat {
    pub lower: T,
    pub upper: T,
    pub resolution: i32
}

pub struct RootSearchOptions<T> where T: DualNumFloat {
    pub patience: i32,
    pub tolerance: T,
    pub lower: T,
    pub upper: T,
    pub resolution: i32
}

pub struct NewtonResult<T> where T: DualNumFloat {
    pub root: Option<T>,
    pub iterations: i32
}

pub struct BisectionResult<T> where T: DualNumFloat {
    pub lower: T,
    pub upper: T,
}

pub struct RootSearchResult<T> where T: DualNumFloat {
    pub roots: Vec<T>,
    pub bisections: Vec<BisectionResult<T>>,
}

fn newton<'a, F, N, T>(f: F, opts: NewtonOptions<T>) -> NewtonResult<T>
where
    F: Fn(N) -> N + Send + Sync + 'a,
    N: Derivable<T> + Coerceable<T> + Display + Copy,
    T: DualNumFloat
{
    let mut current: T = opts.guess;
    let mut count = 0;
    let debug = env::var("DEBUG").unwrap() == "true";
    loop {
        count += 1;
        let x = N::coerce_from(current).execute_derivative();
        let z = f(x);
        let next = x.zeroth_derivative() - z.zeroth_derivative() / z.first_derivative();
        let diff = next - current;
        if diff.abs() < opts.tolerance {
            if debug {
                println!("Found root at: {}", next);
            }
            return NewtonResult{
                root: Some(next),
                iterations: count
            };
        } else {
            if count > opts.patience {
                if debug {
                    println!("Failed to find root with initial guess of {}", opts.guess);
                    println!("Last iteration was: {}", current);
                    println!("Try updating the initial guess or increasing the tolerance or patience");
                }
                return NewtonResult{
                    root: None,
                    iterations: count
                };
            }
            current = next;
        }
    }
}

fn find_bisections<F, N, T>(f: F, opts: BisectionOptions<T>) -> Vec<BisectionResult<T>>
where
    F: Fn(N) -> N + Sync + Send + Copy,
    N: Derivable<T> + Coerceable<T> + Display + Copy + Sub + Div,
    T: DualNumFloat
{
    let step = (opts.upper - opts.lower) / T::from(opts.resolution).unwrap() + T::epsilon();
    // Add off-set to step to deal with roots at middle of lower and upper range
    let mut values: Vec<BisectionResult<T>> = Vec::new();

    for i in 0..opts.resolution {
        let a = opts.lower + step * T::from(i).unwrap();
        let b = opts.lower + step * T::from(i+1).unwrap();
        let fa = f(N::coerce_from(a));
        let fb = f(N::coerce_from(b));
        let pos2neg = fa.zeroth_derivative() > T::zero() && fb.zeroth_derivative() < T::zero();
        let neg2pos = fa.zeroth_derivative() < T::zero() && fb.zeroth_derivative() > T::zero();
        if pos2neg || neg2pos {
            values.push(BisectionResult{lower: a, upper: b});
        }
    };
    values
}

pub fn root_search<F, N, T>(f: F, opts: RootSearchOptions<T>) -> RootSearchResult<T>
where
    F: Fn(N) -> N + Sync + Send + Copy,
    N: Derivable<T> + Coerceable<T> + Display + Copy + Sub + Div,
    T: DualNumFloat
{
    if opts.lower > opts.upper {
        panic!("Lower bound must be greater than upper bound")
    }
    if opts.lower == opts.upper {
        panic!("Bounds cannot be the same")
    }
    let bisections = find_bisections(f, BisectionOptions{
        lower: opts.lower,
        upper: opts.upper,
        resolution: opts.resolution
    });
    let mut roots: Vec<T> = Vec::new();
    for bisection in &bisections {
        let res = T::from(100).unwrap();
        let step = (bisection.upper - bisection.lower) / res;
        for i in 0..res.to_i32().unwrap() {
            let guess = bisection.lower + (T::from(i).unwrap() * step);
            let res = newton(f, NewtonOptions{
                guess: guess,
                patience: opts.patience,
                tolerance: opts.tolerance
            });
            if res.root.is_none() {
                break;
            }
            let root = res.root.unwrap();
            if bisection.lower < root && root < bisection.upper {
                roots.push(root);
                break;
            }
        }

    }
    RootSearchResult{roots, bisections}
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_dual::{Dual32, DualNum};

    #[test]
    fn find_sine_root_newton() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let res = newton::<_,Dual32,f32>(&sine, NewtonOptions{
            guess: 2.0,
            patience: 1000,
            tolerance: 0.0001
        });
        assert_eq!(std::f32::consts::PI, res.root.unwrap())
    }

    #[test]
    fn find_cosine_root_newton() {
        fn cosine<D: DualNum<f32>>(x: D) -> D {
            x.cos()
        }
        let res = newton::<_,Dual32,f32>(&cosine, NewtonOptions{
            guess: 2.0,
            patience: 1000,
            tolerance: 0.0001
        });
        assert_eq!(std::f32::consts::PI / 2.0, res.root.unwrap())
    }

    #[test]
    fn find_sine_bisections() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let bisections = find_bisections::<_,Dual32,f32>(&sine, BisectionOptions{
            lower: -5.0, 
            upper: 5.0, 
            resolution: 1000
        });
        for bisection in &bisections {
            println!("bisection: ({},{})", bisection.lower, bisection.upper)
        }
        assert_eq!(bisections.len(), 3)
    }

    #[test]
    fn find_cosine_bisections() {
        fn cosine<D: DualNum<f32>>(x: D) -> D {
            x.cos()
        }
        let bisections = find_bisections::<_,Dual32,f32>(&cosine, BisectionOptions{
            lower: -5.0, 
            upper: 5.0, 
            resolution: 1000
        });
        for bisection in &bisections {
            println!("bisection: ({},{})", bisection.lower, bisection.upper)
        }
        assert_eq!(bisections.len(), 4)
    }

    #[test]
    fn find_sine_roots() {
        fn sine<D: DualNum<f32>>(x: D) -> D {
            x.sin()
        }
        let res = root_search::<_,Dual32,f32>(&sine, RootSearchOptions{
            lower: -5.0,
            upper: 5.0,
            patience: 2000,
            tolerance: 0.0001,
            resolution: 1000
        });
        for root in &res.roots {
            println!("root: {}", root);
        }
        assert_eq!(res.roots.len(), 3);
        assert!(res.roots.contains(&std::f32::consts::PI));
        assert!(res.roots.contains(&(-std::f32::consts::PI)));
        assert!(res.roots.contains(&0.0));
    }

    #[test]
    fn find_cosine_roots() {
        fn cosine<D: DualNum<f32>>(x: D) -> D {
            x.cos()
        }
        let res = root_search::<_,Dual32,f32>(&cosine, RootSearchOptions{
            lower: -5.0,
            upper: 5.0,
            patience: 2000,
            tolerance: 0.0001,
            resolution: 1000
        });
        for root in &res.roots {
            println!("root: {}", root);
        }
        assert_eq!(res.roots.len(), 4);
        assert!(res.roots.contains(&std::f32::consts::FRAC_PI_2));
        assert!(res.roots.contains(&(-std::f32::consts::FRAC_PI_2)));
        assert!(res.roots.contains(&(std::f32::consts::FRAC_PI_2 * 3.0)));
        assert!(res.roots.contains(&(-std::f32::consts::FRAC_PI_2 * 3.0)));
    }


}