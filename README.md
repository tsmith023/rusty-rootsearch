# rusty-rootsearch

This is a small Rust-package that performs Newton's method for finding the roots of a function using the method of bisections.

It stands completely on the shoulders of the [https://docs.rs/num-dual/latest/num_dual/](num_dual) crate that implements
auto-differentation for floating point numbers.

In an effort to generalise the implementation, the traits `Derivable` and `Coercable` are defined that allow the user to
define their own types that can be used within the `root_search`, `newton`, and `find_bisections` functions. The `Derivable` trait is used to
define the derivative of a number, which in the current implementation works well with the `num_dual` crate. The `Coercable`
trait is used to convert a floating point number into the type that is derivable. This is useful for example when using
the `f64` type, which is not derivable, but can be converted to a `Dual` type exposed by the `num_dual` crate.

## Example

```rust
use rusty_rootsearch::{root_search};
use num_dual::*;

fn find_sine_roots() {
    fn sine<D: DualNum<f32>>(x: D) -> D {
        x.sin()
    }
    let roots = root_search::<_,Dual32,f32>(&sine, -5.0, 5.0, 2000, 1000, 0.0001);
    for root in &roots.0 {
        println!("root: {}", root);
    }
    assert_eq!(roots.0.len(), 3);
    assert!(roots.0.contains(&std::f32::consts::PI));
    assert!(roots.0.contains(&(-std::f32::consts::PI)));
    assert!(roots.0.contains(&0.0));
}
```