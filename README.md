# rusty-search

This is a small Rust-package that performs Newton's method for finding the roots of a function using the method of bisections.

It stands completely on the shoulders of the [https://docs.rs/num-dual/latest/num_dual/](num_dual) crate that implements
auto-differentation for floating point numbers.

This package is currently heavily coupled to the underlying implementation of
DualVec32. However, this will be abstracted away into generics in future interations.
