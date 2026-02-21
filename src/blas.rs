// Forces the BLAS provider crate to be linked when BLAS features are enabled.
#[cfg(any(feature = "blas-openblas-static", feature = "blas-openblas-system"))]
#[allow(unused_imports)]
use blas_src as _;
