// Forces the BLAS provider crate to be linked when BLAS features are enabled.
#[cfg(any(feature = "blas-openblas-static", feature = "blas-openblas-system"))]
#[allow(unused_imports)]
use blas_src as _;

#[cfg(any(feature = "blas-openblas-static", feature = "blas-openblas-system"))]
unsafe extern "C" {
    fn openblas_set_num_threads(num_threads: ::std::os::raw::c_int);
}

pub fn set_openblas_threads(num_threads: usize) {
    #[cfg(any(feature = "blas-openblas-static", feature = "blas-openblas-system"))]
    unsafe {
        openblas_set_num_threads(num_threads as i32);
    }
}
