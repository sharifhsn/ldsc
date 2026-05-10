fn main() {
    println!("cargo:rerun-if-env-changed=CUDARC_CUDA_VERSION");
    let version = std::env::var("CUDARC_CUDA_VERSION").unwrap_or_else(|_| "auto".to_string());
    println!("cargo:rustc-env=LDSC_CUDA_VERSION={version}");
}
