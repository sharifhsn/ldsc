#![cfg(feature = "gpu")]

use anyhow::{Result, anyhow};
use cubecl::client::ComputeClient;
use cubecl::cuda::{CudaDevice, CudaRuntime};
use cubecl::prelude::*;
use cubecl::std::tensor::TensorHandle;
use cubek_matmul::definition::MatmulElems;
use cubek_matmul::launch::{MatmulInputHandleRef, Strategy};

type R = CudaRuntime;

/// GPU context holding a CUDA device and compute client.
pub struct GpuContext {
    client: ComputeClient<R>,
    strategy: Strategy,
}

/// Extract a contiguous row-major (m × n) result from a potentially pitched GPU output buffer.
fn extract_result(raw: &[f32], m: usize, n: usize, row_stride: usize) -> Vec<f32> {
    let mut out = Vec::with_capacity(m * n);
    for i in 0..m {
        let start = i * row_stride;
        out.extend_from_slice(&raw[start..start + n]);
    }
    out
}

impl GpuContext {
    /// Create a new GPU context on CUDA device 0.
    pub fn new() -> Result<Self> {
        let device = CudaDevice::new(0);
        let client = R::client(&device);
        Ok(Self {
            client,
            strategy: Strategy::Auto,
        })
    }

    /// Perform A^T × B on GPU.
    ///
    /// `a_data` is column-major f32 (n_rows × a_cols), `b_data` is column-major f32 (n_rows × b_cols).
    /// Returns the result matrix (a_cols × b_cols) as a flat row-major f32 vec.
    pub fn matmul_tn(
        &self,
        a_data: &[f32],
        n_rows: usize,
        a_cols: usize,
        b_data: &[f32],
        b_cols: usize,
    ) -> Result<Vec<f32>> {
        let m = a_cols;
        let k = n_rows;
        let n = b_cols;
        let f32_storage = f32::as_type_native_unchecked();

        // Upload A and B to GPU
        let a_handle = self.client.create_from_slice(f32::as_bytes(a_data));
        let b_handle = self.client.create_from_slice(f32::as_bytes(b_data));

        // A is column-major (k × m): element (i,j) at offset j*k + i.
        // For A^T (m × k): element (j,i) maps to j*k + i, so shape [m,k], strides [k,1].
        let a_tensor = TensorHandle::<R>::new(a_handle, vec![m, k], vec![k, 1], f32_storage);

        // B is column-major (k × n): element (i,j) at offset j*k + i, strides [1, k].
        let b_tensor = TensorHandle::<R>::new(b_handle, vec![k, n], vec![1, k], f32_storage);

        // Output: (m × n) — may be pitched/padded by the runtime
        let out_tensor = TensorHandle::<R>::empty(&self.client, vec![m, n], f32_storage);
        let out_row_stride = out_tensor.strides()[0];

        let mut dtypes = MatmulElems::from_single_dtype(f32_storage);

        cubek_matmul::launch::launch_ref::<R>(
            &self.strategy,
            &self.client,
            &MatmulInputHandleRef::new(a_tensor.as_ref(), f32_storage),
            &MatmulInputHandleRef::new(b_tensor.as_ref(), f32_storage),
            &out_tensor.as_ref(),
            &mut dtypes,
        )
        .map_err(|e| anyhow!("GPU matmul failed: {:?}", e))?;

        // Read result back, accounting for potential row padding
        let out_bytes = self.client.read_one(out_tensor.handle);
        let raw: &[f32] = f32::from_bytes(&out_bytes);
        Ok(extract_result(raw, m, n, out_row_stride))
    }

    /// Perform A^T × B on GPU with tiled A for large windows that don't fit VRAM.
    #[allow(dead_code)]
    pub fn matmul_tn_tiled(
        &self,
        a_data: &[f32],
        n_rows: usize,
        a_cols: usize,
        b_data: &[f32],
        b_cols: usize,
        tile_cols: usize,
    ) -> Result<Vec<f32>> {
        let n = b_cols;
        let f32_storage = f32::as_type_native_unchecked();
        let mut result = vec![0.0f32; a_cols * n];

        for tile_start in (0..a_cols).step_by(tile_cols) {
            let tile_end = (tile_start + tile_cols).min(a_cols);
            let tile_m = tile_end - tile_start;

            // Extract tile columns from column-major A
            let tile_data: Vec<f32> = (tile_start..tile_end)
                .flat_map(|col| {
                    let offset = col * n_rows;
                    a_data[offset..offset + n_rows].iter().copied()
                })
                .collect();

            // Re-upload B each tile (handle is consumed by read_one)
            let b_handle = self.client.create_from_slice(f32::as_bytes(b_data));
            let b_tensor =
                TensorHandle::<R>::new(b_handle, vec![n_rows, n], vec![1, n_rows], f32_storage);

            let a_handle = self.client.create_from_slice(f32::as_bytes(&tile_data));
            let a_tensor = TensorHandle::<R>::new(
                a_handle,
                vec![tile_m, n_rows],
                vec![n_rows, 1],
                f32_storage,
            );

            let out_tensor = TensorHandle::<R>::empty(&self.client, vec![tile_m, n], f32_storage);
            let out_row_stride = out_tensor.strides()[0];

            let mut dtypes = MatmulElems::from_single_dtype(f32_storage);

            cubek_matmul::launch::launch_ref::<R>(
                &self.strategy,
                &self.client,
                &MatmulInputHandleRef::new(a_tensor.as_ref(), f32_storage),
                &MatmulInputHandleRef::new(b_tensor.as_ref(), f32_storage),
                &out_tensor.as_ref(),
                &mut dtypes,
            )
            .map_err(|e| anyhow!("GPU matmul tile failed: {:?}", e))?;

            let tile_bytes = self.client.read_one(out_tensor.handle);
            let raw: &[f32] = f32::from_bytes(&tile_bytes);

            for i in 0..tile_m {
                for j in 0..n {
                    result[(tile_start + i) * n + j] = raw[i * out_row_stride + j];
                }
            }
        }

        Ok(result)
    }

    /// Perform A^T × B on GPU using flex32 precision (f32 upload, f16 in GPU shared memory,
    /// f32 accumulation and output). Zero CPU-side conversion overhead.
    pub fn matmul_tn_flex32(
        &self,
        a_data: &[f32],
        n_rows: usize,
        a_cols: usize,
        b_data: &[f32],
        b_cols: usize,
    ) -> Result<Vec<f32>> {
        let m = a_cols;
        let k = n_rows;
        let n = b_cols;
        let f32_storage = f32::as_type_native_unchecked();

        // Upload f32 directly — GPU converts to f16 in shared memory
        let a_handle = self.client.create_from_slice(f32::as_bytes(a_data));
        let b_handle = self.client.create_from_slice(f32::as_bytes(b_data));

        let a_tensor = TensorHandle::<R>::new(a_handle, vec![m, k], vec![k, 1], f32_storage);
        let b_tensor = TensorHandle::<R>::new(b_handle, vec![k, n], vec![1, k], f32_storage);

        let out_tensor = TensorHandle::<R>::empty(&self.client, vec![m, n], f32_storage);
        let out_row_stride = out_tensor.strides()[0];

        // flex32: global f32, stage/register f16, accumulator f32
        let mut dtypes = MatmulElems::new_deprecated::<flex32>();

        cubek_matmul::launch::launch_ref::<R>(
            &self.strategy,
            &self.client,
            &MatmulInputHandleRef::new(a_tensor.as_ref(), f32_storage),
            &MatmulInputHandleRef::new(b_tensor.as_ref(), f32_storage),
            &out_tensor.as_ref(),
            &mut dtypes,
        )
        .map_err(|e| anyhow!("GPU flex32 matmul failed: {:?}", e))?;

        let out_bytes = self.client.read_one(out_tensor.handle);
        let raw: &[f32] = f32::from_bytes(&out_bytes);
        Ok(extract_result(raw, m, n, out_row_stride))
    }

    /// Perform A^T × B on GPU using flex32 precision with tiled A.
    #[allow(dead_code)]
    pub fn matmul_tn_tiled_flex32(
        &self,
        a_data: &[f32],
        n_rows: usize,
        a_cols: usize,
        b_data: &[f32],
        b_cols: usize,
        tile_cols: usize,
    ) -> Result<Vec<f32>> {
        let n = b_cols;
        let f32_storage = f32::as_type_native_unchecked();
        let mut result = vec![0.0f32; a_cols * n];

        for tile_start in (0..a_cols).step_by(tile_cols) {
            let tile_end = (tile_start + tile_cols).min(a_cols);
            let tile_m = tile_end - tile_start;

            let tile_data: Vec<f32> = (tile_start..tile_end)
                .flat_map(|col| {
                    let offset = col * n_rows;
                    a_data[offset..offset + n_rows].iter().copied()
                })
                .collect();

            let b_handle = self.client.create_from_slice(f32::as_bytes(b_data));
            let b_tensor =
                TensorHandle::<R>::new(b_handle, vec![n_rows, n], vec![1, n_rows], f32_storage);

            let a_handle = self.client.create_from_slice(f32::as_bytes(&tile_data));
            let a_tensor = TensorHandle::<R>::new(
                a_handle,
                vec![tile_m, n_rows],
                vec![n_rows, 1],
                f32_storage,
            );

            let out_tensor = TensorHandle::<R>::empty(&self.client, vec![tile_m, n], f32_storage);
            let out_row_stride = out_tensor.strides()[0];

            let mut dtypes = MatmulElems::new_deprecated::<flex32>();

            cubek_matmul::launch::launch_ref::<R>(
                &self.strategy,
                &self.client,
                &MatmulInputHandleRef::new(a_tensor.as_ref(), f32_storage),
                &MatmulInputHandleRef::new(b_tensor.as_ref(), f32_storage),
                &out_tensor.as_ref(),
                &mut dtypes,
            )
            .map_err(|e| anyhow!("GPU flex32 matmul tile failed: {:?}", e))?;

            let tile_bytes = self.client.read_one(out_tensor.handle);
            let raw: &[f32] = f32::from_bytes(&tile_bytes);

            for i in 0..tile_m {
                for j in 0..n {
                    result[(tile_start + i) * n + j] = raw[i * out_row_stride + j];
                }
            }
        }

        Ok(result)
    }
}
