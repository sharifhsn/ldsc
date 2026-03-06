use crate::la::{ColF, MatF};

#[derive(Debug)]
#[allow(dead_code)]
pub struct IrwlsResult {
    pub est: ColF,
    pub jknife_se: Option<ColF>,
    pub jknife_var: Option<ColF>,
    pub jknife_cov: Option<MatF>,
    pub delete_values: Option<MatF>,
}

