# AWS Serverless Benchmarking Pipeline Spec

## 1. Project Context & Objective

* **Application:** A highly optimized Rust rewrite of the `ldsc` (LD score regression) library.
* **Current State:** The application is containerized. However, local benchmarking via `hyperfine` on the current workstation is hitting thermal throttling limits, corrupting the timing distributions.
* **Target Deployment:** Bare-metal HPC.
* **Goal:** Build an automated, cost-effective AWS compute pipeline that closely replicates high-core HPC performance to execute isolated, reproducible benchmarking runs.
* **Execution Environment:** The deployment scripts will be run from a macOS terminal, interacting primarily with the AWS CLI or `boto3`.

## 2. Infrastructure Architecture Options & Pricing

The pipeline must utilize AWS Batch to queue and execute the containerized jobs. Depending on the scale of the benchmarking matrix, configure the AWS Batch Compute Environment for one of two paths:

* **Path A: Maximum Convenience (AWS Fargate Spot)**
* **Limits:** Maximum of 16 vCPUs and 120 GB RAM per task.
* **Cost Estimate (us-east-1):** Fargate Spot rates are roughly $0.013 per vCPU/hour and $0.0014 per GB/hour. A maxed-out 16 vCPU / 64 GB task costs approximately **$0.30 per hour**.
* **Trade-offs:** Zero infrastructure management and near-instant startup, but it has a lower compute ceiling and slightly higher hypervisor variance.


* **Path B: HPC-Grade Power (AWS Batch with EC2 Spot)**
* **Target Hardware:** Compute-optimized AMD EPYC instances (`c7a` family) or Intel (`c7i`). For example, a `c7a.24xlarge` provides 96 vCPUs.
* **Cost Estimate (us-east-1):** On-demand pricing is roughly $4.40 per hour. Utilizing EC2 Spot instances typically reduces this by 60% to 70%, bringing the cost to approximately **$1.30 to $1.50 per hour**.
* **Trade-offs:** Delivers true HPC-tier compute with dedicated physical cores, but the Compute Environment requires an extra 1-2 minutes to provision the underlying EC2 instance before the container starts.



## 3. Agent Implementation Tasks

Please generate the code and infrastructure configuration to fulfill the following three stages:

### A. Dockerfile Enhancements

* Ensure `hyperfine` is installed in the final container layer.
* Set the `ENTRYPOINT` to execute the benchmarking matrix.
* **Variance Mitigation:** Because virtualized cloud environments introduce microarchitectural noise, configure the `hyperfine` commands to run aggressively (e.g., `--runs 20` or higher) and utilize the `--warmup` flag to establish a stable mean and reduce standard error.

### B. Deployment Automation Script

Write a unified Bash or Python (`boto3`) script that performs the entire loop:

1. Authenticates with AWS ECR.
2. Builds the updated Docker image locally and pushes it to an ECR repository.
3. Registers a new AWS Batch Job Definition that points to the latest ECR image tag.
4. Submits the Job to the AWS Batch Job Queue.

### C. Log Streaming & Result Retrieval

* The script must not terminate after job submission. It should poll the AWS Batch job status, and once the container is in the `RUNNING` state, automatically fetch and stream the stdout from AWS CloudWatch Logs back to the local terminal so the `hyperfine` results can be viewed in real-time.
