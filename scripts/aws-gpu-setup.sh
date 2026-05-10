#!/bin/bash
#
# One-time AWS infrastructure setup for GPU benchmarking.
# Creates: GPU compute environment (g5.2xlarge SPOT), job queue, job definition.
# Reuses existing ECR repo, S3 bucket, IAM roles, and CloudWatch log group.
#
# Usage: ./scripts/aws-gpu-setup.sh
#
set -euo pipefail

PROFILE="AdministratorAccess-270497617191"
REGION="us-east-1"
ACCOUNT_ID="270497617191"
PREFIX="ldsc-bench"

AWS="aws --profile $PROFILE --region $REGION --no-cli-pager"

# Reuse existing networking from the CPU setup
DEFAULT_VPC=$($AWS ec2 describe-vpcs --filters "Name=isDefault,Values=true" \
    --query "Vpcs[0].VpcId" --output text)
SUBNETS=$($AWS ec2 describe-subnets --filters "Name=vpc-id,Values=$DEFAULT_VPC" \
    --query "Subnets[*].SubnetId" --output json)
SG=$($AWS ec2 describe-security-groups \
    --filters "Name=vpc-id,Values=$DEFAULT_VPC" "Name=group-name,Values=default" \
    --query "SecurityGroups[0].GroupId" --output text)

echo "VPC: $DEFAULT_VPC, SG: $SG"

echo ""
echo "=== 1. GPU Compute Environment ==="

# Reuse the existing launch template for 50GB gp3 EBS
LT_ID="lt-05105c54819b20bfb"

cat > /tmp/ldsc-gpu-ce.json << EOF
{
  "type": "SPOT",
  "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
  "minvCpus": 0,
  "maxvCpus": 16,
  "desiredvCpus": 0,
  "instanceTypes": ["g5.2xlarge", "g5.4xlarge"],
  "subnets": ${SUBNETS},
  "securityGroupIds": ["${SG}"],
  "instanceRole": "arn:aws:iam::${ACCOUNT_ID}:instance-profile/${PREFIX}-instance-profile",
  "ec2Configuration": [{"imageType": "ECS_AL2_NVIDIA"}],
  "launchTemplate": {"launchTemplateId": "${LT_ID}"}
}
EOF

$AWS batch create-compute-environment \
    --compute-environment-name "${PREFIX}-gpu-ce" \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file:///tmp/ldsc-gpu-ce.json \
    2>/dev/null && echo "Created GPU compute environment" \
    || echo "GPU compute environment already exists"

echo "Waiting for GPU compute environment to become VALID..."
for i in $(seq 1 30); do
    STATUS=$($AWS batch describe-compute-environments \
        --compute-environments "${PREFIX}-gpu-ce" \
        --query "computeEnvironments[0].status" --output text 2>/dev/null || echo "UNKNOWN")
    if [ "$STATUS" = "VALID" ]; then
        echo "GPU compute environment is VALID"
        break
    fi
    echo "  Status: $STATUS (${i}/30)"
    sleep 10
done

echo ""
echo "=== 2. GPU Job Queue ==="
$AWS batch create-job-queue \
    --job-queue-name "${PREFIX}-gpu-queue" \
    --state ENABLED \
    --priority 1 \
    --compute-environment-order "order=1,computeEnvironment=${PREFIX}-gpu-ce" \
    2>/dev/null && echo "Created GPU job queue" \
    || echo "GPU job queue already exists"

echo ""
echo "=== 3. GPU Job Definition ==="
cat > /tmp/ldsc-gpu-jobdef.json << EOF
{
  "image": "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}:gpu-latest",
  "resourceRequirements": [
    {"type": "VCPU", "value": "8"},
    {"type": "MEMORY", "value": "28672"},
    {"type": "GPU", "value": "1"}
  ],
  "logConfiguration": {
    "logDriver": "awslogs",
    "options": {
      "awslogs-group": "/aws/batch/${PREFIX}",
      "awslogs-region": "${REGION}",
      "awslogs-stream-prefix": "gpu-bench"
    }
  }
}
EOF

$AWS batch register-job-definition \
    --job-definition-name "${PREFIX}-gpu-job" \
    --type container \
    --container-properties file:///tmp/ldsc-gpu-jobdef.json

echo ""
echo "=== GPU Setup Complete ==="
echo ""
echo "Queue:   ${PREFIX}-gpu-queue"
echo "Job Def: ${PREFIX}-gpu-job"
echo "Image:   ${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}:gpu-latest"
echo ""
echo "Next: ./scripts/aws-gpu-bench.sh"
