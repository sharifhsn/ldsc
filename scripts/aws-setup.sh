#!/bin/bash
#
# One-time AWS infrastructure setup for ldsc benchmarking.
# Creates: ECR repo, S3 bucket, IAM roles, Batch compute environment + queue + job definition.
#
# Usage: ./scripts/aws-setup.sh
#
set -euo pipefail

PROFILE="AdministratorAccess-270497617191"
REGION="us-east-1"
ACCOUNT_ID="270497617191"
PREFIX="ldsc-bench"
S3_BUCKET="${PREFIX}-data-${ACCOUNT_ID}"

AWS="aws --profile $PROFILE --region $REGION --no-cli-pager"

echo "=== 1. ECR Repository ==="
$AWS ecr create-repository \
    --repository-name "$PREFIX" \
    --image-scanning-configuration scanOnPush=false \
    --image-tag-mutability MUTABLE \
    2>/dev/null && echo "Created ECR repo: $PREFIX" \
    || echo "ECR repo '$PREFIX' already exists"

echo ""
echo "=== 2. S3 Bucket for benchmark data ==="
$AWS s3 mb "s3://${S3_BUCKET}" 2>/dev/null \
    && echo "Created S3 bucket: $S3_BUCKET" \
    || echo "S3 bucket '$S3_BUCKET' already exists"

echo "Uploading benchmark data..."
for f in \
    data/1000G_phase3_common_norel.bed \
    data/1000G_phase3_common_norel.bim \
    data/1000G_phase3_common_norel.fam \
    data/bench_5k.bed \
    data/bench_5k.bim \
    data/bench_5k.fam; do
    BASENAME=$(basename "$f")
    echo "  $BASENAME..."
    $AWS s3 cp "$f" "s3://${S3_BUCKET}/${BASENAME}" --quiet
done
echo "Data upload complete."

echo ""
echo "=== 3. IAM Roles ==="

# ECS instance role (for EC2 instances in Batch compute environment)
cat > /tmp/ldsc-ec2-trust.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "ec2.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF

$AWS iam create-role \
    --role-name "${PREFIX}-instance-role" \
    --assume-role-policy-document file:///tmp/ldsc-ec2-trust.json \
    2>/dev/null && echo "Created instance role" \
    || echo "Instance role already exists"

$AWS iam attach-role-policy \
    --role-name "${PREFIX}-instance-role" \
    --policy-arn arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role \
    2>/dev/null || true

$AWS iam attach-role-policy \
    --role-name "${PREFIX}-instance-role" \
    --policy-arn arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess \
    2>/dev/null || true

echo "Attached policies to instance role"

# Instance profile
$AWS iam create-instance-profile \
    --instance-profile-name "${PREFIX}-instance-profile" \
    2>/dev/null && echo "Created instance profile" \
    || echo "Instance profile already exists"

$AWS iam add-role-to-instance-profile \
    --instance-profile-name "${PREFIX}-instance-profile" \
    --role-name "${PREFIX}-instance-role" \
    2>/dev/null || true

echo "Waiting 15s for IAM propagation..."
sleep 15

echo ""
echo "=== 4. Batch Compute Environment ==="

# Get default VPC networking
DEFAULT_VPC=$($AWS ec2 describe-vpcs --filters "Name=isDefault,Values=true" \
    --query "Vpcs[0].VpcId" --output text)
SUBNETS=$($AWS ec2 describe-subnets --filters "Name=vpc-id,Values=$DEFAULT_VPC" \
    --query "Subnets[*].SubnetId" --output json)
SG=$($AWS ec2 describe-security-groups \
    --filters "Name=vpc-id,Values=$DEFAULT_VPC" "Name=group-name,Values=default" \
    --query "SecurityGroups[0].GroupId" --output text)

echo "VPC: $DEFAULT_VPC, SG: $SG"

# Create compute environment config
cat > /tmp/ldsc-batch-ce.json << EOF
{
  "type": "SPOT",
  "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
  "minvCpus": 0,
  "maxvCpus": 64,
  "desiredvCpus": 0,
  "instanceTypes": ["c7a.4xlarge", "c7a.2xlarge", "c6a.4xlarge", "c6a.2xlarge"],
  "subnets": ${SUBNETS},
  "securityGroupIds": ["${SG}"],
  "instanceRole": "arn:aws:iam::${ACCOUNT_ID}:instance-profile/${PREFIX}-instance-profile"
}
EOF

$AWS batch create-compute-environment \
    --compute-environment-name "${PREFIX}-ce" \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file:///tmp/ldsc-batch-ce.json \
    2>/dev/null && echo "Created compute environment" \
    || echo "Compute environment already exists"

echo "Waiting for compute environment to become VALID..."
for i in $(seq 1 30); do
    STATUS=$($AWS batch describe-compute-environments \
        --compute-environments "${PREFIX}-ce" \
        --query "computeEnvironments[0].status" --output text 2>/dev/null || echo "UNKNOWN")
    if [ "$STATUS" = "VALID" ]; then
        echo "Compute environment is VALID"
        break
    fi
    echo "  Status: $STATUS (${i}/30)"
    sleep 10
done

echo ""
echo "=== 5. Batch Job Queue ==="
$AWS batch create-job-queue \
    --job-queue-name "${PREFIX}-queue" \
    --state ENABLED \
    --priority 1 \
    --compute-environment-order "order=1,computeEnvironment=${PREFIX}-ce" \
    2>/dev/null && echo "Created job queue" \
    || echo "Job queue already exists"

echo ""
echo "=== 6. CloudWatch Log Group ==="
$AWS logs create-log-group \
    --log-group-name "/aws/batch/${PREFIX}" \
    2>/dev/null && echo "Created log group" \
    || echo "Log group already exists"

# Set 7-day retention to minimize costs
$AWS logs put-retention-policy \
    --log-group-name "/aws/batch/${PREFIX}" \
    --retention-in-days 7

echo ""
echo "=== 7. Batch Job Definition ==="
cat > /tmp/ldsc-batch-jobdef.json << EOF
{
  "image": "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}:latest",
  "resourceRequirements": [
    {"type": "VCPU", "value": "16"},
    {"type": "MEMORY", "value": "8192"}
  ],
  "logConfiguration": {
    "logDriver": "awslogs",
    "options": {
      "awslogs-group": "/aws/batch/${PREFIX}",
      "awslogs-region": "${REGION}",
      "awslogs-stream-prefix": "bench"
    }
  }
}
EOF

$AWS batch register-job-definition \
    --job-definition-name "${PREFIX}-job" \
    --type container \
    --container-properties file:///tmp/ldsc-batch-jobdef.json

echo ""
echo "=== Setup Complete ==="
echo ""
echo "ECR:     ${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${PREFIX}"
echo "S3:      s3://${S3_BUCKET}/"
echo "Queue:   ${PREFIX}-queue"
echo "Job Def: ${PREFIX}-job"
echo ""
echo "Next: ./scripts/aws-bench.sh"
