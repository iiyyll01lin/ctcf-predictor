# Quick PWM Testing Script
# This script runs a subset of PWM improvements for quick validation

echo "========================================"
echo "Quick PWM Improvement Test"
echo "========================================"

# Test 1: Simple Aligned PWM
echo "Test 1: Simple Aligned PWM"
echo "-------------------------"
./run-in-docker.sh Rscript scripts/simple_aligned_pwm.R data/aligned_sequences.fasta results/test_simple_aligned.rds 0.1

if [ $? -eq 0 ]; then
    echo "✓ Simple aligned PWM test passed"
else
    echo "✗ Simple aligned PWM test failed"
fi

echo ""

# Test 2: Subset PWM (small size for quick test)
echo "Test 2: High-Quality Subset PWM"
echo "------------------------------"
./run-in-docker.sh Rscript scripts/build_subset_pwm.R data/training_sequences.fasta results/test_subset 1000 0.01

if [ $? -eq 0 ]; then
    echo "✓ Subset PWM test passed"
else
    echo "✗ Subset PWM test failed"
fi

echo ""

# Test 3: PWM Comparison
echo "Test 3: PWM Comparison"
echo "---------------------"
./run-in-docker.sh Rscript scripts/compare_pwms.R results results/test_comparison.html

if [ $? -eq 0 ]; then
    echo "✓ PWM comparison test passed"
else
    echo "✗ PWM comparison test failed"
fi

echo ""
echo "Quick test completed!"
echo "Check results/ directory for generated files"
