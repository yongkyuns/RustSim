#!/usr/bin/env bash

# Benchmark runner script for RustSim vs PathSim (Python) comparison
#
# This script:
# 1. Compiles the RustSim benchmarks in release mode
# 2. Runs both Python and Rust benchmarks
# 3. Collects and compares the results
# 4. Generates a summary report

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Banner
print_banner() {
    echo -e "${BOLD}${BLUE}"
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                                                                ║"
    echo "║        RustSim vs PathSim (Python) Benchmark Suite            ║"
    echo "║                                                                ║"
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

# Print section header
print_section() {
    echo -e "\n${BOLD}${BLUE}==== $1 ====${NC}\n"
}

# Print error message and exit
error_exit() {
    echo -e "${RED}ERROR: $1${NC}" >&2
    exit 1
}

# Print warning message
warn() {
    echo -e "${YELLOW}WARNING: $1${NC}" >&2
}

# Print success message
success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check dependencies
check_dependencies() {
    print_section "Checking Dependencies"

    local all_ok=true

    # Check for Rust/Cargo
    if command_exists cargo; then
        local rust_version=$(cargo --version)
        success "Cargo found: $rust_version"
    else
        error_exit "cargo not found. Please install Rust from https://rustup.rs/"
    fi

    # Check for Python 3
    if command_exists python3; then
        local python_version=$(python3 --version)
        success "Python found: $python_version"
    else
        error_exit "python3 not found. Please install Python 3.8 or later"
    fi

    # Check for jq (optional, for JSON parsing)
    if command_exists jq; then
        success "jq found (for JSON parsing)"
    else
        warn "jq not found. Install for better JSON output (optional)"
    fi

    echo ""
}

# Compile Rust benchmarks
compile_rust_benchmarks() {
    print_section "Compiling Rust Benchmarks"

    cd "$PROJECT_ROOT"

    echo "Building in release mode with optimizations..."
    echo "  - LTO enabled"
    echo "  - Single codegen unit"
    echo "  - Optimization level 3"
    echo ""

    # Build the benchmark binary
    if cargo build --release --bin harmonic_oscillator_rust 2>&1 | grep -v "Compiling\|Finished"; then
        success "Rust benchmark compiled successfully"
    else
        # Try adding the binary to Cargo.toml if it doesn't exist
        if ! grep -q "harmonic_oscillator_rust" "$PROJECT_ROOT/Cargo.toml"; then
            echo ""
            echo "Adding benchmark binary to Cargo.toml..."

            # Create a temporary benchmark binary entry
            cat >> "$PROJECT_ROOT/Cargo.toml" << 'EOF'

[[bin]]
name = "harmonic_oscillator_rust"
path = "comparison/harmonic_oscillator_rust.rs"
EOF
            success "Added benchmark binary to Cargo.toml"

            # Try building again
            echo "Retrying compilation..."
            cargo build --release --bin harmonic_oscillator_rust || error_exit "Failed to compile Rust benchmark"
            success "Rust benchmark compiled successfully"
        else
            error_exit "Failed to compile Rust benchmark"
        fi
    fi

    echo ""
}

# Run Python benchmark
run_python_benchmark() {
    print_section "Running Python Benchmark"

    cd "$SCRIPT_DIR"

    echo "Executing Python harmonic oscillator simulation..."
    echo ""

    # Make Python script executable
    chmod +x harmonic_oscillator_python.py

    # Run and capture output
    python3 harmonic_oscillator_python.py --json | tee /tmp/python_output.txt

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        success "Python benchmark completed successfully"
    else
        error_exit "Python benchmark failed"
    fi

    echo ""
}

# Run Rust benchmark
run_rust_benchmark() {
    print_section "Running Rust Benchmark"

    cd "$PROJECT_ROOT"

    echo "Executing Rust harmonic oscillator simulation..."
    echo ""

    # Run and capture output
    ./target/release/harmonic_oscillator_rust --json | tee /tmp/rust_output.txt

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        success "Rust benchmark completed successfully"
    else
        error_exit "Rust benchmark failed"
    fi

    echo ""
}

# Compare results
compare_results() {
    print_section "Benchmark Comparison"

    # Check if JSON files exist
    if [ ! -f /tmp/python_benchmark_results.json ] || [ ! -f /tmp/rust_benchmark_results.json ]; then
        warn "JSON result files not found. Skipping detailed comparison."
        return
    fi

    # Extract results using jq if available, otherwise use grep/sed
    if command_exists jq; then
        local python_time=$(jq -r '.elapsed_time' /tmp/python_benchmark_results.json)
        local rust_time=$(jq -r '.elapsed_time' /tmp/rust_benchmark_results.json)
        local python_steps=$(jq -r '.steps_per_second' /tmp/python_benchmark_results.json)
        local rust_steps=$(jq -r '.steps_per_second' /tmp/rust_benchmark_results.json)
        local python_pos_err=$(jq -r '.errors.position' /tmp/python_benchmark_results.json)
        local rust_pos_err=$(jq -r '.errors.position' /tmp/rust_benchmark_results.json)
        local python_vel_err=$(jq -r '.errors.velocity' /tmp/python_benchmark_results.json)
        local rust_vel_err=$(jq -r '.errors.velocity' /tmp/rust_benchmark_results.json)
        local python_energy_err=$(jq -r '.errors.energy' /tmp/python_benchmark_results.json)
        local rust_energy_err=$(jq -r '.errors.energy' /tmp/rust_benchmark_results.json)
    else
        local python_time=$(grep "elapsed_time" /tmp/python_benchmark_results.json | sed 's/.*: \([0-9.]*\).*/\1/')
        local rust_time=$(grep "elapsed_time" /tmp/rust_benchmark_results.json | sed 's/.*: \([0-9.]*\).*/\1/')
        local python_steps="N/A"
        local rust_steps="N/A"
        local python_pos_err="N/A"
        local rust_pos_err="N/A"
        local python_vel_err="N/A"
        local rust_vel_err="N/A"
        local python_energy_err="N/A"
        local rust_energy_err="N/A"
    fi

    # Calculate speedup
    local speedup=$(python3 -c "print(f'{${python_time}/${rust_time}:.2f}')" 2>/dev/null || echo "N/A")

    # Print comparison table
    echo -e "${BOLD}Performance Comparison:${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    printf "%-25s %15s %15s\n" "Metric" "Python" "Rust"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    printf "%-25s %15.6f %15.6f\n" "Execution Time (s)" "$python_time" "$rust_time"

    if [ "$python_steps" != "N/A" ]; then
        printf "%-25s %15.0f %15.0f\n" "Steps/Second" "$python_steps" "$rust_steps"
    fi

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    echo -e "${BOLD}Speedup Analysis:${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    if [ "$speedup" != "N/A" ]; then
        echo -e "  ${GREEN}${BOLD}Rust is ${speedup}x faster than Python${NC}"

        local percentage=$(python3 -c "print(f'{($speedup - 1) * 100:.1f}')" 2>/dev/null || echo "N/A")
        if [ "$percentage" != "N/A" ]; then
            echo -e "  (${percentage}% performance improvement)"
        fi
    else
        echo "  Unable to calculate speedup"
    fi
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    if [ "$python_pos_err" != "N/A" ]; then
        echo -e "${BOLD}Accuracy Comparison:${NC}"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        printf "%-25s %15s %15s\n" "Error Metric" "Python" "Rust"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        printf "%-25s %15.2e %15.2e\n" "Position Error" "$python_pos_err" "$rust_pos_err"
        printf "%-25s %15.2e %15.2e\n" "Velocity Error" "$python_vel_err" "$rust_vel_err"
        printf "%-25s %15.2e %15.2e\n" "Energy Error" "$python_energy_err" "$rust_energy_err"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "Note: Both implementations use Euler integration, so accuracy should be similar."
        echo ""
    fi
}

# Generate summary report
generate_report() {
    print_section "Summary Report"

    local report_file="$SCRIPT_DIR/benchmark_report.txt"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    {
        echo "RustSim vs PathSim (Python) Benchmark Report"
        echo "Generated: $timestamp"
        echo "==========================================="
        echo ""
        echo "Test: Harmonic Oscillator Simulation"
        echo "  - System: d²x/dt² = -x"
        echo "  - Initial conditions: x(0) = 1.0, v(0) = 0.0"
        echo "  - Duration: 2π seconds (one period)"
        echo "  - Time step: 0.001 seconds"
        echo "  - Total steps: ~6283"
        echo ""

        if [ -f /tmp/python_output.txt ]; then
            echo "Python Results:"
            echo "---------------"
            cat /tmp/python_output.txt
            echo ""
        fi

        if [ -f /tmp/rust_output.txt ]; then
            echo "Rust Results:"
            echo "-------------"
            cat /tmp/rust_output.txt
            echo ""
        fi
    } > "$report_file"

    success "Report saved to: $report_file"
    echo ""
}

# Cleanup temporary files
cleanup() {
    if [ "$1" != "--keep-temp" ]; then
        rm -f /tmp/python_output.txt /tmp/rust_output.txt
        rm -f /tmp/python_benchmark_results.json /tmp/rust_benchmark_results.json
    fi
}

# Main execution
main() {
    print_banner

    # Parse arguments
    local skip_compile=false
    local python_only=false
    local rust_only=false
    local keep_temp=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            --skip-compile)
                skip_compile=true
                shift
                ;;
            --python-only)
                python_only=true
                shift
                ;;
            --rust-only)
                rust_only=true
                shift
                ;;
            --keep-temp)
                keep_temp=true
                shift
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --skip-compile    Skip Rust compilation step"
                echo "  --python-only     Run only Python benchmark"
                echo "  --rust-only       Run only Rust benchmark"
                echo "  --keep-temp       Keep temporary files"
                echo "  --help, -h        Show this help message"
                echo ""
                exit 0
                ;;
            *)
                error_exit "Unknown option: $1. Use --help for usage."
                ;;
        esac
    done

    # Check dependencies
    check_dependencies

    # Compile Rust benchmarks
    if [ "$skip_compile" = false ] && [ "$python_only" = false ]; then
        compile_rust_benchmarks
    fi

    # Run benchmarks
    if [ "$rust_only" = false ]; then
        run_python_benchmark
    fi

    if [ "$python_only" = false ]; then
        run_rust_benchmark
    fi

    # Compare results
    if [ "$python_only" = false ] && [ "$rust_only" = false ]; then
        compare_results
    fi

    # Generate report
    generate_report

    # Cleanup
    if [ "$keep_temp" = true ]; then
        cleanup --keep-temp
    else
        cleanup
    fi

    echo -e "${BOLD}${GREEN}Benchmark suite completed successfully!${NC}"
    echo ""
}

# Run main function
main "$@"
