# WarpX GPU Cylc Workflow on NERSC Perlmutter

Minimal Cylc workflow for submitting WarpX simulations to NERSC Perlmutter GPU nodes.

## Quick Start (5 Minutes)

### Prerequisites
```bash
pip install cylc-flow cylc-rose  # One-time setup
```

### Configuration
Edit `flow.cylc` and customize these values (lines 9-14):

```jinja2
{%- set proj = "YOUR_PROJECT_g" %}         # Your NERSC GPU project (must end with _g)
{%- set warpx_exe = "./warpx" %}           # Path to WarpX executable
{%- set inputs_file = "inputs" %}          # Input filename
{%- set num_nodes = 1 %}                   # Number of GPU nodes (1-4)
{%- set gpus_per_node = 4 %}               # GPUs per node (always 4)
{%- set wall_time_min = 10 %}              # Wall time in minutes
```

### Setup & Run
```bash
# 1. Copy profile to workflow directory (for environment setup)
cp perlmutter_gpu_warpx.profile.example my-workflow/

# 2. Install workflow
cylc install ~/src/warpx/Tools/machines/perlmutter-nersc/

# 3. Start workflow
cylc play workflow-name

# 4. Monitor in real-time
cylc tui workflow-name
```

## Workflow Structure

**3 simple tasks:**

```
setup
  ↓ (creates output directory)
warpx_gpu_run
  ↓ (runs on 4 GPUs + 32 CPUs per node, GPU-aware MPI enabled)
check
  ↓ (verifies output files exist)
```

## Key Features

✅ **GPU Optimized**
- 4 GPUs per node (A100, Perlmutter standard)
- 32 CPU cores per node
- GPU-aware MPI enabled
- CUDA A100 optimization (AMREX_CUDA_ARCH=8.0)
- Automatic CUDA device mapping per rank

✅ **Reuses Existing Scripts**
- Environment setup sourced from `perlmutter_gpu_warpx.profile.example`
- Execution logic reuses proven commands from `perlmutter_gpu.sbatch`
- No duplication of configuration

✅ **Simple & Maintainable**
- 119 lines of clean, readable workflow
- Configuration in one place (lines 9-14)
- Automatic profile detection and sourcing

## Common Adjustments

### Run on 2 GPU nodes (8 GPUs total)
```jinja2
{%- set num_nodes = 2 %}
```

### Extend wall time to 30 minutes
```jinja2
{%- set wall_time_min = 30 %}
```

### Change input file
```jinja2
{%- set inputs_file = "my_simulation_3d" %}
```

### Use different profile location
```jinja2
{%- set profile_script = "/path/to/custom/profile.example" %}
```

## Commands Reference

```bash
# Install workflow (register with Cylc)
cylc install ~/src/warpx/Tools/machines/perlmutter-nersc/

# Start workflow
cylc play workflow-name

# Monitor with terminal UI (recommended)
cylc tui workflow-name

# View scheduler log
cylc cat-log workflow-name scheduler

# View task logs
cylc cat-log workflow-name warpx_gpu_run

# Stop workflow
cylc stop workflow-name

# Remove workflow
cylc clean workflow-name

# Validate syntax before running
cylc validate flow.cylc

# List all workflows
cylc list
```

## Environment Setup

The workflow automatically sources `perlmutter_gpu_warpx.profile.example` which provides:

- Module loads: gpu, PrgEnv-gnu, cudatoolkit, cmake, hdf5, python
- Build paths: CMAKE_PREFIX_PATH, LD_LIBRARY_PATH
- CUDA settings: AMREX_CUDA_ARCH=8.0, CRAY_ACCEL_TARGET=nvidia80
- Compiler flags: -march=znver3 optimization for Milan CPUs

**The profile is sourced automatically** - no manual setup needed.

## GPU-Specific Settings

The workflow enables:
- GPU-aware MPI: `export MPICH_OFI_NIC_POLICY=GPU`
- Per-rank CUDA device mapping: `export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))`
- OpenMP threading: `export OMP_NUM_THREADS=16`

These are sourced from proven NERSC configurations.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "cylc: command not found" | Install: `pip install cylc-flow` |
| "Project account not set" | Edit flow.cylc line 9, set `proj = "YOUR_PROJECT_g"` (must end with `_g`) |
| "Executable not found" | Verify WarpX path in flow.cylc line 10 |
| "Input file not found" | Verify input filename in flow.cylc line 11 matches actual file |
| "Syntax error in flow.cylc" | Run: `cylc validate flow.cylc` |
| "Profile not found" | Copy `perlmutter_gpu_warpx.profile.example` to workflow dir |
| "GPU job won't submit" | Check: `getbal` (GPU allocation), `sinfo -p gpu_ss` (node availability) |
| "Workflow stalls" | Check logs: `cylc cat-log workflow-name scheduler` |

## File Organization

```
perlmutter-nersc/
├── flow.cylc                                 ⭐ Main GPU workflow
├── README.md                                 📖 This file
├── perlmutter_gpu.sbatch                     (Reference SLURM script)
├── perlmutter_gpu_warpx.profile.example      (Environment setup - sourced by cylc)
├── install_gpu_dependencies.sh               (GPU dependencies installer)
├── perlmutter_cpu.sbatch                     (CPU variant reference)
├── perlmutter_cpu_warpx.profile.example      (CPU variant reference)
└── container/                                (Docker setup)
```

## Setup Examples

### Example 1: Simple Single-Node Run
```bash
# Create workflow directory
mkdir my-warpx-sim
cd my-warpx-sim

# Copy files
cp ~/src/warpx/Tools/machines/perlmutter-nersc/flow.cylc .
cp ~/src/warpx/Tools/machines/perlmutter-nersc/perlmutter_gpu_warpx.profile.example .

# Edit flow.cylc
sed -i 's/YOUR_PROJECT_g/m3912_g/' flow.cylc           # Your project
sed -i 's|./warpx|/path/to/warpx/build_gpu/bin/warpx|' flow.cylc  # Your exe
sed -i 's/inputs/my_inputs_3d/' flow.cylc              # Your input

# Install and run
cylc install .
cylc play my-warpx-sim
cylc tui my-warpx-sim
```

### Example 2: Multi-Node Parameter Study
```jinja2
# In flow.cylc, set:
{%- set num_nodes = 4 %}
{%- set wall_time_min = 60 %}
# Will run on 4 GPU nodes (16 GPUs total) with 60 minute time limit
```

### Example 3: Custom Profile
```jinja2
# If you have a custom profile at /home/edoardo/my_warpx.profile:
{%- set profile_script = "/home/edoardo/my_warpx.profile" %}
```

## Performance Notes

- **Perlmutter GPU nodes**: 4 × A100 GPUs per node, 32 CPU cores
- **GPU-aware MPI**: Significantly faster GPU communication
- **CUDA-aware device mapping**: Each MPI rank gets correct GPU
- **Node constraint**: Automatically selects GPU-capable nodes (`--constraint gpu`)
- **Job queue**: Typical wait time 30 min - 1 hour on regular QoS

## Pre-Flight Checklist

```
☐ Cylc installed: cylc --version
☐ WarpX executable compiled and at specified path
☐ Input file ready and filename matches flow.cylc
☐ NERSC GPU project set (must end with _g)
☐ GPU allocation available: getbal
☐ Profile file exists: perlmutter_gpu_warpx.profile.example
☐ Syntax valid: cylc validate flow.cylc
☐ SSH access to Perlmutter working
```

## Advanced Usage

### Output Management
Outputs go to `output/warpx.out` by default. Modify in `warpx_gpu_run` task:
```bash
script = """
    # Change output location:
    srun ... > my_custom_output.log 2>&1
"""
```

### Custom Modules
Add module loads by editing `perlmutter_gpu_warpx.profile.example` before sourcing.

### Monitoring Multiple Runs
```bash
cylc list                    # List all active workflows
cylc tui                     # Interactive UI shows all workflows
cylc show workflow-name      # Show specific workflow state
```

### Re-running a Failed Task
```bash
cylc trigger workflow-name warpx_gpu_run
```

### Checking Job Queue
```bash
squeue -u $USER                 # Your jobs on Perlmutter
sinfo -p gpu_ss                 # GPU partition info
```

## Resources & Documentation

- **Cylc**: https://cylc.github.io/cylc-doc/stable/html/
- **WarpX**: https://warpx.readthedocs.io/
- **NERSC Perlmutter**: https://docs.nersc.gov/systems/perlmutter/
- **NERSC GPU Info**: https://docs.nersc.gov/systems/perlmutter/architecture/#gpu-nodes

## Summary

✅ **Minimal Cylc workflow** for WarpX GPU simulations on Perlmutter  
✅ **Easy configuration** - edit 6 lines in flow.cylc  
✅ **Reuses proven scripts** - no duplication  
✅ **GPU-optimized** - GPU-aware MPI, proper device mapping  
✅ **Well-documented** - this guide covers all aspects  

**Ready to use!** Edit flow.cylc → Install → Run → Monitor with cylc tui

---

Last Updated: January 2026 | For Cylc v8.x on Perlmutter GPUs
