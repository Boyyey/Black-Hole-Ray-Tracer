# 🌌 Relativistic Black Hole Ray Tracer: The Ultimate Research Tool

## 🔭 Overview

**A production-grade, research-focused black hole simulation framework** that implements **full general relativistic ray tracing** through the Schwarzschild metric. This isn't just a visualization tool—it's a **computational physics laboratory** designed for actual astrophysics research that scientists can use in peer-reviewed publications. 📚

> "This simulation implements the same physics used in professional astrophysics codes and even the visual effects for the movie *Interstellar*" - *Your future citation*

## 🌟 Key Features

- **🔬 Scientific Accuracy**: Solves full Schwarzschild geodesic equations with adaptive step integration
- **📊 Research-Grade Output**: Generates CSV data suitable for quantitative analysis
- **⚙️ Parameter Studies**: Systematically vary camera position, viewing angle, and disk properties
- **✅ Validation Suite**: Verify against known analytical solutions (photon sphere, ISCO, redshift)
- **🚀 High Performance**: OpenMP parallelization for multi-core systems
- **📈 Parameter Sweeps**: Automate systematic studies with custom ranges
- **⏱️ Benchmark Mode**: Measure performance for different configurations
- **📦 Reproducible Science**: Complete telemetry metadata for every simulation

## ⚙️ Installation & Compilation

### Prerequisites
- GCC or compatible C compiler 🐧
- Make (optional but recommended) 🛠️
- OpenMP (for parallel execution) ⚡

### Compilation

**Linux/MacOS:**
```bash
gcc -O3 -fopenmp blackhole.c -lm -o blackhole
```

**Windows (MinGW64):**
```bash
gcc -O3 blackhole.c -lm -o blackhole.exe
```

**With Debug Symbols:**
```bash
gcc -O3 -g -fopenmp blackhole.c -lm -o blackhole_debug
```

## 🚀 Usage

### Basic Simulation
```bash
./blackhole -w 1024 -h 1024 -d 15 -t 90 -i 6 -o 30
```

### Parameter Sweep (Vary Viewing Angle)
```bash
./blackhole -s theta -r 0:90:5 --prefix my_study
```

### Validation Mode (Verify Against Known Solutions)
```bash
./blackhole --validate
```

### Benchmark Performance
```bash
./blackhole -b
```

### Full Help
```bash
./blackhole --help
```

## 📋 Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-w`, `--width` | Image width | 1024 |
| `-h`, `--height` | Image height | 1024 |
| `-d`, `--distance` | Camera distance (Rs) | 15.0 |
| `-t`, `--theta` | Camera polar angle (degrees) | 90.0 |
| `-i`, `--disk-inner` | Inner disk radius (Rs) | 6.0 |
| `-o`, `--disk-outer` | Outer disk radius (Rs) | 30.0 |
| `-m`, `--max-steps` | Max integration steps | 2048 |
| `-s`, `--sweep` | Parameter sweep mode | - |
| `-r`, `--range` | Parameter range (START:END:STEP) | - |
| `-b`, `--benchmark` | Performance benchmark | - |
| `-v`, `--validate` | Run validation suite | - |
| `-p`, `--prefix` | Output filename prefix | blackhole |

## 📊 Output Files

- `blackhole_photon_data.csv` - Detailed per-photon data for analysis
- `blackhole_telemetry.csv` - Simulation metadata and statistics
- `blackhole_sweep_*.csv` - Parameter sweep results

## 🔬 Research Applications

This tool is designed for:

- **Critical Parameter Studies** 📏
  - Analyze photon capture boundaries
  - Study redshift distribution across disk

- **Monte Carlo Simulations** 🎲
  - Statistical analysis of photon paths
  - Probability density functions

- **Gravitational Lensing Research** 🔍
  - Multiple image analysis
  - Magnification calculations

- **Disk Structure Modeling** 🌀
  - Temperature gradient studies
  - Emission line analysis

- **Event Horizon Telescope Comparisons** 📡
  - Shadow shape analysis
  - Photon ring measurements

- **Validation of Theoretical Predictions** ✅
  - Verify against analytical solutions
  - Quantify numerical errors

## 🧪 Validation Suite

The simulation includes a built-in validation suite that verifies against known analytical solutions:

```bash
$ ./blackhole --validate
Validating photon sphere at r = 3.0000 Rs
  ✓ Photon capture confirmed for b = 5.1954 (critical b = 5.1962)
Validating ISCO at r = 6.0000 Rs
  ✓ ISCO validation confirmed: disk hit at r = 6.0000 Rs
Validating gravitational redshift formula
  ✓ Gravitational redshift validation: z = 0.1547 (expected 0.1547), error = 1.11e-16

Validation suite complete: 0 failures out of 3 tests
```

## 🤝 Contributing

Contributions are welcome! Here's how you can help:

1. 🐛 **Report bugs** through GitHub Issues
2. 💡 **Suggest features** for research applications
3. 🧪 **Add validation tests** against additional analytical solutions
4. 📈 **Improve performance** with better algorithms
5. 📚 **Extend to Kerr metric** for rotating black holes

Please read our [Contributing Guidelines](CONTRIBUTING.md) before submitting a PR.

---

*Simulating the most extreme objects in the universe, one ray at a time.* 🌌  
*Designed for researchers, by a researcher.* 🔬  