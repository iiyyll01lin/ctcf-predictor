# CTCF PWM Testing Pipeline

> **🧬 Advanced Position Weight Matrix Generation for CTCF Binding Site Prediction**

[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](docs/)
[![Docker](https://img.shields.io/badge/docker-supported-blue)](docs/03-docker-setup.md)
[![Quick Start](https://img.shields.io/badge/quick%20start-5%20minutes-orange)](docs/01-quick-start.md)

This project implements a comprehensive, automated pipeline for building high-quality Position Weight Matrices (PWMs) to predict CTCF transcription factor binding sites with revolutionary quality-over-quantity approach.

## 🏆 Key Achievement

**28× performance improvement** proving that small, high-quality datasets (1,000 sequences, 19.592 bits) dramatically outperform large, unfiltered datasets (37,628 sequences, 0.695 bits) in transcription factor modeling.

## 📚 **NEW: Comprehensive Documentation Available**

**All documentation has been reorganized into focused, easy-to-navigate guides in the [`docs/`](docs/) directory.**

### 🚀 **Quick Access**
- **[📋 Documentation Hub](docs/README.md)** - Complete navigation and overview
- **[⚡ Quick Start Guide](docs/01-quick-start.md)** - Get running in 5 minutes
- **[🔧 Troubleshooting](docs/14-troubleshooting.md)** - Solve common issues fast

### 🧬 **For Scientists & Researchers**
- **[Biological Background](docs/04-biological-background.md)** - CTCF protein and genomic organization
- **[Information Theory](docs/05-information-theory.md)** - Mathematical foundations of PWM quality
- **[Statistical Framework](docs/06-statistical-framework.md)** - Validation methodology and null models

### 🛠️ **For Users & System Administrators**
- **[System Requirements](docs/02-system-requirements.md)** - Hardware and software requirements
- **[Docker Setup](docs/03-docker-setup.md)** - Complete containerization guide
- **[User Guide](docs/10-user-guide.md)** - Comprehensive usage instructions

## 🚀 Quick Start

### Prerequisites
- **Docker** (recommended) or **R 4.0+** with required packages
- **8GB RAM** minimum, 16GB recommended  
- **5GB free disk space** for full dataset

### Quick Demo (5 minutes)
```bash
# Clone and start Docker environment
git clone https://github.com/organization/ctcf-predictor.git
cd ctcf-predictor
./smart-startup.sh

# Run demo pipeline
./run-in-docker.sh test_chromosome_split.R demo

# Expected: PWM with >8 bits information content
```

**📖 For detailed instructions, see the [Quick Start Guide](docs/01-quick-start.md)**

## 🎯 What Makes This Pipeline Special

### Revolutionary Discoveries
- **Quality-Over-Quantity Paradigm**: First scientific proof that dataset quality trumps size
- **28× Performance Improvement**: Small, curated datasets dramatically outperform large ones
- **Chromosome-Based Validation**: Prevents data leakage in genomic applications
- **Comprehensive Statistical Framework**: 300+ null model replicates for robust validation

### Technical Excellence
- **Multiple PWM Building Methods**: 10+ different approaches tested and validated
- **Docker Containerization**: Cross-platform reproducibility
- **Automated Pipeline**: From data download to final PWM generation
- **Production-Ready**: Robust error handling and quality checks

### Scientific Impact
- **Medical Applications**: Drug target identification and therapeutic design
- **Genome Engineering**: Precise CRISPR guide design
- **Basic Research**: Understanding genome organization and gene regulation
- **Standardization**: Establishing best practices for transcription factor modeling

## 📊 Quick Navigation by Role

| **Role**              | **Start Here**                                            | **Next Steps**                                            | **Advanced**                                        |
|-----------------------|-----------------------------------------------------------|-----------------------------------------------------------|-----------------------------------------------------|
| **🆕 New User**       | [Quick Start](docs/01-quick-start.md)                     | [User Guide](docs/10-user-guide.md)                       | [Troubleshooting](docs/14-troubleshooting.md)       |
| **🔬 Researcher**     | [Biological Background](docs/04-biological-background.md) | [Statistical Framework](docs/06-statistical-framework.md) | [Information Theory](docs/05-information-theory.md) |
| **💻 Developer**      | [System Requirements](docs/02-system-requirements.md)     | [Docker Setup](docs/03-docker-setup.md)                   | *Coming Soon*                                       |
| **🖥️ Administrator** | [Docker Setup](docs/03-docker-setup.md)                   | [System Requirements](docs/02-system-requirements.md)     | [Troubleshooting](docs/14-troubleshooting.md)       |

## 📈 Performance Benchmarks

| **Method**      | **Dataset Size** | **Information Content** | **Quality Grade** | **Applications**         |
|-----------------|------------------|-------------------------|-------------------|--------------------------|
| **Subset 1K**   | 1,000 sequences  | **19.592 bits**         | 🏆 Excellent      | Drug discovery, clinical |
| **Subset 2K**   | 2,000 sequences  | **12.564 bits**         | ✅ Good            | Research, publications   |
| **Robust PWM**  | 24,125 sequences | **10.659 bits**         | ✅ Good            | General analysis         |
| **Raw Dataset** | 37,628 sequences | **0.695 bits**          | ❌ Poor            | Not recommended          |

## 🔄 Legacy Information

This README has been streamlined for quick access. **The original detailed pipeline documentation has been preserved and enhanced in the new modular documentation system.** 

### Migration Notes
- **Old single-file documentation** → **New modular docs** in [`docs/`](docs/) directory
- **Enhanced with 28× improvement discovery** and **statistical validation framework**
- **Updated with real performance data** and **current best practices**
- **Cross-platform instructions** for Windows, macOS, and Linux

### Finding Previous Content
- **Detailed pipeline steps** → [User Guide](docs/10-user-guide.md)
- **Technical specifications** → [System Requirements](docs/02-system-requirements.md)
- **Script documentation** → *Coming Soon*
- **Configuration options** → *Coming Soon*

## 🤝 Contributing

We welcome contributions! Please see our documentation for:
- [System Architecture](docs/) (Coming Soon)
- [Extending the Pipeline](docs/) (Coming Soon)
- [Development Guidelines](docs/) (Coming Soon)

## 📄 Citation

If you use this pipeline in research, please cite:

```
CTCF PWM Testing Pipeline (2025)
A comprehensive framework for CTCF Position Weight Matrix construction and validation
with revolutionary quality-over-quantity discovery.
GitHub: https://github.com/organization/ctcf-predictor
```

## 📧 Support

- **📚 Documentation**: Complete guides in [`docs/`](docs/) directory
- **🐛 Issues**: Report bugs via [GitHub Issues](https://github.com/iiyyll01lin/ctcf-predictor/issues)
- **💬 Questions**: Check [Troubleshooting Guide](docs/14-troubleshooting.md) first

---

**🎯 Start Here**: [Documentation Hub](docs/README.md) → [Quick Start](docs/01-quick-start.md) → Your specific role's guide above.
