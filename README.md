# CDCC Emulator: Human-AI Collaborative Development

This repository contains the source code and development logs for a reduced-basis emulator for Continuum-Discretized Coupled-Channel (CDCC) nuclear reaction calculations, developed through human-AI collaboration.

## Overview

This emulator was developed in **7 days** (December 16-19, 2025) using Claude Opus 4.5 via the Claude Code CLI, compared to an estimated 3-6 months for conventional development. The project demonstrates a reproducible workflow for AI-assisted computational physics research.

### Key Results
- **220× speedup** over direct CDCC calculations
- **<0.1% accuracy** compared to exact solutions
- **37 coupled channels** with 18 optical potential parameters
- Complete development from design to working code in one week

## Repository Contents

```
AI_emulator/
├── README.md                    # This file
├── Makefile                     # Build configuration
├── setup.sh                     # Environment setup script
│
├── Source Code/
│   ├── emulator_train.F90       # Training module (POD basis construction)
│   ├── emulator_predict.F90     # Prediction module (Galerkin projection)
│   ├── emulator_io.F90          # I/O utilities
│   ├── param_mapping.F90        # Parameter space mapping
│   └── test_*.F90               # Test and validation programs
│
├── Documentation/
│   ├── EMULATOR_README.md       # Technical documentation
│   ├── EMULATOR_INPUT_DESIGN.md # Input file specification
│   └── UNIFIED_BASIS_EMULATOR.md # Design notes
│
├── conversation_logs/           # Claude Code session logs
│   ├── 2025-12-16_session*.jsonl
│   ├── 2025-12-17_session*.jsonl
│   └── 2025-12-19_session*.jsonl
│
└── git_history/                 # Development history
    ├── commits_detailed.log     # Full commit messages
    ├── commits_oneline.log      # Concise commit list
    └── commits_per_day.log      # Commit statistics
```

## Development Timeline

| Date | Commits | Focus |
|------|---------|-------|
| Dec 16 | 2 | Initial design, core implementation |
| Dec 17 | 4 | Multi-J support, EIM integration |
| Dec 18 | 11 | Parameter handling, I/O, validation |
| Dec 19 | 7 | BLAS optimization, accuracy testing |

## Requirements

- Fortran compiler (gfortran recommended)
- LAPACK/BLAS libraries
- OpenMP support (optional, for parallelization)

## Building

```bash
# Set up environment (installs OpenBLAS if needed)
./setup.sh

# Build emulator programs
make
```

## Related Publications

1. **Technical paper**: J. Lei, "Reduced-basis emulator for continuum-discretized coupled-channel calculations," arXiv:2512.17687 (submitted to Physical Review C)

2. **Methodology paper**: J. Lei, "Human-AI collaboration compresses computational physics development cycle from months to days" (submitted to Nature Computational Science)

## Conversation Logs

The `conversation_logs/` directory contains complete Claude Code session logs in JSON Lines format. Each file represents a development session with:
- User prompts and instructions
- AI responses and code generation
- Tool calls (file reads, writes, bash commands)
- Iterative debugging and refinement

These logs provide full transparency into the human-AI collaborative development process.

## License

This code is provided for research and educational purposes. Please cite the relevant publications if you use this work.

## Repository

https://github.com/jinleiphys/AI_emulator

## Contact

Jin Lei

## Acknowledgments

This work was developed using Claude Opus 4.5 (Anthropic) via the Claude Code CLI interface.
