# Plan: Generate Cursor Rules for YClade

## Overview

Generate comprehensive `.cursor/rules/*.mdc` files in MDC format for the YClade repository, emphasizing TDD-first development and Python best practices.

## Analysis

### Current State
- Python 3.10+ bioinformatics project
- Uses: Black (100 chars), Ruff, MyPy (strict), pytest
- TDD workflow established per HANDOFF.md
- src-layout with modular design
- Dataclasses for data structures
- Context managers for resources
- Google-style docstrings

### Requirements
- MDC format with YAML frontmatter
- Rules should cover: general conventions, Python standards, TDD, code style, error handling, git, bioinformatics domain
- Emphasize "fail hard" philosophy (no fallbacks)
- Surgical changes principle

## Tasks

- [x] Create `.cursor/rules/` directory
- [x] Create `general.mdc` - Core project conventions
- [x] Create `python.mdc` - Python best practices (types, docstrings, imports)
- [x] Create `tdd.mdc` - TDD workflow and testing standards
- [x] Create `code-style.mdc` - Black, Ruff, MyPy formatting rules
- [x] Create `error-handling.mdc` - Fail hard philosophy
- [x] Create `git.mdc` - Commit conventions
- [x] Create `bioinformatics.mdc` - Domain-specific conventions

## Files Created

| File | Description | Key Content |
|------|-------------|-------------|
| `general.mdc` | Core conventions | Architecture, structure, dependencies, patterns |
| `python.mdc` | Python standards | Types, docstrings, imports, naming, paths |
| `tdd.mdc` | Testing standards | TDD cycle, fixtures, assertions, markers |
| `code-style.mdc` | Formatting | Black, Ruff, MyPy configs, line length |
| `error-handling.mdc` | Error philosophy | Fail hard, no fallbacks, exit codes |
| `git.mdc` | Git workflow | Conventional commits, branch naming |
| `bioinformatics.mdc` | Domain rules | Haplogroups, coordinates, VCF, ancient DNA |

## Status

✅ **Complete** - All 7 MDC rule files generated and written to `.cursor/rules/`

