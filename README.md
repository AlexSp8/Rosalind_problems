# Rosalind Problems

A comprehensive Python solution library for [Rosalind](http://rosalind.info/) bioinformatics problems.

## Project Structure

```
Rosalind_problems/
├── Bioinformatics_stronghold/      # Advanced bioinformatics problems
│   ├── biosequence.py              # Core DNA/RNA sequence class
│   ├── biostronghold.py            # Problem solvers for stronghold track
│   ├── Datasets/                   # Input and output files for problems
│   └── Notes/                      # Problem-specific notes
├── Python_village/                 # Introductory Python problems
│   ├── pythonvillage.py            # Simple Python problem solvers
│   └── Datasets/                   # Input and output files
├── config.py                       # Centralized configuration and constants
├── core.py                         # Shared I/O utilities and helpers
└── __init__.py                     # Package initialization
```

## Features

### Bioinformatics Stronghold

Advanced bioinformatics problem solutions including:

- **DNA/RNA Operations**: Transcription, complement, reverse complement
- **Sequence Analysis**: GC content, nucleotide counting, point mutations
- **Translation**: RNA to protein, codon tables
- **Pattern Matching**: Motif finding, shared motifs
- **Sequence Analysis**: Consensus strings, profile matrices
- **Genetics**: Fibonacci sequences, Mendelian probability, allele inheritance
- **Protein Matching**: Open reading frames, protein motif matching, mRNA inference

### Python Village

Introductory Python problem solutions:

- Variables and arithmetic
- String manipulation and list operations
- Control flow (conditions and loops)
- File I/O
- Dictionary operations

## Usage

### Running Bioinformatics Problems

```python
from Bioinformatics_stronghold import BioStronghold

# Create solver instance
solver = BioStronghold(seq="ACGTACGTACGT")

# DNA transcription
rna = solver.transcribing_DNA_into_RNA()
print(rna)  # ACGUACGUACGU

# Count nucleotides
counts = solver.counting_DNA_nucleotides()
print(counts)

# Reverse complement
rev_comp = solver.complementing_a_strand_of_DNA()
print(rev_comp)
```

### Running Python Village Problems

```python
from Python_village import PythonVillage

solver = PythonVillage()

# Calculate sum of squares
result = solver.variables_and_some_arithmetics(3, 4)
print(result)  # 25

# Count words
word_count = solver.dictionaries("hello world hello")
print(word_count)  # {'hello': 2, 'world': 1}
```

## Core Classes

### Biosequence

Base class for DNA/RNA sequence manipulation:

- `transcribe_sequence()` - Convert DNA ↔ RNA
- `complement_sequence()` - Get nucleotide complement
- `reverse_complement_sequence()` - Get reverse complement
- `get_gc_content()` - Calculate GC percentage
- `get_point_mutations()` - Count differences
- `translate_sequence()` - Convert to amino acids
- `motif_locations_in_sequence()` - Find pattern matches
- And many more...

### BioStronghold

Extends Biosequence with specific problem solvers for Rosalind challenges.

## Configuration

Constants are centralized in `config.py`:

```python
from config import (
    BASE_NUCLEOTIDES,           # Valid nucleotides for DNA/RNA
    RNA_CODON_TO_AMINO,         # RNA codon translation table
    DNA_CODON_TO_AMINO,         # DNA codon translation table
)
```

## Type Hints

The codebase uses Python type hints for better IDE support and type checking:

```python
def translate_sequence(self, istart: int = 0, seq: Optional[str] = None) -> str:
    """Translate DNA/RNA sequence to protein"""
```

## Requirements

- Python 3.7+
- Standard library only (no external dependencies)

## File I/O Utilities

The `core.py` module provides utilities for working with FASTA files:

```python
from core import fasta_path_to_dict, file_to_list, write_file

# Read FASTA file
sequences = fasta_path_to_dict("path/to/sequences.fasta")

# Read text file
lines = file_to_list("path/to/file.txt")

# Write output
write_file("output.txt", result)
```

## Notes

- All sequence inputs are automatically converted to uppercase
- Positions in output are 1-indexed (following Rosalind conventions)
- Multiple solutions use dynamic programming and efficient algorithms

## Future Improvements

- [ ] Integrate more Rosalind problems solutions
- [ ] Add visualization tools for sequences
- [ ] Performance optimization for large datasets
- [ ] Create jupyter notebook tutorials
