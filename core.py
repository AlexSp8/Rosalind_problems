"""
Core utilities for file I/O and data processing.
"""

from pathlib import Path
from typing import Dict, List


def file_to_list(file_path: str | Path) -> List[str]:
    """
    Read a file and return its lines as a list.

    Args:
        file_path: Path to the file to read

    Returns:
        List of stripped lines from the file
    """
    with open(file_path, 'r') as f:
        file_list = [l.strip() for l in f.readlines()]
    return file_list


def file_to_int_list(file_path: str | Path) -> List[List[int]]:
    """Returns a list of lists of integers from a file"""
    with open(file_path, 'r') as f:
        int_list = [ [int(x) for x in line.split()] for line in f]
    return int_list


def write_file(file_path: str | Path, content: str, mode: str = 'w') -> None:
    """
    Write content to a file.

    Args:
        file_path: Path to the file to write
        content: Content to write
        mode: File mode ('w' for write, 'a' for append)
    """
    with open(file_path, mode) as f:
        f.write(content + '\n')


def fasta_list_to_dict(fasta_list: List[str]) -> Dict[str, str]:
    """
    Convert a FASTA list to a dictionary.

    Args:
        fasta_list: List of lines from a FASTA file

    Returns:
        Dictionary with FASTA labels as keys and sequences as values
    """
    fasta_dict = {}
    fasta_label = ""
    for line in fasta_list:
        if '>' in line:
            fasta_label = line
            fasta_dict[fasta_label] = ""
        else:
            fasta_dict[fasta_label] += line
    return fasta_dict


def fasta_path_to_dict(fasta_path: str | Path) -> Dict[str, str]:
    """
    Read a FASTA file from a path and return as a dictionary.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        Dictionary with FASTA labels as keys and sequences as values
    """
    fasta_list = file_to_list(fasta_path)
    return fasta_list_to_dict(fasta_list)

def str_is_sorted(s: str) -> bool:
    """Check if a string is sorted"""
    return s == ''.join(sorted(s))
