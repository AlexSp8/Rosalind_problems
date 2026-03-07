
"""Python Village - Introductory Python problem solutions."""

from typing import Dict
from pathlib import Path


class PythonVillage:
    """Class that contains functions to solve all Python Village problems in Rosalind"""

    def __init__(self) -> None:
        """Initialize PythonVillage solver"""
        pass

    @staticmethod
    def installing_python() -> None:
        """Print the Zen of Python"""
        import this

    @staticmethod
    def variables_and_some_arithmetics(a: int, b: int) -> int:
        """Calculate sum of squares"""
        return a ** 2 + b ** 2

    @staticmethod
    def strings_and_lists(text: str, a: int, b: int, c: int, d: int) -> str:
        """Extract and return two slices from a string"""
        slice1 = text[a:b + 1]
        slice2 = text[c:d + 1]
        return f"{slice1} {slice2}"

    @staticmethod
    def conditions_and_loops(a: int, b: int) -> int:
        """Sum odd numbers in range [a, b]"""
        return sum(i for i in range(a, b+1) if i % 2 == 1)

    @staticmethod
    def working_with_files(input_path: str | Path, output_path: str | Path) -> str:
        """Copy odd-indexed lines from input to output file"""
        with open(input_path, "r") as f_in, open(output_path, "w") as f_out:
            f_out.writelines(line for i, line in enumerate(f_in) if i % 2 == 1)
        return str(output_path)

    @staticmethod
    def dictionaries(text: str) -> Dict[str, int]:
        """Count word frequencies in text"""
        word_dict = {}
        for word in text.split():
            word_dict[word] = word_dict.get(word, 0) + 1
        return word_dict
