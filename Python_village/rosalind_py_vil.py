"""Rosalind Python Village - Problem solver runner"""

from pythonvillage import PythonVillage
from pathlib import Path

# Get the path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / "Datasets"

py_vil = PythonVillage()

print("Installing Python")
py_vil.installing_python()

print("\nVariables and Some Arithmetics")
print(py_vil.variables_and_some_arithmetics(a=992,b=894))

print("\nStrings and Lists")
input_path = DATASETS_DIR / "rosalind_ini3.txt"
with open(input_path, "r") as f:
    text = f.read().strip()
a, b = int(text.split()[1]), int(text.split()[2])
c, d = int(text.split()[3]), int(text.split()[4])
print(py_vil.strings_and_lists(text, a, b, c, d))

print("\nConditions and Loops")
print(py_vil.conditions_and_loops(a=4774,b=9459))

print("\nWorking with Files")
input_path = DATASETS_DIR / "rosalind_ini5.txt"
output_path = DATASETS_DIR / "rosalind_output5.txt"
print(py_vil.working_with_files(input_path, output_path))

print("\nDictionaries")
input_path = DATASETS_DIR / "rosalind_ini6.txt"
with open(input_path, "r") as f:
    text = f.read()
word_dict = py_vil.dictionaries(text)
for word, count in word_dict.items():
    print(word, count)
