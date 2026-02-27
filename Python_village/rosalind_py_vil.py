from pythonvillage import PythonVillage

py_vil = PythonVillage()

print("Installing Python")
py_vil.installing_python()

print("\nVariables and Some Arithmetics")
print(py_vil.variables_and_some_arithmetics(a=992,b=894))

print("\nStrings and Lists")
text = "iQ1xFPLPKS3zVBtqkPTjXkDKNPH9GnlEevaJw53tTSZyUcyRNfohk9" \
"C6aalwb6uJo3fLycaenopsisRDi4GjieTCnyON9Z85ZdyIvgZ0CZW4lILDrByT1" \
"LmB1f28CoCucoelestinuslBBRuQqu5pREhvvWlwta"
print(py_vil.strings_and_lists(text, a=67, b=77, c=128, d=138))

print("\nConditions and Loops")
print(py_vil.conditions_and_loops(a=4774,b=9459))

print("\nWorking with Files")
input_path = 'Python_village/Datasets/rosalind_ini5.txt'
output_path = 'Python_village/Datasets/rosalind_output5.txt'
print(py_vil.working_with_files(input_path, output_path))

print("\nDictionaries")
text = "When I find myself in times of trouble Mother Mary comes to me " \
"Speaking words of wisdom let it be And in my hour of darkness she is " \
"standing right in front of me Speaking words of wisdom let it be Let it " \
"be let it be let it be let it be Whisper words of wisdom let it be And when " \
"the broken hearted people living in the world agree There will be an answer " \
"let it be For though they may be parted there is still a chance that they " \
"will see There will be an answer let it be Let it be let it be let it be let " \
"it be There will be an answer let it be Let it be let it be let it be let it " \
"be Whisper words of wisdom let it be Let it be let it be let it be let it be " \
"Whisper words of wisdom let it be And when the night is cloudy there is still " \
"a light that shines on me Shine until tomorrow let it be I wake up to the sound " \
"of music Mother Mary comes to me Speaking words of wisdom let it be Let it be " \
"let it be let it be yeah let it be There will be an answer let it be Let it be " \
"let it be let it be yeah let it be Whisper words of wisdom let it be"
word_dict = py_vil.dictionaries(text)
for word, count in word_dict.items():
    print(word, count)
