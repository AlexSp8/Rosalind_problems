
class PythonVillage():

    """Class that contains functions to solve all Python Village problems in Rosalind"""

    def __init__(self):
        pass

    @staticmethod
    def installing_python():
        import this

    @staticmethod
    def variables_and_some_arithmetics(a, b):
        return a**2 + b**2

    @staticmethod
    def strings_and_lists(text, a, b, c, d):
        slice1 = text[a:b+1]
        slice2 = text[c:d+1]
        return f"{slice1} {slice2}"

    @staticmethod
    def conditions_and_loops(a, b):
        return sum(i for i in range(a, b+1) if i % 2 == 1)

    @staticmethod
    def working_with_files(input_path, output_path):
        with open(input_path,"r") as f_in, open(output_path, "w") as f_out:
            f_out.writelines(line for i, line in enumerate(f_in) if i % 2 == 1)
        return output_path

    @staticmethod
    def dictionaries(text):
        # from collections import Counter
        # return Counter(text.split())
        word_dict = {}
        for word in text.split():
            word_dict[word] = word_dict.get(word, 0) + 1
        return word_dict
