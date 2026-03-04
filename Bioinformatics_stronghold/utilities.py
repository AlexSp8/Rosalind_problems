
def file_to_list(file_path):
    with open(file_path, 'r') as f:
        file_list = [l.strip() for l in f.readlines()]
    return file_list

def write_file(filePath, seq, mode='w'):
    #'w' write, 'a' append
    with open(filePath, mode) as f:
        f.write(seq + '\n')

def fasta_list_to_dict(fasta_file):
    fasta_dict = {}
    fasta_label = ""
    for line in fasta_file:
        if '>' in line:
            fasta_label = line
            fasta_dict[fasta_label] = ""
        else:
            fasta_dict[fasta_label] += line
    return fasta_dict

def fasta_path_to_dict(fasta_path):
    fasta_list = file_to_list(fasta_path)
    return fasta_list_to_dict(fasta_list)
