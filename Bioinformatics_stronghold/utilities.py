
def read_file(file_path):
    with open(file_path, 'r') as f:
        list = [l.strip() for l in f.readlines()]
    return list

def writeFile(filePath, seq, mode ='w'):
    #'w' write, 'a' append
    with open(filePath, mode) as f:
        f.write(seq + '\n')

def fasta_to_dict(fasta_path):
    fasta_file = read_file(fasta_path)
    # print(fasta_file)
    fasta_dict = {}
    fasta_label = ""
    for line in fasta_file:
        if '>' in line:
            fasta_label = line
            fasta_dict[fasta_label] = ""
        else:
            fasta_dict[fasta_label] += line
    return fasta_dict
