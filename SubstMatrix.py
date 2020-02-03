class SubstMatrix:
    '''
    implementa uma matriz de substituição
    '''

    def __init__(self):
        self.alphabet = ''
        self.sm = {}
        
    def loadFromFile(self, filename, sep):
        '''
        ler matrizes de substituiçao de um ficheiro
        '''
        f = open(filename, 'r')
        line = f.readline()
        tokens = line.split(sep)
        ns = len(tokens)
        self.alphabet = ''                 # _
        for i in range(ns):             #  |- 1a linha: alfabeto
            self.alphabet += tokens[i][0]  # _|
        for i in range(ns):                            # _
            line = f.readline()                           #  |
            tokens = line.split(sep)                      #  |
            for j in range(len(tokens)):               #  |- outras linhas: scores
                k = self.alphabet[i] + self.alphabet[j]   #  |
                self.sm[k] = int(tokens[j])               # _|
        f.close()
        return None
                
    def createFromMatchPars(self, match, mismatch, alphabet):
        self.alphabet = alphabet
        for c1 in alphabet:
            for c2 in alphabet:
                if c1 == c2:
                    self.sm[c1+c2] = match
                else:
                    self.sm[c1+c2] = mismatch
        return None
    
    def scorePair(self, c1, c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1+c2]
    
    def __getitem__(self, ij):
        '''
        permite acesso com listas []
        '''
        i, j = ij
        return self.scorePair(i, j)


# --------------------------\\-----------------------------

# TESTS

def test1():
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    print(sm.alphabet)
    print(sm.scorePair('G', 'M'))
    print(sm.scorePair('W', 'W'))
    print(sm.scorePair('A', 'S'))
    print(sm.scorePair('X', 'X'))
    print(sm['G', 'K'])
    print(sm['T', 'T'])

# Result test 1
# ARNDCQEGHILKMFPSTWYV
# -3
# 11
# 1
# None
# -2
# 5


def test2():
    sm = SubstMatrix()
    sm.createFromMatchPars(3, -1, 'ACGU')
    print(sm.alphabet)
    print(sm.scorePair('A', 'A'))
    print(sm.scorePair('A', 'U'))
    print(sm.scorePair('T', 'T'))
    print(sm['G', 'G'])

# Result test 2   
# ACGU
# 3
# -1
# None
# 3


if __name__ == '__main__':   
    test1()
    print()
    test2()
