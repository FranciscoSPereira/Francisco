
class MyAlign:
    '''
    serve para guardar e manipular alinhamentos de sequências (duas ou mais sequências)
    '''

    def __init__(self, lseqs, tipo='protein'):
        self.listseqs = lseqs               # seria a list res (AlignSeq), resultante do traceback
        self.tipo = tipo
    
    def __len__(self):  # number of columns
        return len(self.listseqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) == 2:
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int:
            return self.listseqs[n]
        return None
    
    def __str__(self):
        res = ''
        for seq in self.listseqs:
            res += '\n' + seq 
        return res
    
    def numSeqs(self):
        return len(self.listseqs)
    
    def column(self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res

    def consensus(self):
        '''
        alinhamento representado com lista de strings
        retorna o caractere mais comum em cada coluna, ignorando gaps
        '''
        cons = ''
        for i in range(len(self)):         # self significa len(self.listseqs[0]) - que é o tamanho da seq
            count = {}
            for j in range(len(self.listseqs)):
                c = self.listseqs[j][i]
                if c in count:
                    count[c] = count[c] + 1
                else:
                    count[c] = 1
            maximum = 0
            cmax = None
            for ke in count.keys():
                if ke != '-' and count[ke] > maximum:
                    maximum = count[ke]
                    cmax = ke
            cons = cons + cmax
            #print(count)
        return cons

    def consensus_column(self, indice_col):
        listacar = self.column(indice_col)
        res = True
        i = 1
        while res and i < len(listacar):
            if listacar[i] != listacar[0]:
                res = False
            else:
                i = i+1
        return res

    def consensus_list(self):
        res = []
        for c in range(len(self)):              # so self ?
            if self.consensus_column(c):
                res.append(c)
        return res


# --------------------------\\-----------------------------

# TESTS

if __name__ == '__main__': 
    alig = MyAlign(['ATGA-A', 'AA-AT-'], 'dna')
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1, 1])
    print(alig[0, 2])
    print(alig.consensus())


# Results
# ATGA-A
# AA-AT-
# 6
# ['G', '-']
# A
# G
