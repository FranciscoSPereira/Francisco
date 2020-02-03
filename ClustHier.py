from BinaryTrees import BinaryTree
from MatrixNum import MatrixNum


class ClustHier:

    def __init__(self, matdists):
        self.matdists = matdists    # matriz de distancias (MatrixNum) - dimensão N*N
        
    def distance(self, tree1, tree2):
        '''
        calcula a distância entre duas árvores
        média da distância entre os elementos do cluster (método dado)
        '''
        c1 = tree1.getCluster()     # dá lista de valores(das folhas/seqs) da tree1
        c2 = tree2.getCluster()     # dá lista de valores da tree2
        sd = 0.0
        for i in range(len(c1)):
            for j in range(len(c2)):
                sd += self.matdists.getValue(c1[i], c2[j])  # soma a sd, os valores correspondentes na matriz de distâncias
        return sd / len(c1) * len(c2)
    
    def executeClustering(self):
        '''
        corre o algoritmo criando uma árvore a partir de uma matriz de distâncias
        '''
        trees = []
        tableDist = self.matdists.copy()       # falta-nos o método copy na MatrixNum (feito)
        for i in range(self.matdists.numRows()):
            t = BinaryTree(i)   # cria árvores(apenas folhas de value i), conforme o número de linhas da matriz de distâncias; 1ªseq ->value0, 2ªseq ->value1...
            trees.append(t)     # adiciona as folhas à lista trees
        for k in range(self.matdists.numRows(), 1, -1):     # começa no número de Rows(folhas) e vai até 2
            # choose minimum value from matrix not in diagonal
            mins = tableDist.minDistIndexes()       # falta método na MatrixNum, para retornar os índices correspondentes ao valor mais pequeno da matriz (feito)
            i = mins[0]
            j = mins[1]
            # pega nas folhas da lista tree, que correspondem ao valor mais pequeno na MD (?), e cria uma nova árvore a partir delas
            n = BinaryTree(-1, tableDist.getValue(i, j)/2.0, trees[i], trees[j])
            if k > 2:
                trees.pop(i)        # retira os valores das folhas que usamos para fazer o cluster novo (binary tree anterior)
                trees.pop(j)
                tableDist.removeRow(i)
                tableDist.removeRow(j)
                tableDist.removeCol(i)
                tableDist.removeCol(j)      # fuck, também não temos os métodos para remover colunas na MD (hmm yes we do)
                dists = []
                for x in range(len(trees)):
                    dists.append(self.distance(n, trees[x]))     # calcula distâncias entre a nova árvore e as restantes folhas, acrescenta na lista dists
                tableDist.addRow(dists)     # adicionamos linha à MD (matdist?) com essas distâncias novas
                cdists = []
                for y in range(len(dists)):     # o mesmo é feito para adicionar colunas
                    cdists.append(dists[y])
                cdists.append(0.0)
                tableDist.addCol(cdists)
                trees.append(n)         # adiciona tmb a nova árvore à lista trees
            else:
                # retorna árvore final
                return n
            tableDist.printmat()


# --------------------------\\-----------------------------

# TESTS

def test():
    m = MatrixNum(5,5)
    m.setValue(1, 0, 2)
    m.setValue(0, 1, 2)
    m.setValue(2, 0, 5)
    m.setValue(0, 2, 5)
    m.setValue(3, 0, 7)
    m.setValue(0, 3, 7)
    m.setValue(4, 0, 9)
    m.setValue(0, 4, 9)
    m.setValue(2, 1, 4)
    m.setValue(1, 2, 4)
    m.setValue(3, 1, 6)
    m.setValue(1, 3, 6)
    m.setValue(4, 1, 7)
    m.setValue(1, 4, 7)
    m.setValue(3, 2, 4)
    m.setValue(2, 3, 4)
    m.setValue(4, 2, 6)
    m.setValue(2, 4, 6)
    m.setValue(4, 3, 3)
    m.setValue(3, 4, 3)
    hc = ClustHier(m)
    arv = hc.executeClustering()
    arv.printtree()


if __name__ == '__main__': 
    test()
    
# Result for test
# Root:[1, 0, 2, 4, 3] Dist.: 3.25
#     Left:[1, 0, 2] Dist.: 2.25
#         Left:[1, 0] Dist.: 1.0
#             Left:[1] Dist.: 0
#             Right:[0] Dist.: 0
#         Right:[2] Dist.: 0
#     Right:[4, 3] Dist.: 1.5
#         Left:[4] Dist.: 0
#         Right:[3] Dist.: 0
