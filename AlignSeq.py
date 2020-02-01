from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix
from MatrixNum import MatrixNum
# from matrizes import maxMat


class AlignSeq:
    '''
    esta classe implementa algoritmos de programação dinamica (PD) para alinhamentos
    de pares de sequencias tendo como atributos:
    - Scoring matrix (sm): objeto da classe SubstMatrix
    - Penalizaçao de gaps (g)
    - Seq a alinhar
    - Matrizes S e T
    '''

    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None
        self.seq2 = None
        
    def scorePos(self, c1, c2):
        '''
        Score de uma posição: par de letras
        '''
        if c1 == '-' or c2 == '-':
            return self.g
        else:
            return self.sm[c1, c2]
        
    def scoreAlin(self, alin):
        '''
        score de um alinhamento (objecto MyAlign)
        '''
        res = 0
        for i in range(len(alin)):
            res += self.scorePos(alin[0][i], alin[1][i])
        return res
    
    def needlemanWunsch(self, seq1, seq2):
        '''
        corre algoritmo de Needleman Wunsch, atualizando as matrizes S e T
        retorna score do melhor alinhamento global
        matriz s -> elementos da seq1 ao colocados nas linhas e os da seq2 nas colunas
        matriz T -> serve para fazer o traceback
          |
          |__  1 -> diagonal
          |__  2 -> vertical
          |__  3 -> horizontal
        '''
        if seq1.tipo != seq2.tipo:
            return None
        self.S = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.T = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.seq1 = seq1
        self.seq2 = seq2
        # nunca toca no ponto (0,0)
        for j in range(1, len(seq2)+1):
            self.S[0, j] = self.g * j
            self.T[0, j] = 3
        for i in range(1, len(seq1)+1):
            self.S[i, 0] = self.g * i
            self.T[i, 0] = 2
        for i in range(0, len(seq1)):           # pq o zero?
            for j in range(len(seq2)):
                s1 = self.S[i, j] + self.scorePos(seq1[i], seq2[j])     # diagonal
                s2 = self.S[i, j+1] + self.g            # gaps
                s3 = self.S[i+1, j] + self.g            # gaps
                self.S[i+1, j+1] = max(s1, s2, s3)      # ve qual o maior score para por no ponto
                self.T[i+1, j+1] = max3t(s1, s2, s3)
        # faz return das coordenadas do score perfeito (ultimo elem da matrix)
        return self.S[len(seq1), len(seq2)]
    
    def recoverAlignment(self):
        '''
        retorna melhor alinhamento global
        faz traceback para ter o alinhamento otimo
        '''
        res = ['', '']          # lista com as duas seq do alinhamento
        i = len(self.seq1)
        j = len(self.seq2)
        while i > 0 or j > 0:
            if self.T[i, j] == 1:
                res[0] = self.seq1[i-1] + res[0]        # o i é o len e nao o indice, dai o -1
                res[1] = self.seq2[j-1] + res[1]        # o j é o len e nao o indice, dai o -1
                i -= 1
                j -= 1
            elif self.T[i, j] == 3:
                res[0] = '-' + res[0]
                res[1] = self.seq2[j-1] + res[1] 
                j -= 1
            else:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = '-' + res[1]
                i -= 1
        return MyAlign(res, self.seq1.tipo)

    def recoverAlignment_with_ties(self):
        '''
        explicar...
        '''
        i = len(self.seq1)
        j = len(self.seq2)
        alins = [['', '', i, j]]
        res = []
        while alins:
            al = alins.pop(0)
            i = al[2]
            j = al[3]
            if i == 0 and j == 0:
                res.append(al[:2])
            else:
                for t in self.T[i][j]:
                    p = []
                    if t == 1:
                        p.append(self.seq1[i-1] + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i-1)
                        p.append(j-1)
                    elif t == 3:
                        p.append('-' + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i)
                        p.append(j-1)
                    else:
                        p.append(self.seq1[i-1] + al[0])
                        p.append('-' + al[1])
                        p.append(i-1)
                        p.append(j)
                    alins.append(p)
        return res
 
    def smithWaterman(self, seq1, seq2):
        '''
        corre Smith Waterman, atualizando matrizes S e T
        retorna score do melhor alinhamento local
        '''
        if seq1.tipo != seq2.tipo:   # se as seq forem de tipos diferentes da None (rna, dna, proteina)
            return None
        else:
            self.S = MatrixNum(len(seq1) + 1, len(seq2) + 1)
            self.T = MatrixNum(len(seq1) + 1, len(seq2) + 1)
            self.seq1 = seq1
            self.seq2 = seq2

            maxscore = 0
            # nunca toca no ponto (0,0)
            for i in range(1, len(seq1) + 1):
                self.S[i, 0] = 0
                self.T[i, 0] = 0
            for j in range(1, len(seq2) + 1):
                self.S[0, j] = 0
                self.T[0, j] = 0
            for i in range(0, len(self.seq1)):
                for j in range(len(self.seq2)):
                    s1 = self.S[i, j] + self.scorePos(seq1[i], seq2[j])
                    s2 = self.S[i, j + 1] + self.g
                    s3 = self.S[i + 1, j] + self.g
                    b = max(s1, s2, s3)
                    if b <= 0:
                        self.S[i + 1, j + 1] = 0
                        self.T[i + 1, j + 1] = 0
                    else:
                        self.S[i + 1, j + 1] = b
                        self.T[i + 1, j + 1] = max3t(s1, s2, s3)
                        if b > maxscore:
                            maxscore = b
            return maxscore

    def smithWatermanTies (self, seq1, seq2, ties=False):
        '''
        explicar...
        onde se aplica? como se sabe que e um empate?
        '''
        if seq1.tipo != seq2.tipo:
            return None
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        maxscore = 0
        for j in range(1, len(seq2)+1):
            self.S[0].append(0)
            if ties:
                self.T[0].append([0])
            else:
                self.T[0].append(0)
        for i in range(1, len(seq1)+1):
            self.S.append([0])
            if ties:
                self.T.append([[0]])
            else:
                self.T.append([0])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos(seq1[i], seq2[j])
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                b = max(s1, s2, s3)
                if b <= 0:
                    self.S[i+1].append(0)
                    self.T[i+1].append(0)
                else:
                    self.S[i+1].append(b)
                    if ties:
                        self.T[i+1].append(max3t_with_ties(s1, s2, s3))
                    else:
                        self.T[i+1].append(max3t(s1, s2, s3))
                    if b > maxscore:
                        maxscore = b
        return maxscore

    def recoverLocalAlignment(self):
        '''
        retorna melhor alinhamento local
        '''
        # from matrizes import maxMat ???? donde vem isto?
        res = ['', '']
        i, j = self.S.maxIndexes()
        while self.T[i, j] > 0:
            if self.T[i, j] == 1:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = self.seq2[j - 1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i, j] == 2:
                res[0] = self.seq1[i - 1] + res[0]
                res[1] = '-' + res[1]
                i -= 1
            elif self.T[i, j] == 3:
                res[0] = '-' + res[0]
                res[1] = self.seq2[j - 1] + res[0]
                j -= 1
        return MyAlign(res, self.seq1.tipo)

    def recoverAlignLocal_with_ties (self):
        '''
        explicar...
        '''
        maxval = self.S[0][0]
        maxtups = []
        for i in range(0, len(self.S)):
            for j in range(0, len(self.S[i])):
                if self.S[i][j] > maxval:
                    maxval = self.S[i][j]
                    maxtups = [(i, j)]
                elif self.S[i][j] == maxval:
                    maxtups.append((i, j))
        alins = []
        for (i, j) in maxtups:
            alins.append(['', '', i, j])
        res = []
        while alins:
            al = alins.pop(0)
            i = al[2]
            j = al[3]
            if (i == 0 and j == 0) or (0 in self.T[i][j]):
                res.append(al[:2])
            else:
                for t in self.T[i][j]:
                    p = []
                    if t == 1:
                        p.append(self.seq1[i-1] + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i-1)
                        p.append(j-1)
                    elif t == 3:
                        p.append('-' + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i)
                        p.append(j-1)
                    else:
                        p.append(self.seq1[i-1] + al[0])
                        p.append('-' + al[1])
                        p.append(i-1)
                        p.append(j)
                    alins.append(p)
        return res

    def alignQuery(self, query, setSeqs):
        '''
        explicar...
        '''
        bestScore = -1
        bestAlin = None
        for seq in setSeqs:
            sc = self.smithWaterman(query, seq)
            if sc > bestScore:
                bestScore = sc
                bestAlin = self.recoverLocalAlignment()
        return bestAlin, bestScore

    def alignSeqWithSet(self, query, setSeqs):
        '''
        explicar...
        '''
        scores = []
        res = [[0, 0]] * len(setSeqs)
        k = 0
        for s in setSeqs:
            scores[k] = self.smithWaterman(query, s)
            k = k + 1
        for i in range(len(setSeqs)):
            j = 0
            while j < i and scores[res[j][0]] > scores[i]:
                j = j + 1
            for k in range(i, j, -1):
                res[k][0] = res[k - 1][0]
            res[j][0] = i
        for i in range(len(setSeqs)):
            res[i][1] = scores[res[i][0]]
        return res

    def alignSets(self, c1, c2):
        '''
        explicar...
        '''
        res = []
        for i in range(len(c1)):
            maxscore = -1
            bestseq = None
            for j in range(len(c2)):
                score = self.smithWaterman(c1[i], c2[j])
                if score > maxscore:
                    maxscore = score
                    bestseq = c2[j]
            res.append(bestseq)
        return res

    def alignProtDNA(self, prot, dna):
        '''
        explicar...
        '''
        res = None
        orfs = dna.orfs()
        prots = []
        for i in range(len(orfs)):
            prots.append(orfs[i].extraiTodasProts())
        mx = 0
        for p in prots:
            score = self.smithWaterman(prot.seq, p.seq)
            if score > mx:
                mx = score
                res = self.recoverLocalAlignment()          # ?? estava self.getLocalAlignment() - nao existe
        return res


# funcoes externa a classe
def max3t(s1, s2, s3):
    if s1 > s2:
        if s1 > s3:
            return 1
        else:
            return 3
    else:
        if s2 > s3:
            return 2
        else:
            return 3


def max3t_with_ties(v1, v2, v3):
    if v1 > v2:
        if v1 > v3:
            return [1]
        elif v1 == v3:
            return [1, 3]
        else:
            return [3]
    elif v1 == v2:
        if v1 > v3:
            return [1, 2]
        elif v1 == v3:
            return [1, 2, 3]
        else:
            return [3]
    else:
        if v2 > v3:
            return [2]
        elif v2 == v3:
            return [2, 3]
        else:
            return [3]


def printMat(mat):
    for i in range(0, len(mat)):
        print(mat[i])


def lcs(seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(1, 0, seq1.alfabeto())
    aseq = AlignSeq(sm, 0)
    aseq.needlemanWunsch(seq1, seq2)
    alin = aseq.recoverAlignment()
    sizeal = len(alin[0])
    lcs = ''
    for i in range(sizeal):
        if alin[0][i] != '-' and alin[0][i] == alin[1][i]:
            lcs += alin[0][i]
    return lcs


def edit_distance(seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(0, -1, seq1.alfabeto())
    aseq = AlignSeq(sm, -1)
    sc = aseq.needlemanWunsch(seq1, seq2)
    return -sc


# --------------------------\\-----------------------------

# TESTS


def testScoresProt():
    s1 = MySeq('LGPS-GCASGIWTKSA', 'protein')
    s2 = MySeq('TGPSGG--SRIWTKSG', 'protein')
    alin = MyAlign([s1, s2], 'protein')
    
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    
    aseq = AlignSeq(sm, -8) 
    score = aseq.scoreAlin(alin)
    
    print('Score', score)
 
# 29


def testScoresDNA():
    s1 = MySeq('ATGA-AGGT', 'dna')
    s2 = MySeq('A-GAGAGGC', 'dna')
    alin = MyAlign([s1, s2], 'dna')
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1,0, 'ACGT')

    aseq = AlignSeq(sm, 0) 
    score = aseq.scoreAlin(alin)

    print('Score', score)

# 6


def test1():
    seq1 = MySeq('PHSWG','protein')
    seq2 = MySeq('HGWAG','protein')
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    alin = AlignSeq(sm, -8)

    print('Alinhamento global:')

    print('Score:', alin.needlemanWunsch(seq1, seq2))
    print('Matriz S:')
    alin.S.printmat()
    print('Matriz T:')
    alin.T.printmat()
    print('Alinhamento otimo:')
    print(alin.recoverAlignment())

    # ------------------------------

    print('Alinhamento local:')
    
    print('Score:', alin.smithWaterman(seq1, seq2))
    print('Matriz S:')
    alin.S.printmat()
    print('Matriz T:')
    alin.T.printmat()
    print('Alinhamento otimo:')
    print(alin.recoverLocalAlignment())


# Alinhamento global
# Score: 9

# Matriz S:
# [0,    -8, -16, -24, -32, -40]
# [-8,   -2, -10, -18, -25, -33]
# [-16,   0,  -4, -12, -20, -27]
# [-24,  -8,   0,  -7, -11, -19]
# [-32, -16,  -8,  11,   3,  -5]
# [-40, -24, -10,   3,  11,   9]

# Matriz T:
# [0, 3, 3, 3, 3, 3]
# [2, 1, 3, 3, 1, 3]
# [2, 1, 1, 3, 3, 1]
# [2, 2, 1, 1, 1, 3]
# [2, 2, 2, 1, 3, 3]
# [2, 2, 1, 2, 1, 1]

# Alinhamento otimo:
# PHSW-G
# -HGWAG

# Alinhamento local
# Score: 19

# Matriz S:
# [0, 0, 0,  0,  0,  0]
# [0, 0, 0,  0,  0,  0]
# [0, 8, 0,  0,  0,  0]
# [0, 0, 8,  0,  1,  0]
# [0, 0, 0, 19, 11,  3]
# [0, 0, 6, 11, 19, 17]

# Matriz T:
# [0, 0, 0, 0, 0, 0]
# [0, 0, 0, 0, 0, 0]
# [0, 1, 0, 0, 0, 0]
# [0, 0, 1, 0, 1, 0]
# [0, 0, 0, 1, 3, 3]
# [0, 0, 1, 2, 1, 1]

# Alinhamento otimo:
# HSW
# HGW


if __name__ == '__main__':   
    testScoresProt()
    print()
    testScoresDNA()
    print()
    test1()
