
from MyAlign import MyAlign
from MySeq import MySeq
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix


class MultipleAlign(object):
    '''
    AM - procedimento de comparação de mais do que duas sequências, procurando
    séries de caracteres individuais que se encontrem na mesma ordem

    Permite implementar o algoritmo de alinhamento múltiplo referido

    Tem como atributos lista de seq a alinhar (lista de objetos da classe MySeq), parametros do alinhamento (dados por
    um objeto da classe AlignSeq - inclui a matriz de substituiçao e o valor de penalizaçao de gaps)
    '''

    def __init__(self, seqs, alignseq):
        self.seqs = seqs                # lista de seqs a alinhar - lista de objetos da classe MySeq
        self.alignpars = alignseq       # parametros dos alinhamentos dados por objeto AlignSeq - lista de 2 seqs
    
    def numSeqs(self):
        return len(self.seqs)
    
    def scoreColumn(self, charsCol):
        sc = 0
        # ...
        return sc
     
    def scoreSP (self, alinhamento):
        sp = 0
        # ...
        return sp
    
    def addSeqAlignment(self, alignment, seq):
        '''
        adiciona uma nova sequência a um alinhamento existente (fazendo o consenso do alinhamento)
        cria alinhamento multiplo inculindo a nova seq
        '''
        res = []
        for i in range(len(alignment.listseqs)+1):  # cria alinhamento vazio
            res.append('')
        cons = MySeq(alignment.consensus(), alignment.tipo)  # consenso do alinhamento anterior --- .consensus (MyAlign)
        self.alignpars.needlemanWunsch(cons, seq)           # alinhamento global
        align2 = self.alignpars.recoverAlignment()  # alinha consenso com nova seq --- sera um objeto MyAlign
        orig = 0
        for i in range(len(align2)):  # cada coluna
            if align2[0][i] == '-':                         # align2[0, i] ==> align2[0][i]  --- é o mesmo?
                for k in range(len(alignment.listseqs)):
                    res[k] += '-'                           # once a gap, always a gap
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k][orig]
                orig += 1
        res[len(alignment.listseqs)] = align2.listseqs[1]      # acrescenta a ultima pos da lista res o align da seq
        return MyAlign(res, alignment.tipo)

    def alignConsensus(self):
        '''
        implementa o algoritmo completo de alinhamento
        alinha duas primeiras seq da lista e depois adiciona uma a uma as restantes seq com a funçao anterior
        '''
        self.alignpars.needlemanWunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recoverAlignment()
        for i in range(2, len(self.seqs)):
            res = self.addSeqAlignment(res, self.seqs[i])
        return res


# funçao externa a classe
def printMat(mat):
    for i in range(0, len(mat)):
        print(mat[i])


# --------------------------\\-----------------------------

# TESTS

def test():  
    s1 = MySeq('PHWAS', 'protein')
    s2 = MySeq('HWASW', 'protein')
    s3 = MySeq('HPHWA', 'protein')
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    aseq = AlignSeq(sm, -8)
    ma = MultipleAlign([s1, s2, s3], aseq)
    alinm = ma.alignConsensus()
    print(alinm)
    print(ma.scoreSP(alinm))        # por fazer

# -PHWAS-
# --HWASW
# HPHWA--


def test2():  
    s1 = MySeq('SWSSKLMKKIM', 'protein')
    s2 = MySeq('SYSLMKLKSWK', 'protein')
    s3 = MySeq('SWSSLMKLILS', 'protein')
    s4 = MySeq('SWSLMKLISSW', 'protein')
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    aseq = AlignSeq(sm, -8)
    ma = MultipleAlign([s1, s2, s3, s4], aseq)
    alinm = ma.alignConsensus()
    print(alinm) 
    print(ma.scoreSP(alinm))

# SWSSKLMKKIM--
# SYSLMKLKSWK--
# SWSS-LMKLILS-
# SWS--LMKLISSW
# 30


def testBiom():
    # s1 = MySeq("NRCD","protein")
    # s2 = MySeq("RADC","protein")
    # s3 = MySeq("NRDC","protein")
    
    s1 = MySeq('RADN', 'protein')
    s2 = MySeq('ARCD', 'protein')
    s3 = MySeq('ARDN', 'protein')

    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    aseq = AlignSeq(sm, -2)
    ma = MultipleAlign([s1, s2, s3], aseq)
    alin = MyAlign(['-RADN', 'ARCD-', 'AR-DN'], 'protein')
    print(alin)
    print(ma.scoreSP(alin))

    alinm = ma.alignConsensus()
    print(alinm) 
    
    print(ma.highQuality(alin))


def testAASB():
    s1 = MySeq('SWSSKLMKKIM', 'protein')
    s2 = MySeq('SYSLMKLKSWK', 'protein')
    s3 = MySeq('SWSSLMKLILS', 'protein')
    s4 = MySeq('SWSLMKLISSW', 'protein')
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    aseq = AlignSeq(sm, -2)
    ma = MultipleAlign([s1, s2, s3, s4], aseq)

    alin = MyAlign(['SWSSKLMKKIM-', 'SYSLMKLKSWK-', 'SWS-SLMKLILS', 'SWSL-MKLISSW'], 'protein')
    print(alin)
    print(ma.scoreSP(alin))
    
    x = ma.improveAlignment(alin)
    print(x)
    print(ma.scoreSP(x))


def testAASB2017 ():
    s1 = MySeq('ACATATCAT')
    s2 = MySeq('ACTAGATCT')
    s3 = MySeq('AGATATTAG')
    s4 = MySeq('GCATCGATT')
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1, -1, 'ACGT')
    aseq = AlignSeq(sm, -1)
    
    print(aseq.needlemanWunsch(s1, s2, True))
    printMat(aseq.S)
    printMat(aseq.T)
    print(aseq.recoverAlignment_with_ties())
    
    ma = MultipleAlign([s1, s2, s3, s4], aseq)
    al = ma.alignConsensus()
    print(al)
    al.distPairs()


def testBook ():
    s1 = MySeq('ATAGC')
    s2 = MySeq('AACC')
    s3 = MySeq('ATGAC')
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1, -1, 'ACGT')
    aseq = AlignSeq(sm, -1)
    ma = MultipleAlign([s1, s2, s3], aseq)
    al = ma.alignConsensus()
    print(al)
    

if __name__ == '__main__':
    test()
    test2()
    #testBiom()
    #testAASB2017()
    #test2()
    testBook()
