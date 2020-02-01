
from MatrixNum import MatrixNum
from ClustHier import ClustHier
from MySeq import MySeq
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix


class UPGMA:

    def __init__(self, seqs, alseq):
        '''
        recebe um conjunto de sequências e um objeto AlignSeq com os parâmetros para realizar alinhamentos
        '''
        self.seqs = seqs
        self.alseq = alseq      # objeto AlignSeq com os parâmetros para realizar alinhamentos(sm e g)
        self.matdist = MatrixNum(len(seqs), len(seqs))      # template para a matriz de distâncias
        self.criaMatDists()     # resultado do método abaixo, completa automaticamente a matdist ao criar uma instância UPGMA
        
    def criaMatDists(self):
        '''
        calcula a matriz de distâncias (objeto MatrixNum) fazendo alinhamento global de cada par de sequências
        distância = nº caracteres diferentes entre as sequências (após alinhamento)
        a matriz fica guardada como atributo da classe
        '''
        for i in range(len(self.seqs)): 
            self.matdist.setValue(i, i, 0.0)        # cria diagonal com zeros
        
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):    # esta forma de iteração impede que se realizem alinhamentos repetidos
                # ex: alinhar seq0 com seq1 e na próxima iteração i, evita que haja alinhamento entre seq1 e seq0
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needlemanWunsch(s1, s2)
                alin = self.alseq.recoverAlignment()        # alin --> objeto MyAlign
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)    # col --> lista de nucleotidos/proteinas na coluna k do alinhamento
                    if col[0] != col[1]:
                        ncd += 1    # ncd --> nª de carateres que diferem entre as sequencias do alinhamento, vai corresponder à distância posta na matriz de distâncias
                self.matdist.setValue(i, j, ncd)    # valores acima da diagonal de zeros
                self.matdist.setValue(j, i, ncd)    # simetria abaixo da diagonal de zeros
                # ahhh daí não termos deixado fazer alinhamentos repetidos nos 2 ciclos for de cima, fascinante.
                    
    def run(self):
        '''
        aplica método executeClustering de ClustHier
        0 resultado é a árvore filogenética
        '''
        ch = ClustHier(self.matdist)
        arv = ch.executeClustering()
        # retorna a árvore(cluster) final, que podemos posteriormente imprimir e ver a filogenia
        return arv
        

# --------------------------\\-----------------------------

# TESTS

def test1():                # o resultado da diferente
    seq1 = MySeq('ATAGCGAT')    
    seq2 = MySeq('ATAGGCCT')    
    seq3 = MySeq('CTAGGCCC')    
    sm = SubstMatrix()    
    sm.createFromMatchPars(3, -1, 'ACGT')    
    alseq = AlignSeq(sm, -8)    
    up = UPGMA([seq1, seq2, seq3], alseq)
    arv = up.run()    
    arv.printtree() 
    
  
# Root:[2, 1, 0] Dist.: 2.0
#        Left:[2, 1] Dist.: 1.0
#                Left:[2] Dist.: 0
#                Right:[1] Dist.: 0
#        Right:[0] Dist.: 0    
    

def test2():        # o resultado da diferente
    seq1 = MySeq('AGTAG')
    seq2 = MySeq('TATA')
    seq3 = MySeq('TAT')
    seq4 = MySeq('TAGTG')    
    sm = SubstMatrix()
    sm.createFromMatchPars(0, -1, 'ACGT')
    alseq = AlignSeq(sm, -1)
    up = UPGMA([seq1, seq2, seq3, seq4], alseq)
    up.matdist.printmat()
    arv = up.run()
    arv.printtree()

# [0, 3, 3, 2]
# [3, 0, 1, 2]
# [3, 1, 0, 2]
# [2, 2, 2, 0]
#
# Root:[3, 0, 2, 1] Dist.: 1.25
#        Left:[3, 0] Dist.: 1.0
#                Left:[3] Dist.: 0
#                Right:[0] Dist.: 0
#        Right:[2, 1] Dist.: 0.5
#                Left:[2] Dist.: 0
#                Right:[1] Dist.: 0


def test3():
    seq1 = MySeq('AGTAG')
    seq2 = MySeq('TATA')
    seq3 = MySeq('TAT')
    seq4 = MySeq('TAGTG')
    sm = SubstMatrix()
    sm.loadFromFile('blosum62.mat', '\t')
    alseq = AlignSeq(sm, -1)
    up = UPGMA([seq1, seq2, seq3, seq4], alseq)
    up.matdist.printmat()
    arv = up.run()
    arv.printtree()


if __name__ == '__main__': 
    test1()
    print()
    test2()
    print()
    test3()
