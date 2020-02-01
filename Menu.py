from MySeq import MySeq
from MyBlast import MyBlast 
from upgma import UPGMA
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix
from Bio import Entrez, SeqIO, Phylo
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Align.Applications import ClustalwCommandline
import re
import os

cwd=os.getcwd().strip('src')

class GestorAA:
   
    def __init__(self):
        self.dic={}
        self.seqs_alfabeto='ACTG'
        self.seqs_tipo='dna'
        self.name=None
        self.seq=None
        self.listaIds=[]

    def loadFromFile(self,nome):
        name=nome
        self.name=name
        self.listIds()
    
    def createNew(self,nome):
        name=nome
        self.name=name
    
    def addSeq(self,ident,seq):
        Id=ident
        sequen=seq
        if self.valSeq(seq):
            self.listaIds.append(Id)
            if Id not in self.dic.keys():
                self.dic[Id]={'Sequência':sequen}
                self.dic[Id]['Anotações']={}
            else:
                print('Identificador já utilizado')
        else:
            print('Sequência proteica inválida')
            
    def addAnnotation(self,tipo,anno,Id):
        if tipo not in self.dic[Id]['Anotações'].keys():
            self.dic[Id]['Anotações'][tipo]=[]
            self.dic[Id]['Anotações'][tipo].append(anno)
        else:
            self.dic[Id]['Anotações'][tipo].append(anno)
            
    def valSeq(self,seq):
        vs=MySeq(seq,'dna')
        if vs.validaER():
            return True
        else:
            return False
        
    def printbyid(self,Id):
        for key,val in self.dic[Id].items():
                print (key+':')
                if type(val)==dict:
                    self.printSubDic(val)
                elif type(val)==list:
                    i=0
                    for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
                else:
                    print(self.dic[Id][key])
                    print()     

    def printDic(self,dic=None):
        if dic==None:
            dic=self.dic
        else:
            dic=dic
        l=self.listaIds
        tid=None
        for x in range(len(l)):
            tid=l[x]
            print('>'+tid)
            for key,val in dic[tid].items():
                print (key+':')
                if type(val)==dict:
                    self.printSubDic(val)
                elif type(val)==list:
                    i=0
                    for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
                else:
                    print(dic[tid][key])
                    print()
                    
    def __printSubDic(self,subdic):
        for key, val in subdic.items():
            print(key+':')
            if type(val)==dict:
                self.printSubDic(val)
            elif type(val)==list:
                i=0
                for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
            else:
                print(subdic[key][val])
                print()
        
    def listIds(self):
        for key in self.dic.keys():
            self.listaIds.append(key)
                
    def editAnnotations(self,Id,tipo,nr,new):
        self.dic[Id]['Anotações'][tipo][nr]=new
        
    def editSeq(self,Id,new):
        self.dic[Id]['Sequência']=new
        
    def search_seqs(self,seq1):
        db1=[]
        for Id in self.listaIds:
            db1.append(self.dic[Id]['Sequência'])
        blast=MyBlast(db=db1)
        query=seq1
        res=blast.bestAlignment(query)
        bestmatch=db1[res[4]]
        print(bestmatch)
        for Id in self.listaIds:
            if self.dic[Id]['Sequência']==bestmatch:
                print('Melhor match encontrado na sequência com Id %s' %Id)
    
    def search_pattern(self,pat,Id):
        '''
        Procura ocorrências de padrões definidos pelo utilizador nas sequências da base de dados.
        '''
        seq=self.dic[Id]['Sequência']
        matches = re.finditer(pat,seq)
        res=[]
        for match in matches:
            res.append(match.span())
        if len(res)!=0:
            print('Foram encontradas %d ocorrências do padrão na sequência %s' % (len(res),Id))
            i=1
            for match in res:
                print('%d ª ocorrência, posição: '%i,match)
                i+=1
        else:
            print('Não foram encontradas ocorrências do padrâo na sequência %s' %Id)
            return None

    def search_db(self,pat):
        for Id in self.listaIds:
            self.search_pattern(pat,Id)

    def freq(self,sim=None,Id=None):
        if Id != None:
            seq=self.dic[Id]['Sequência']
            res=[]
            matches=re.finditer(sim,seq)
            for match in matches:
                res.append(match.span())
            freq=len(res)/len(seq)
            print('O simbolo apresenta uma frequência de %f'%freq)
        else:
            for Id in self.listaIds:
                seq=self.dic[Id]['Sequência']
                res=[]
                matches=re.finditer(sim,seq)
                for match in matches:
                    res.append(match.span())
                freq=len(res)/len(seq)
                print('O simbolo apresenta uma frequência de %f, na sequência %s'%(freq,Id))
        
        
    def arvore_filogenetica(self):
        '''
        Imprime uma árvore filogenética das diferentes sequências de DNA com as distâncias entre as folhas e nodos.
        '''
        bd=[]
        for Id in self.listaIds:
            bd.append(MySeq(self.dic[Id]['Sequência'],'dna'))
        sm = SubstMatrix()
        sm.createFromMatchPars(3, -1, self.seqs_alfabeto)
        alin=AlignSeq(sm,-1)
        up=UPGMA(bd,alin)
        arv=up.run()
        arv.printtree()
        
    def traducao(self,Id):
        seq=self.dic[Id]['Sequência']
        ms=MySeq(seq,self.seqs_tipo)
        res=ms.traduzSeq()
        print(res)
    
    def orfs(self,Id):
        seq=self.dic[Id]['Sequência']
        ms=MySeq(seq,self.seqs_tipo)
        res=ms.orfs()
        for x in res:
            x.printseq()
    
    def saveOrfs(self,Id):
        try:
            seq=self.dic[Id]['Sequência']
            ms=MySeq(seq,self.seqs_tipo)
            res=ms.orfs()
            file=open(cwd+'Orfs\\' + Id + '_orfs.txt','w')
            for orf in res:
                file.write(orf.seq+'\n')
            
            file.close()
        except:
            os.mkdir(cwd+'Orfs\\')
            self.saveOrfs(Id)
            
    def searchonline(self,Id,email):
        '''
        Permite acesso ao NCBI através da web e retirar sequências desta.
        '''
        try:
            Entrez.email=email
            if type(Id)==str:
                handle=Entrez.efetch(db='nucleotide',id=Id,rettype='fasta',retmode='text')
                record=SeqIO.read(handle,'fasta')
                handle.close()
                self.addSeq(Id,str(record.seq))
            elif type(Id)==list:
                for id_ in Id:
                    handle=Entrez.efetch(db='nucleotide',id=id_,rettype='fasta',retmode='text')
                    record=SeqIO.read(handle,'fasta')
                    self.addSeq(id_,str(record.seq))
                    handle.close()
        except:
            print('Erro ao executar procura online')
            
    def blast(self,seq,program='blastn',base='nt',evalue=0.05,align=50,matrix='BLOSUM62'):
        '''
        Função do Biopython que corre o blast de sequencias de dna contra a base de dados pretendida, tendo como parametros:
        o programa-blastn
        a matrix- BLOSUM62
        o e-value=0.001
        o nº de alinhamentos maximos=100
        '''
        
        result_handle=NCBIWWW.qblast(program,'nt',seq) 
        hit_id=[]
        blast_records = NCBIXML.parse(result_handle)
        e_value_thresh = 0.001
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        hit_id.append(alignment.hit_id)
                        print("****Alignment****")
                        print("sequence:", alignment.title)
                        print("length:", alignment.length)
                        print("e value:", hsp.expect)
                        print(hsp.query[0:75] + "...")
                        print(hsp.match[0:75] + "...")
                        print(hsp.sbjct[0:75] + "...")
        return hit_id
                     
    def addBlastSeqs(self,ids):
        '''
        Permite retirar sequências do NCBI.
        '''
        hit_id = list(set(ids))                  
        for ind in hit_id:               
            handle = Entrez.efetch(db ="nucleotide", rettype ='fasta', retmode ="fasta", id=ind)
            seq_record = SeqIO.read(handle, 'fasta')
            self.addSeq(ind,str(seq_record.seq))
            handle.close()
    
    def make_seqs(self,filename,Ids):
        '''
        Permite a escrita de sequências num ficheiro.                                                                               
        '''
        try:
            file=open(cwd+'AlinMuExe\\'+filename+'.fasta','w')
            for id_ in Ids:
                file.write('>'+id_+'\n')
                file.write(self.dic[id_]['Sequência']+'\n')
            file.close()
        except:
            os.mkdir(cwd+'AlinMuExe\\')
            self.make_seqs()
    
    def multiple_alignment(self,filename):
        '''
        Função que realiza o alinhamento múltiplo com as melhores sequências 
        do blast corrido anteriormente, e escolhidas pela função anterior 
        '''
        cmdline=ClustalwCommandline(r'C:\Program Files (x86)\ClustalW2\clustalw2',infile=cwd+'AlinMuExe\\'+filename+'.fasta')
        cmdline()
        
    def create_tree(self,filename):
        '''
        Função que cria a árvore filogenética com base no alinhamento múltiplo
        '''
        infile=(cwd+'AlinMuExe\\'+filename+'_seqs.dnd')
        tree = Phylo.read(infile, "newick")
        Phylo.draw_ascii(tree)
        return tree
    
    