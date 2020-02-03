from MySeq import MySeq
from MyBlast import MyBlast 
from upgma import UPGMA
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix
from Bio import Entrez, SeqIO, Phylo, AlignIO
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Align.Applications import ClustalwCommandline
from MultipleAlign import MultipleAlign
import pickle
import re
import os

cwd=os.getcwd().strip('src')

class GestorSeq:
   
    def __init__(self):
        self.dic={}
        self.seqs_alfabeto='ACTG'
        self.seqs_tipo='dna'
        self.nome=None
        self.seq=None
        self.listaids=[]

    def loadFromFile(self,name):
        nome=name
        self.nome=nome
        self.listaids()
    
    def criar(self,name):
        nome=name
        self.nome=nome
    
    def adicionar(self,identificador,seq):
        '''
        Adiciona sequências por identificador.
        '''
        Id=identificador
        sequence=seq
        if self.valida(seq):
            self.listaids.append(Id)
            if Id not in self.dic.keys():
                self.dic[Id]={'Sequência':sequence}
                self.dic[Id]['Anotações']={}
            else:
                print('Identificador já utilizado! Tente outro.')
        else:
            print('Sequência de proteínas inválida!')
            
    def add_anotacoes(self,tipo,anno,Id):
        if tipo not in self.dic[Id]['Anotações'].keys():
            self.dic[Id]['Anotações'][tipo]=[]
            self.dic[Id]['Anotações'][tipo].append(anno)
        else:
            self.dic[Id]['Anotações'][tipo].append(anno)
            
    def valida(self,seq):
        vs=MySeq(seq,'dna')
        if vs.validaER():
            return True
        else:
            return False
        
    def printbyid(self,Id):
        for key,val in self.dic[Id].items():
                print (key+':')
                if type(val)==dict:
                    self.printsubdicionario(val)
                elif type(val)==list:
                    i=0
                    for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
                else:
                    print(self.dic[Id][key])
                    print()     

    def dicionario(self,dic=None):
        if dic==None:
            dic=self.dic
        else:
            dic=dic
        lista=self.listaIds
        tid=None
        for x in range(len(lista)):
            tid=lista[x]
            print('>'+tid)
            for key,val in dic[tid].items():
                print (key+':')
                if type(val)==dict:
                    self.subdicionario(val)
                elif type(val)==list:
                    i=0
                    for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
                else:
                    print(dic[tid][key])
                    print()
                    
    def subdicionario(self,subdic):
        for key, val in subdic.items():
            print(key+':')
            if type(val)==dict:
                self.subdicionario(val)
            elif type(val)==list:
                i=0
                for x in range(len(val)):
                        print(str(i)+':'+val[x])
                        i+=1
            else:
                print(subdic[key][val])
                print()
        
    def list_id(self):
        for key in self.dic.keys():
            self.list_id.append(key)
                
    def edit_anotacoes(self,Id,tipo,nr,new):
        self.dic[Id]['Anotações'][tipo][nr]=new
        
    def edit_sequencias(self,Id,new):
        self.dic[Id]['Sequência']=new
        
    def multialign(self,filename):
        alignments = AlignIO.parse(filename, "fasta")
        for alignment in alignments:
            print ("Alinhamento:", alignment)
        alignments = list(AlignIO.parse (filename, "fasta"))
        last_align = alignments[-1]
        first_align = alignments[0]
        
    def procura_seqs(self,seq1):
        db1=[]
        for Id in self.listaids:
            db1.append(self.dic[Id]['Sequência'])
        blast=MyBlast(db=db1)
        query=seq1
        res=blast.bestAlignment(query)
        bestmatch=db1[res[4]]
        print(bestmatch)
        for Id in self.listaIds:
            if self.dic[Id]['Sequência']==bestmatch:
                print('Melhor match encontrado na sequência com Id %s' %Id)
    
    def procura_pad(self,pad,Id):
        '''
        Procura ocorrências de padrões definidos pelo utilizador nas sequências por identificador.
        '''
        seq=self.dic[Id]['Sequência']
        matches = re.finditer(pad,seq)
        res=[]
        for match in matches:
            res.append(match.span())
        if len(res)!=0:
            print('Foram encontradas %d ocorrências do padrão na sequência %s' % (len(res),Id))
            i=1
            for x in res:
                print('%d ª ocorrência, posição: '%i,x)
                i+=1
        else:
            print('Não foram encontradas ocorrências do padrâo na sequência %s' %Id)
            return None

    def search_db(self,pad):
        '''
        Procura ocorrências de padrões de sequências numa base de dados.
        '''
        for Id in self.listaids:
            self.search_pattern(pad,Id)

    def freq_simbolos(self,s=None,Id=None):
        '''
        Determina a frequência de símbolos de uma sequência selecionada.
        '''
        if Id != None:
            seq=self.dic[Id]['Sequência']
            res=[]
            matches=re.finditer(s,seq)
            for match in matches:
                res.append(match.span())
            freq=len(res)/len(seq)
            print('O simbolo apresenta uma frequência de %f'%freq)
        else:
            for Id in self.listaids:
                seq=self.dic[Id]['Sequência']
                res=[]
                matches=re.finditer(s,seq)
                for match in matches:
                    res.append(match.span())
                freq=len(res)/len(seq)
                print('O simbolo apresenta uma frequência de %f, na sequência %s'%(freq,Id))
        
        
    def arvore_filogenetica(self):
        '''
        Imprime uma árvore filogenética das diferentes sequências de DNA com as distâncias entre as folhas e nodos.
        '''
        bd=[]
        for Id in self.listaids:
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
    
    def write_seqs(self,filename,ids):
        '''
        Permite a escrita de sequências num ficheiro.                                                                               
        '''
        try:
            file=open(cwd+'AlinMuExe\\'+filename+'.fasta','w')
            for Id in ids:
                file.write('>'+Id+'\n')
                file.write(self.dic[Id]['Sequência']+'\n')
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

    def RemindSeq(self):
        for Id in self.listaids:
            print(Id + '->' + self.dic[Id]["Sequência"])
    
    def import_(self,filename):
        fasta_seq=SeqIO.parse(open(filename),"fasta")
        for seq in fasta_seq:
            if self.adicionar(seq.id,str(seq.seq).upper()):
                self.add_anotacoes("nome",seq.nome,seq.id)
                
    def guardar_dic(self):
        file=open(self.name+'.txt','wb')
        pickle.dump(self.dic,file)
        file.close()