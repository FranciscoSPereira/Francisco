from MySeq import MySeq
from cmd import *
import re
from Menu import GestorAA
from ABC import Seq

class GestorAAShell(Cmd):
    intro = '''Interpretador de comandos para o gestor de sequências proteicas.\n
    1. Criar sequências: inserir criar
    2. Adicionar sequências: inserir addSeq
    3. Procurar por identificador: inserir searchid
    4.
    5.
    6. 
    7. Sair do programa: inserir sair
    '''
    def do_criar(self,arg):
        try:
            nome=arg
            if re.match('^[0-9]*$',nome):
                eng.createNew(nome)
        except:
            print('Erro ao criar')
    
    def do_addSeq(self,arg):
        try:
            seq= input("Seq:")
            a.valida(seq)

        except:
            print('Erro ao adicionar sequência')

    def do_searchid(self,arg):
        try:
            if arg in eng.listIds():
                eng.printDic()
            else:
                print('Identificador não encontrado no gestor')
        except:
            print('Erro ao procurar por identificador')
            
    def do_editar(self,arg):
        pass
    
    def do_print(self,arg):
        eng.printDic()
        
    def do_getSeqsOnline(self, arg):
        pass
        
    def do_sair(self, arg):
        "Sair do programa"
        print('Obrigado por ter utilizado o gestor de sequências!')
        global janela
        if janela is not None:
            del janela
        return True

def adicionar():
    seq_dna = input("Sequencia:")
    s1 = MySeq(seq_dna)
    s1.printseq()
    
def teste_ficheiro():
    ficheiro = input("Nome ficheiro: ")
    #por mensagem de erro caso o ficheiro não existir
    
    with open(f"{ficheiro}.txt", 'r') as readf:
        readf.readline()
        seq_ficheiro = readf.read()
        seq_ficheiro2 = seq_ficheiro.replace(" ", "")
        seq_ficheiro3 = seq_ficheiro2.replace("\n", "")
    
    s2 = MySeq(seq_ficheiro3)
    db.adicionar(s2)
    s2.guardar(seq_ficheiro3)
    s2.get_guardar()
    s2.printseq()
    print(s2.valida())

if __name__ == '__main__':
    eng = GestorAA()
    a=Seq()
    janela = None
    sh = GestorAAShell()
    sh.cmdloop()
    