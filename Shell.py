# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 17:53:13 2020

@author: Francisco
"""
from MySeq import MySeq
import pickle
from cmd import *
import sys
from Menu import GestorSeq

class Shell(Cmd):
    introd= """ Bem vindo ao gestor de sequências de DNA. Para iniciar introduzir "entrar". """
    prompt= "Gestor de sequências -> "
    
    def do_menu(self, arg=None):
        while True:
            print("""
            ======== Teste das funções trabalho algoritmos ============
            1. Criar uma base de dados
            2. Importar uma base de dados
            3. Adicionar sequências manualmente
            4. Adicionar anotações
            5. Editar sequências
            6. Imprimir validade da sequência
            7. Escrever sequências num ficheiro
            8. Procurar a sequência mais semelhante na base de dados
            9. Procurar um padrão
            10. Determinar frequência de símbolos
            11. Procurar padrões em bases de dados
            12. Imprimir uma árvore filogenética
            13. Fazer a tradução da sequência
            14. Imprimir as ORFs
            15. Fazer alinhamento múltiplo
            16. Sair
            """)
            op = input("Escolha um teste: ")
            if op == "1":
                self.chamaCriar()
            elif op == "2":
                self.chamaFasta()
            elif op == "3":
                self.chamaAdicionar()
            elif op == "4":
                self.chamaAnotacoes()
            elif op == "5":
                self.chamaEditar()
            elif op == "6":
                self.chamaValida()
            elif op == "7":
                self.chamaEscritaSeqs()
            elif op == "8":
                self.chamaProcuraSeqs()
            elif op == "9":
                self.chamaProcuraPad()
            elif op == "10":
                self.chamaFreqSimb()
            elif op == "11":
                self.chamaSearchDB() 
            elif op == "12":
                self.chamaFilogenia()
            elif op == "13":
                self.chamaTraducao()
            elif op == "14":
                self.chamaOrfs()
            elif op == "15":
                self.chamaMultiAlign()
            elif op == "16":
                self.sair_do_menu()
                break
            else:
                print("\n Opção não valida. Tente de novo!")
        print("\n Obrigado por ter usado o gestor de sequências!") 
    
    def chamaCriar(self):
        try:
            nome=input("insira o nome da base de dados ")
            gest.criar(nome)
        except:
            print("Erro ao dar o nome à base de dados!")
    
    def chamaFasta(self):
        try:
            filename=input("Escreva a base de dados a importar ")
            gest._import(filename)
        except:
            print("Erro ao importar!")
            
    def chamaAdicionar(self):
        try:
            Id=input("Insira um identificador ")
            seq=input("Introduza uma sequência ")
            gest.adicionar(Id,seq)
        except:
            print("Erro ao adicionar a sequência!")
    
    def chamaAnotacoes(self):
        try:
            Id=input("Insira um identificador da sequência ")
            anotacao=input("Introduza a anotação ")
            tipo=input("Introduza um tipo de anotação ")
            gest.add_anotacoes(tipo,anotacao,Id)
        except:
            print("Erro ao adicionar anotação!")
    
    def chamaEditar(self):
        try:
            Id=input("Insira um identificador da sequência ")
            new=input("Nova sequência:")
            gest.edit_sequencias(Id,new)
        except:
            print("Erro ao editar sequências!")
    
    def chamaValida(self):
        try:
            seq=input("Escreva uma sequência ")
            gest.valida(seq)
        except:
            print("Erro ao validar a sequência!")
    
    def chamaEscritaSeqs(self):
        try:
            filename=input("Ficheiro ")
            ids=input("Identificador ")
            gest.write_seqs(filename,ids)
        except:
            print("Erro ao escrever!")
    
    def chamaProcuraSeqs(self):
        try:
            seq1=input("Escreva a sequência a procurar ")
            gest.procura_seqs(seq1)
        except:
            print("Erro ao procurar!")
    
    def chamaProcuraPad(self):
        try:
            Id=input("Escreva um identificador ")
            pad=input("Introduza um padrão que pretenda procurar ")
            gest.procura_pad(pad,Id)
        except:
            print("Erro ao procurar!")
    
    def chamaFreqSimb(self):
        try:
            s=input("Escreva o símbolo ")
            Id=input("Introduza o identificador ")
            gest.freq_simbolos(s,Id)
        except:
            print("Erro ao procurar!")
    
    def chamaSearchDB(self):
        try:
            pad=input("Introduza um padrão ")
            gest.search_db(pad)
        except:
            print("Erro ao procurar!")
    
    def chamaFilogenia(self):
        try:
            gest.arvore_filogenetica()
        except:
            print("Erro ao criar árvore filogenética!")
    
    def chamaMultiAlign(self):
        try:
            file=input("Nome do ficheiro ")
            gest.multialign(file)
        except:
            print("Erro ao fazer alinhamento múltiplo!")
    
    def chamaTraducao(self):
        try:
            Id=input("Escreva um identificador da sequência ")
            gest.traducao(Id)
        except:
            print("Erro ao efetuar a tradução!")
    
    def chamaOrfs(self):
        try:
            Id=input("Introduza um identificador da sequência ")
            gest.orfs(Id)
        except:
            print("Erro!")
    
    def sair_do_menu(self):
        "Sair do menu"
        print("Saída do menu feita com sucesso!")
        global janela
        if janela is not None:
            del janela
        return True
    
    def do_sair(self,arg):
        "Sair do gestor de sequências"
        print("Obrigado por usar o gestor de sequências!")
        global janela
        if janela is not None:
            del janela
        return True

            
if __name__ == "__main__":
    gest = GestorSeq()
    janela=None
    sh= Shell()
    sh.cmdloop()