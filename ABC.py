# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 16:29:07 2020

@author: Francisco
"""

from abc import ABC

class Seq(ABC):
    
    def __init__(self,seq):
        '''
        Construtor de uma sequência
        '''
        self.seq=seq.upper()
        self.alfabeto='ACTG'
        self.propriedades={}
         
    def __eq__(self, obj):
        if isinstance(obj,Seq):
            return self.seq == obj.seq
        else:
            return False
        
    def __str__(self):
        return self.seq
    
    def __getitem__(self,n):
        return self.seq[n]
    
    def __getslice__(self,i,j):
        return self.seq[i,j]
    
    def __len__(self):
        return len(self.seq)

    def add_prop(self,key,value):
        self.propriedades[key]=value
        
    def get_prop(self,key):
        if key not in self.propriedades.keys():
            raise IndexError()
        return self.propriedades[key]

    def getlength(self):
        '''
        Retorna o tamanho da sequência de DNA.
        '''
        return len(self.seq)
    
    def traduzCodao(self, cod):
        '''
        O dicionário tc mostra a correspondência entre o tripleto (a chave) e o respetivo aminoácido.
        Se o codão (cod) escrito for igual a um que se encontra no dicionário, o mesmo é adicionado à string aa, que irá representar uma proteína.
        Caso contrário, é adicionado um X para indicar ausência de um aminoácido.
        '''
        tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCC": "A", "TGT": "C", "TGC": "C",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
              "ATA": "I", "ATT": "I", "ATC": "I",
              "AAA": "K", "AAG": "K",
              "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "ATG": "M", "AAT": "N", "AAC": "N",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "CAA": "Q", "CAG": "Q",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
              "TGG": "W",
              "TAT": "Y", "TAC": "Y",
              "TAA": "_", "TAG": "_", "TGA": "_"}
        if cod in tc:
            aa = tc[cod]
        else:
            aa = "X"  
        return aa
    
    def traducao(self, iniPos=0):
        '''
        Retorna a sequência de aminoácidos codificada pela sequência de DNA
        '''
        seqM = self.sequencia
        seqAA = ""
        for pos in range(iniPos, len(seqM)-2, 3):
            cod = seqM[pos:pos+3]
            aa = self.traduzCodao(cod)
            seqAA = seqAA + aa
        return seqAA
    
    
    def valida(self):
        '''
        Valida a sequência criada como uma sequêncía de DNA
        '''
        for i in range(len(self.sequencia)):
            if self.seq[i] not in self.alfabeto:
                return ('Sequencia de DNA inválida')
        return ('Sequencia de DNA válida')
    
    def printsequencia(self):
        '''
        Imprime a sequência no terminal

        '''
        print(self.seq)
        
        
if __name__ == '__main__':
    sequence = Seq('ACGGGTCACACGT')
    print(sequence.valida()) 
    print(sequence.traducao())