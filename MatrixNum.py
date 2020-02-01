import math
import numpy as np


class MatrixNum:

    def __init__(self, rows, cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for j in range(cols):
                self.mat[i].append(0.0)         # pq 0.0? pcausa dos scores?

    def __getitem__(self, ij):      # replica o modo de utilizaçao das listas com indices?
        i, j = ij
        return self.mat[i][j]
    
    def __setitem__(self, ij, value):
        i, j = ij
        self.mat[i][j] = value

    def numRows(self):
        return len(self.mat)
    
    def numCols(self):
        return len(self.mat[0])
    
    def getValue(self, i, j):           # equivalente a __getitem__()
        return self.mat[i][j]           # opcional?
    
    def setValue(self, i, j, value):
        self.mat[i][j] = value
    
    def printmat(self):
        for r in self.mat:
            print(r)
        print()

    def sumMat(self):
        s = 0.0
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                s += self.mat[i][j]
        return s

    def mean(self):
        return self.sumMat()/(self.numRows() * self.numCols())

    def maximum(self):
        m = self.mat[0][0]
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if self.mat[i][j] > m:
                    m = self.mat[i][j]
        return m

    def minimum(self):
        m = self.mat[0][0]
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
        return m

    def maxIndexes(self):      # return os indices do ponto (tupulo) que tem o maior valor
        m = self.mat[0][0]          # valor do ponto
        res = (0, 0)                # indices desse ponto
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if self.mat[i][j] > m:
                    m = self.mat[i][j]
                    res = (i, j)
        return res

    def minIndexes(self):
        m = self.mat[0][0]
        res = (0, 0)
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
                    res = (i, j)
        return res

    def minDistIndexes (self):
        m = self.mat[1][0]
        res = (1, 0)
        for i in range(self.numCols()):
            for j in range(i+1, self.numRows()):
                if self.mat[j][i] < m:
                    m = self.mat[j][i]
                    res = (j, i)
        return res

    def sumRow(self, li):
        return sum(self.mat[li])

    def meanRow(self, li):      # faz a media de todos os elementos de uma linha, que e dada
        s = 0
        for k in range(len(self.mat[li])):
            s += self.mat[li][k]
        return s/len(self.mat[li])

    def maxRow(self, li):
        m = self.mat[li][0]
        for k in range(1, len(self.mat[li])):
            if self.mat[li][k] > m:
                m = self.mat[li][k]
        return m

    def minRow(self, li):
        m = self.mat[li][0]
        for k in range(1, len(self.mat[li])):
            if self.mat[li][k] < m:
                m = self.mat[li][k]
        return m

    def meanRows(self):         # faz media de todos os elementos de cada linha (media por linha)
        res = []
        for r in range(self.numRows()):
            res.append(self.meanRow(r))
        return res

    def minRows(self):
        res = []
        for r in range(self.numRows()):
            res.append(self.minRow(r))
        return res

    def maxRows(self):
        res = []
        for r in range(self.numRows()):
            res.append(self.maxRow(r))
        return res

    def sumCol(self, lc):
        s = 0
        for k in range(self.numRows()):
            s += self.mat[k][lc]
        return s

    def meanCol(self, lc):
        s = 0
        for k in range(self.numRows()):
            s += self.mat[k][lc]
        return s / self.numRows()

    def maxCol(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.numRows()):
            if self.mat[c][lc] > m:
                m = self.mat[c][lc]
        return m

    def minCol(self, lc):
        m = self.mat[0][lc]
        for c in range(1, self.numRows()):
            if self.mat[c][lc] < m:
                m = self.mat[c][lc]
        return m

    def sumCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.sumCol(c))
        return res

    def meanCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.meanCol(c))
        return res

    def maxCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.maxCol(c))
        return res

    def minCols(self):
        res = []
        for c in range(self.numCols()):
            res.append(self.minCol(c))
        return res

    def addRow(self, newrow):
        self.mat.append(newrow)

    def addCol(self, newcol):
        for r in range(self.numRows()):
            self.mat[r].append(newcol[r])

    def removeRow(self, ind):
        del self.mat[ind]

    def removeCol(self, ind):
        for r in range(self.numRows()):
            del self.mat[r][ind]

    def copy(self):
        newm = MatrixNum(self.numRows(), self.numCols())
        for i in range(self.numRows()):
            for j in range(self.numCols()):
                newm.mat[i][j] = self.mat[i][j]
        return newm

    def csv_data(self, file):           # le a lista de listas com os valores do csv e passa-os para o self.mat
        i = 0
        for row in file:
            j = 0
            for col in file:
                self.setValue(i, j, file[i][j])
                j += 1
            i += 1


def readcsv(path_filename):
    file_array = np.genfromtxt(path_filename, delimiter=',')
    file = file_array.tolist()
    m = MatrixNum(len(file), len(file[0]))
    m.csv_data(file)
    return m


# fazer uma deepcopy da matriz onde se pode aplicar funçoes em todos os elementos (ex: tipo func log)?
# somar e multiplicar matrizes (numpy?)


# --------------------------\\-----------------------------

# TESTS

def test():
    m = MatrixNum(3, 2)
    m.setValue(2, 1, 3.0)
    m.setValue(1, 0, 2.0)
    m.setValue(0, 1, 5.0)
    m.setValue(1, 1, -2.0)
    m.printmat()
    print(m.mean())
    print(m.maximum())
    print(m.maxIndexes())
    print('Mean rows')
    print(m.meanRows())
    m.printmat()


def test2():
    m = MatrixNum(3, 2)
    m.setValue(2, 1, 3.0)
    m.setValue(1, 0, 2.0)
    m.setValue(0, 1, 5.0)
    m.setValue(1, 1, -2.0)
    m.printmat()
    print(m.mean())
    print(m.maximum())
    print(m.minIndexes())
    print("Mean rows")
    print(m.meanRows())
    print("Mean cols")
    print(m.meanCols())
    print("")
    m.removeRow(1)
    m.printmat()
    m.addRow([2.0, 4.0])
    m.printmat()
    nm = m.copy()
    nm.printmat()
    nm.addCol([3.0, 4.0, 5.0])
    nm.printmat()
    nm.removeCol(0)
    nm.printmat()
    print(nm.meanRows())
    print()


def test3():
    m = readcsv('./teste.csv')
    m.printmat()


if __name__ == '__main__':  # testa se esta script e lançada
    test()
    test2()
    test3()