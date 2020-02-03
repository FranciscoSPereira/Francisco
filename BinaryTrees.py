class BinaryTree:
	'''
	Estruturas de dados recursivas
	'''

	def __init__(self, val, dist=0.0, left=None, right=None):
		'''
		Construtor que recebe o valor, a distância e as duas sub-árvores
		None (a árvore esta vazia)
		'''
		self.value = val
		self.distance = dist
		self.left = left
		self.right = right

	def getCluster(self):
		'''
		Chamar o método para a sub-árvore esquerda, se esta não for nula
		Chamar o método para a sub-árvore direita, se esta não for nula
		Juntar os dois conjuntos
		'''
		res = []
		if self.value >= 0:
			res.append(self.value)		# se for folha, acrescenta a lista res
		else:
			if self.left is not None:
				res.extend(self.left.getCluster())		# .extend - adiciona à lista res
			if self.right is not None:
				res.extend(self.right.getCluster())		# o extend faz com que não adicione uma lista res nova à lista res original, assim fica tudo na mesma lista
		return res
    
	def printtree(self):
		self.printtree_rec(0, 'Root')

	def printtree_rec(self, level, side):
		tabs = ''				# lista que guarda tabs
		for i in range(level):
			tabs += '\t'		# imprime um tab na lista tabs, se level for 0 não imprime
			if i is level-1:
				tabs += '|_'
		if self.value >= 0:			# >= 0 ou > 0 ? da dist ou value?
			print(tabs, side, '-> value:', self.value)		# se for uma folha, imprime a lista de tabs, o side e o value da folha
		else:
			print(tabs, side, '-> Dist.:', self.distance)	 # se for um nó, faz o mesmo com o valor da distância desse nó
		if self.left is not None:
			self.left.printtree_rec(level + 1, 'Left')		# se existir sub-árvore esquerda volta  fazer de inicio adicionando 1 tab
		if self.right is not None:
			self.right.printtree_rec(level + 1, 'Right')	 # se existir sub-árvore direita volta  fazer de inicio adicionando 1 tab


# --------------------------\\-----------------------------

# TESTS

def test():
	a = BinaryTree(4)
	b = BinaryTree(3)
	c = BinaryTree(2)
	d = BinaryTree(1)
	e = BinaryTree(-1, 1.5, a, d)
	f = BinaryTree(-1, 2, c, b)
	g = BinaryTree(-1, 4.5, f, e)
	g.printtree()


if __name__ == '__main__':
	test()