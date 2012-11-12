from sympy import *

nodes_1 = [0,1,1.2,2]
nodes_2 = [0,1,1.2,2]

elements_1 = [[0,1], [1,2], [2,3]]
elements_2 = [[0,1], [1,2], [2,3]]

# Create the inversed ordering
nodes_2.reverse()
[el.reverse() for el in elements_2]
elements_2.reverse()

print "Left-right ordering"
print nodes_1
print elements_1
print "Right-left ordering"
print nodes_2
print elements_2

# Insert node at the correct position. Element elements_1[2] is still correct since it points to the index only
nodes_1.insert(3,1.6)
elements_1.append([3,4])
print "Fixing elements:"
print nodes_1
print elements_1
