from sympy import *

nodes_1 = [0,1,1.2,2]
nodes_2 = [0,1,1.2,2]

elements_1 = [[0,1], [1,2], [2,3]]
elements_2 = [[0,1], [1,2], [2,3]]

# Create the inversed ordering
nodes_2.reverse()

print "Left-right ordering"
print nodes_1
print elements_1
print "Right-left ordering"
print nodes_2
print elements_2

# Insert node at the correct position. Element elements_1[2] is still correct since it points to the index only
print "Fixing left-right ordered elements:"
nodes_1.insert(3,1.6)
elements_1.append([3,4])
print nodes_1
print elements_1

# In this part, I assume that also the local indices should be ordered right to left
print "Fixing right-left ordered elements:"
nodes_2.insert(1,1.6)
elements_2.append([3,4])
print nodes_2
print elements_2