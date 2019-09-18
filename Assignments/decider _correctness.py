
import numpy as np


def midpoint_decider(cell, values, iso):
    
    alpha = (values[0,0] + values[1,0] + values[1,1] + values[0,1]) / 4
    
    if iso < alpha :
        print("K-1")
    else :
        print("K-2")
    
    return alpha
    

def symptotic_decider(cell, values, iso):
    
    value = (values[0,0]*values[1,1] - values[1,0]*values[0,1])
    alpha = value / (values[0,0]  - values[1,0]  - values[0,1] + values[1,1]) 
    
    if iso < alpha :
        print("K-1")
    else :
        print("K-2")

    return alpha
    
    
def check(values) :
    
    A = values[0,0]
    B = values[1,0]
    C = values[1,1]
    D = values[0,1]
    
    return (A**2 - 2*A*C - 6*B*D - B**2 + C**2 - D**2)
    
    
    
# cell_1_values = np.array([[50,-2],[-300,40]])
# cell_1_sign = np.array([[1,-1],[-1,1]])
# print()
# print("Case 1")
# print("+1 ______ -1 ")
# print("  |      |   ")
# print("  |      |   ")
# print("  |______|   ")
# print("-1        +1 ")
# iso = 30
# print("symptotic_decider : ", symptotic_decider(cell_1_sign, cell_1_values, iso))
# print()
# print("midpoint_decider : ", midpoint_decider(cell_1_sign, cell_1_values, iso))

cell_2_values = np.array([[-50,2],[310,-40]]) # ok
# cell_2_values = np.array([[-50,10],[310,-70]])
# cell_2_values = np.array([[5,100],[310000,7]])
cell_2_sign = np.array([[-1,1],[1,-1]])
print()
print("Case 2")
print("-1 ______ +1 ")
print("  |      |   ")
print("  |      |   ")
print("  |______|   ")
print("+1        -1 ")

iso = 0
print("symptotic_decider : ", symptotic_decider(cell_2_sign, cell_2_values, iso))
print()
print("midpoint_decider : ", midpoint_decider(cell_2_sign, cell_2_values, iso))
print()
print(cell_2_values)
print()
print("check(values) : ", check(cell_2_values))
