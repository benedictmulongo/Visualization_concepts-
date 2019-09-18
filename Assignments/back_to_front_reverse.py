
C = [0, 2, 4, 2, 8, 7, 9, 5, 1]
alpha = [1/6, 1/3, 1/6, 1/6, 1/2, 1/6, 1/3, 1/6, 1/6]

def backToFront(C, alpha) :
    
    I = C[0]
    opacity = 0
    
    for i in range(1,len(C)) :
        I = C[i]*alpha[i] + (1 - alpha[i])*I
        
        opacity = (1 - alpha[i])*opacity + alpha[i]
        
    return I, opacity
    
def FrontToBack(C, alpha) :
    
    I = 0 
    T = 1 
    
    index = list(reversed(list( range(1,len(C)) )))
    for i in index :
        # print(i)
        I = I +  C[i]*alpha[i]*T 
        # print("I : ", I)
        T = (1 - alpha[i]) * T 
        print("Index : ", i ,  "Opacity : ",  1 - T)
    print("I : ", I)
    print("Opacity : ", 1 -T )
    I = I + T*C[0]
    
    return I 
    
def accumulate_opacity(C, alpha) :
    
    
    Acc = 1
    for i in range(1,len(C)) :
        Acc = Acc * (1 - alpha[i])
        
    opacity = 1 - Acc
    
    return opacity
    
    
def accumulate_intensity(C, alpha) :
    
    
    Acc = 1
    Sum = 0
    for i in range(1,len(C)) :
        Acc = 1
        for j in range(1,i-1) :
            Acc = Acc * (1 - alpha[j])
            
        Sum = Sum + C[i]*Acc
        
    
    return Sum
    
    
    
print(len(C))
print(len(alpha))
print("backToFront(C, alpha) : ", backToFront(C, alpha))
print("FrontToBack(C, alpha) : ",FrontToBack(C, alpha))


print()
print("Acc. opacity : ", accumulate_opacity(C, alpha))
print()
print("Acc. intensity : ", accumulate_intensity(C, alpha))
    


