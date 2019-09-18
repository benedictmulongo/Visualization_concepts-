import numpy as np


def barycentric_triangle_function(x,y):
    
    solution = (1 - y - x) * 1 + (x) * 2 + y * 3
    print("Barycentric(", x, ", ", y ,") = ", solution )
    return  solution
    
def shepard_triangle_function(x,y):
    
    w00 = 1 / np.sqrt(x**2 + y**2)
    w10 = 1 / np.sqrt((x-1)**2 + y**2)
    w01 = 1 / np.sqrt(x**2 + (y-1)**2)
    
    N = w00 + w10 + w01
    
    solution = (w00 / N) * (1)  + (w10 / N) * (2)  + (w01 / N) * (3) 
    
    print("Shepard(", x, ", ", y ,") = ", solution )
    return  solution
    
print()
barycentric_triangle_function(1/2, 0)
barycentric_triangle_function(0, 1/2)
barycentric_triangle_function(1/2, 1/2)
barycentric_triangle_function(2/3, 2/3)
print()
shepard_triangle_function(1/2, 0)
shepard_triangle_function(0, 1/2)
shepard_triangle_function(1/2, 1/2)
shepard_triangle_function(2/3, 2/3)

    
