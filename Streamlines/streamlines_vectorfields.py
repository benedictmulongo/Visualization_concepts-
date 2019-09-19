
import numpy as np
import matplotlib.pyplot as plt

def plot_arrow_at_xy(x,y) :
    
    length = 0.5
    width=0.5
    fc="r"
    ec="k"
    angle = np.deg2rad(45.0)
    plt.axis("equal")
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.arrow(x,y, length * np.cos(angle), length * np.sin(angle),
                    fc=fc, ec=ec, head_width=width, head_length=width)
    plt.show()

def angle_xy(x,y) :
    
    return np.arctan2(y,x)
    
def plot_arrow(x, y, length):  

    #length = 0.5
    width=0.5
    fc="r"
    ec="k"
    angle = angle_xy(x,y)
    #angle = np.deg2rad(45.0)
    plt.axis("equal")
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    for i in range(len(x)):
        
        plt.arrow(x[i],y[i], length[i] * np.cos(-angle[i]), length[i] * np.sin(angle[i]),
                    fc=fc, ec=ec, head_width=width, head_length=width)
                    
    x = np.linspace(-5,5,11)
    y = np.linspace(-5,5,11)
    X, Y = np.meshgrid(x, y)
    plt.plot(X, Y, 'k')
    plt.plot(Y, X, 'k')
    plt.show()

def plot_arrow2(plt, x, y, length):  

    width=0.5
    fc="r"
    ec="k"
    angle = angle_xy(x,y)
    #angle = np.deg2rad(45.0)
    plt.axis("equal")
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    for i in range(len(x)):
        
        plt.arrow(x[i],y[i], length[i] * np.cos(angle[i]), length[i] * np.sin(angle[i]),
                    fc=fc, ec=ec, head_width=width, head_length=width)
    

## Case 1 The velocity field is in its analytical forms and the analytical integration is possible
## V(x,y) = (x,y)  and the streamline is y = Cx

def fiedlines_ex(C = 0, N = 10 ):
    
    # N = Number of points
    x = np.linspace(-5,5,N)
    
    #  The streamline is y = Cx
    y = C*x
    
    length = np.sqrt(x**2 + y**2)
    
    # Scale to [0,1]
    length_norm = ( length -  length.min() ) / (length.max() - length.min())
    
    return x,y, length_norm
    
def testCase1() :
        
    # x,y, length = fiedlines_ex(C = 1, N = 10 ) 
    # print(fiedlines_ex(C = 1, N = 10 ))
    # print()
    # plot_arrow(x, y, length)
    
    # Plot Streamlines
    cc = np.linspace(-10,10,30)
    for i in cc.tolist() :
        x,y, length = fiedlines_ex(C = i, N = 10 ) 
        plot_arrow2(plt, x, y, length)
        plt.show()


def euler(fcn, a, b, y0, N):
    """
    Function for generating an Euler solution to
    y ’ = f(x,y) in ‘N‘ steps with initial
    condition y[a] = y0 .
    """
    h = (b - a) / N
    x = a + np.arange(N + 1) * h 
    y = np.zeros(x.size)
    y[0] = y0 
    
    for k in range(N) :
        
        y[k+1] = y[k] + h * fcn(x[k], y[k])
        
    return (x,y)
    
def runge_kutta_methods():
    
    return 0
    
def func(s0):
    
    x = s0[0]
    y = s0[1]
    
    back = np.array([1, (6*(x**2) - 1)*y])
    return back
    
def func2(s0):
    
    x = s0[0]
    y = s0[1]
    
    back = np.array([-y, x/2])
    return back

def euler_vectorfields(func, x0,y0, N):
    
    """
    N : Number of steps 
    dt : The step size 
    
    """
    
    values = []
    
    dt = (y0 - x0) / N
    s0 = np.array([x0,y0])
    #v0 = n
    # Test specific step size comment afterwards
    dt = 1/2
    
    values.append(s0.tolist())
    
    for i in range(N):
        s0 = s0 + func(s0)*dt
        values.append(s0.tolist())

    return np.array(values )
    
def RK4_vectorfields(func, x0,y0, N):
    
    """
    N : Number of steps 
    dt : The step size 
    
    """
    
    values = []
    
    dt = (y0 - x0) / N
    s0 = np.array([x0,y0])
    #v0 = n
    # Test specific step size comment afterwards
    dt = 1/2
    
    values.append(s0.tolist())
    
    for i in range(N):
        print("index : ", i)
        k1 = dt*func(s0)
        k2 = dt*func(s0 + k1/2) 
        k3 = dt*func(s0 + k2/2 )
        k4 = dt*func(s0 + k3)
        s0 = s0 + (k1 + 2*(k2 + k3) + k4)/6
        values.append(s0.tolist())

    return np.array(values )
    
    
def test_eulers() :
    
    N = 30
    x0 = 0
    y0 = -1
    print(euler_vectorfields(func2, x0,y0, N))
    matrix = euler_vectorfields(func2, x0,y0, N)
    x = np.array(matrix[:,0])
    y = np.array(matrix[:,1])
    length = np.sqrt(x**2 + y**2)
    # Scale to [0,1]
    length_norm = ( length -  length.min() ) / (length.max() - length.min())
    
    print("x = ", np.array(x))
    print("y : ", y)
    
    plot_arrow(x, y, length_norm)
    
def test_RK4() :
    
    
    plt.figure()
    
    N = 30
    x0 = 0
    y0 = -1

    #matrix = euler_vectorfields(func2, x0,y0, N)
    matrix = RK4_vectorfields(func2, x0,y0, N)
    x = np.array(matrix[:,0])
    y = np.array(matrix[:,1])
    length = np.sqrt(x**2 + y**2)
    # Scale to [0,1]
    length_norm = ( length -  length.min() ) / (length.max() - length.min())
    

    plot_arrow(x, y, length_norm)

    
    
test_eulers() 
    
test_RK4()
    
    