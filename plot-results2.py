import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
data = np.genfromtxt('results2.csv', delimiter=',')[1:]
plt1 = plt.figure(1)
plt.plot(data[:,0], data[:,1], label='1 core', c='red')
plt.plot(data[:,0], data[:,2], label='2 cores', c='blue')
plt.plot(data[:,0], data[:,3], label='3 cores', c='orange')
plt.plot(data[:,0], data[:,4], label='4 cores', c='green')
plt.legend(loc=2)
plt.title('Time to run 10 timesteps')
plt.xlabel('Number of bodies')
plt.ylabel('Time (s)')
plt1.show()
plt2 = plt.figure(2)
plt.plot(data[:,0], data[:,1] / data[:,2], label='2 cores', c='blue')
plt.plot(data[:,0], data[:,1] / data[:,3], label='3 cores', c='orange')
plt.plot(data[:,0], data[:,1] / data[:,4], label='4 cores', c='green')
plt.legend(loc=2)
plt.title('Speedup')
plt.xlabel('Number of bodies')
plt.ylabel('Speedup')
plt2.show()
plt3 = plt.figure(3)
plt.plot([2,3,4], data[-1,1]/data[-1,2:])
plt.plot([2,3,4], data[-2,1]/data[-2,2:])
plt.plot([2,3,4], data[-3,1]/data[-3,2:])
plt.plot([2,3,4], data[-4,1]/data[-4,2:])
plt.plot([2,3,4], data[-5,1]/data[-5,2:])
plt.plot([2,3,4], data[-6,1]/data[-6,2:])
plt.plot([2,3,4], data[-7,1]/data[-7,2:])
plt.plot([2,3,4], data[-8,1]/data[-8,2:])
plt.plot([2,3,4], data[-9,1]/data[-9,2:])
plt.plot([2,3,4], data[-10,1]/data[-10,2:])
plt.xticks([2,3,4], [2,3,4])
plt.xlabel('Number of cores')
plt.ylabel('Speedup')
plt.title('Speedup')
plt3.show()
def strong_scaling(p, f):
    return 1 / (f + (1-f)/p)
strong_f, _ = curve_fit(
    strong_scaling,
    [2,3,4],
    data[-1,1]/data[-1, 2:]
)
print(f'Strong f: {strong_f}')
def weak_scaling(p, f):
    return f + (1-f)*p
weak_f, _ = curve_fit(
    weak_scaling,
    [2,3,4],
    data[-1,1]/data[-1, 2:]
)
print(f'Weak f: {weak_f}')
input()
