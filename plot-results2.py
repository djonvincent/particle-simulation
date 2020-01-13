import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('results2.csv', delimiter=',')[1:]
plt1 = plt.figure(1)
plt.plot(data[:,0], data[:,1], label='Serial')
plt.plot(data[:,0], data[:,2], label='Parallel')
plt.legend(loc=2)
plt.title('Time to run 10 timesteps')
plt.xlabel('Number of bodies')
plt.ylabel('Time (s)')
plt1.show()
plt2 = plt.figure(2)
plt.plot(data[:,0], data[:,1] / data[:,2])
plt.ylim(2.74,4.75)
plt.title('Time to run 10 timesteps')
plt.xlabel('Number of bodies')
plt.ylabel('Speedup')
plt2.show()
input()
