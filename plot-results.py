import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('results.csv', delimiter=',')[1:]
x = np.log10(data[1:,0])
#y1 = np.log10(np.absolute(data[:-1,1] - data[-1,1]))
#y2 = np.log10(np.absolute(data[:-1,2] - data[-1,2]))
#y3 = np.log10(np.absolute(data[:-1,3] - data[-1,3]))
#y4 = np.log10(np.absolute(data[:-1,4] - data[-1,4]))
y1 = np.log10(np.absolute(data[1:,1] - data[:-1,1]))
y2 = np.log10(np.absolute(data[1:,2] - data[:-1,2]))
y3 = np.log10(np.absolute(data[1:,3] - data[:-1,3]))
y4 = np.log10(np.absolute(data[1:,4] - data[:-1,4]))
plt1 = plt.figure(1)
plt.scatter(x, y1, label='x-coordinate', marker='x', c='orange')
plt.scatter(x, y2, label='y-coordinate', marker='x', c='blue')
best_fit = np.poly1d(np.polyfit(
    np.concatenate((x,x)), np.concatenate((y1,y2)), 1
))
print(best_fit)
plt.plot(x, best_fit(x)) 
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel('log10(h)')
plt.ylabel('log10(change)')
plt.title('Change in position of first collision point')
plt.legend(loc=2)
plt1.show()
plt2 = plt.figure(2)
plt.scatter(x, y3, label='x-coordinate', marker='x', c='green')
plt.scatter(x, y4, label='y-coordinate', marker='x', c='red')
best_fit = np.poly1d(np.polyfit(
    np.concatenate((x,x)), np.concatenate((y3,y4)), 1
))
print(best_fit)
plt.plot(x, best_fit(x)) 
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xlabel('log10(h)')
plt.ylabel('log10(change)')
plt.title('Change in position of second collision point')
plt.legend(loc=2)
plt2.show()
input()
