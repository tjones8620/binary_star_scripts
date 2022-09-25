#!/usr/bin/python
# -*- coding: iso-8859-15 -*-


import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt("/mnt/local/thomas/WR140_test/trajectory.txt")
t = np.array(data[:,0])
x1 = np.array(data[:,2])
y1 = np.array(data[:,3])
x2 = np.array(data[:,6])
y2 = np.array(data[:,7])

plt.figure()
plt.plot(x1,y1, label="star 1")
plt.plot(x2,y2, label="star 2")
plt.grid()
plt.show()
plt.savefig('/mnt/local/thomas/WR140_test/plot_traj1')

plt.figure()
plt.plot(t,x1, label="star 1, $x$")
plt.plot(t,y1, label="star 1, $y$")
plt.plot(t,x2, label="star 2, $x$")
plt.plot(t,y2, label="star 2, $y$")
plt.grid()
plt.show()
plt.savefig('/mnt/local/thomas/WR140_test/plot_traj2')