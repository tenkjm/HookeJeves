import numpy as np
import time
from sys import argv, stdout, exit
from skimage.io import imread
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw


def load_data(path):
    fi = open(path)
    lines = [line if line[-1] != '\n' else line[:-1]
             for line in fi.readlines()]
    fi.close()
    i = 0
    data = []
    gt = []

    while i < len(lines):
        data.append([float(n) for n in lines[i][1: len(lines[i])-1].split(",")])
        i= i+1

    a=data[0][0]

    norma = []
    x =[]
    i = 0
    while i < len(data):
        dist = np.linalg.norm(np.array( data[i]) - np.array( data[len(data)-1]))
        norma.append(dist)
        x.append(i)
        i = i + 1


    print(norma[1])

    #plt.plot(x, norma, 'g--^', label='regional')

    #plt.show()
    return (x, norma)


filename = argv[1]
data1=load_data(filename)

filename = argv[2]
data2=load_data(filename)

filename = argv[3]
data3=load_data(filename)

filename = argv[4]
data4=load_data(filename)



fig = plt.figure(figsize=(6, 4))
#t = np.arange(-5.0, 1.0, 0.1)
sub1 = fig.add_subplot(221) # instead of plt.subplot(2, 2, 1)
sub1.set_title('HJ std') # non OOP: plt.title('The function f')
sub1.plot(data1[0], data1[1])
sub2 = fig.add_subplot(222, axisbg="lightgrey")
sub2.set_title('HJ rnd explores')
sub2.plot(data2[0], data2[1])
#t = np.arange(-3.0, 2.0, 0.02)
sub3 = fig.add_subplot(223)
sub3.set_title('HJ linear')
sub3.plot(data3[0], data3[1])
#t = np.arange(-0.2, 0.2, 0.001)
sub4 = fig.add_subplot(224, axisbg="lightgrey")
sub4.set_title('HJ var explores')
sub4.plot(data4[0], data4[1])
#plt.plot(t, g(t))
plt.tight_layout()
plt.show()
