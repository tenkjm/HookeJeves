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

    plt.plot(x, norma, 'g--^', label='regional')

    plt.show()



filename = argv[1]
load_data(filename)

