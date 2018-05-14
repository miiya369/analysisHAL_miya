#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys, os
import numpy as np

from work_space    import *
from io_mnist_data import *
from activ_funcs   import *

### ======================== Parameters ======================== ###
ifbase  = "/Users/miiya/Dropbox/programs/data/mnist_data"
ifname  = "trained_W.layer0.60000" #None
ofname  = None #"trained_W.layer0.60000"
N_input = 60000
N_test  = 10000
func    = ReLU
dfunc   = d_ReLU

#SizeOfWorkSpace = np.array([28*28, 128, 10])
SizeOfWorkSpace = np.array([28*28,      10])

### =========================== Main =========================== ###
def main ():
    ws = WorkSpace(SizeOfWorkSpace)
    
    if (ifname is None):
        # Larning
        dimg = input_mnist_data_image(ifbase+"/train-images.idx3-ubyte", a_N = N_input)
        dlbl = input_mnist_data_label(ifbase+"/train-labels.idx1-ubyte", a_N = N_input)
        
        for i in range(len(dimg[:,0,0])):
            print("Reading : i =", i)
            print_image(dimg[i,:,:])
            ans = np.zeros(10); ans[dlbl[i]] = 1
            ws.calc(dimg[i,:,:].flatten(), func)
            ws.back_prop(ans, dfunc)
            
        if (ofname is not None):
            ws.save_w(ofname)
    else:
        # Read weights from the file
        ws.load_w(ifname)
    
    # Testing
    dimg_t  = input_mnist_data_image(ifbase+"/t10k-images.idx3-ubyte" , a_N = N_test)
    dlbl_t  = input_mnist_data_label(ifbase+"/t10k-labels.idx1-ubyte" , a_N = N_test)
    
    i_count = np.zeros(10)
    a_count = np.zeros(10)
    for i in range(len(dimg_t[:,0,0])):
        i_count[dlbl_t[i]] += 1
        if (np.argmax(ws.calc(dimg_t[i,:,:].flatten(), func)) == dlbl_t[i]):
            a_count[dlbl_t[i]] += 1
    
    print("")
    for i in range(10):
        if (i_count[i] == 0):
            print("correct rate for %3d = %6d/%6d" % (i, a_count[i], i_count[i]))
        else:
            print("correct rate for %3d = %6d/%6d (= %6.2f %%)" % 
                  (i, a_count[i], i_count[i], float(a_count[i]) / float(i_count[i]) * 100))
    print("")
    print("correct rate for all = %6d/%6d (= %6.2f %%)" %
          (np.sum(a_count), np.sum(i_count), float(np.sum(a_count)) / float(np.sum(i_count)) * 100))
    
    return 0

### ============================================================ ###
def perceptron_1(irate, idlbl, idimg, idlbl_t, idimg_t):
    w = np.zeros(28*28)
    
    for i in range(len(idimg[:,0,0])):
        dtmp = idimg[i,:,:].flatten()
        y    = int(ReLU(np.sum(dtmp * w)))
        if (idlbl[i] != y):
            w += irate * (idlbl[i]-y) * dtmp
    
    i_count = 0
    for i in range(len(idimg_t[:,0,0])):
        y = int(ReLU(np.sum(idimg_t[i,:,:].flatten() * w)))
        if (idlbl_t[i] == y):
            i_count += 1
    
    print("correct rate = %f %%" % (float(i_count) / float(len(idimg_t[:,0,0])) * 100))

def print_img(i, idlbl, idimg):
    print("============================")
    print("index = %d, answer = %d:" % (i, idlbl[i]))
    print_image(idimg[i,:,:])
    print("============================")

### ============================================================ ###
### ============================================================ ###

if (__name__ == "__main__"):
    argv = sys.argv; argc = len(argv)
    
    if (argc != 1):
        sys.exit("usage: python %s" % os.path.basename(argv[0]))
    
    if (main() != 0):
        sys.exit("ERROR EXIT.")
