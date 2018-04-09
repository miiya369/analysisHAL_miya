# -*- coding: utf-8 -*-
# Author: Takaya Miyamoto
# E-mail: takaya.miyamoto@yukawa.kyoto-u.ac.jp
# Since : Sun Feb 18 21:31:38 JST 2018
# 
# Brief : The definitions of the activation function

def step(a_x):
    return (1 if a_x >= 0.0 else 0)

def sigmoid(a_x):
    return 1.0 / (1.0 + np.exp(-a_x))

def ReLU(a_x):
    return max(0, a_x)

def d_sigmoid(a_x):
    return sigmoid(a_x) * (1.0 - sigmoid(a_x))

def d_ReLU(a_x):
    return step(a_x)

