# -*- coding: utf-8 -*-
# Author: Takaya Miyamoto
# E-mail: takaya.miyamoto@yukawa.kyoto-u.ac.jp
# Since : Sun Feb 18 21:31:38 JST 2018
# 
# Brief : The definition of the class of work space for neural network.

import numpy as np

class WorkSpace():
    """
    The class of work space for neural network.
    
    @ The meaning of member variables
    m_L : #.intermediate layer (>= 0)
    m_N : #.node in each layer (>= 1)
    m_V : The value of output for each node (dtype=float)
    m_W : The weight between each node      (dtype=float)
    
    @ The type and dimension of the member variables
    m_L                                                         (Scalar int)
    m_N[#.in-layer + 2]                                      (1-dim ndarray)
    m_V[#.in-layer + 2][#. node]              (1-dim list and 1-dim ndarray)
    m_W[#.in-layer + 1][#.inode + 1, #.onode] (1-dim list and 2-dim ndarray)
    
    @ Notes
    - Bias is not in each x_0, but each x_N.
    """
    m_L = 0
    m_N =  np.array( [1,1] )
    m_V = [np.array( [0.0] ), np.array( [0.0] )]
    m_W = [np.array([[0.0]]), np.array([[0.0]])]
    
    def __init__(self, a_Nnodes):
        """
        Constructor
        
        @ Arguments
        - a_Nnodes[#.layer] (1-dim ndarray)
        """
        self.reset(a_Nnodes)
    
    def reset(self, a_Nnodes):
        """
        Reset & Initialize the work space
        
        @ Arguments
        - a_Nnodes[#.layer] (1-dim ndarray)
        """
        if (len(a_Nnodes) < 2):
            print("\nERROR: #.layer should be 2 or larger.\n")
            return None
        
        if (np.any(a_Nnodes < 1)):
            print("\nERROR: #.node should be 1 or larger.\n")
            return None
        
        self.m_L = len(a_Nnodes) - 2
        self.m_N = a_Nnodes
        self.m_V = [np.zeros( self.m_N[i] )                  for i in range(self.m_L + 2)]
        self.m_W = [np.zeros((self.m_N[i]+1, self.m_N[i+1])) for i in range(self.m_L + 1)]
    
    def calc(self, a_in, a_func):
        """
        The function to calculate the neural network from input to output
        
        @ Arguments
        - a_in [#.node in first layer] (1-dim ndarray) : input values
        - a_func                            (function) : activation function
        
        @ Return
        - m_V[-1][:]                  (1-dim ndarray) : output values
        """
        if (len(self.m_V[0]) != len(a_in)):
            print("\nERROR: #.node of input data is differ.\n")
            return None
        
        self.m_V[0][:] = a_in[:]
        for n in range(1, self.m_L+2):
            self.m_V[n][:] = np.array(
                [a_func(np.sum(np.append(self.m_V[n-1], np.array([1.0])) *
                               self.m_W[n-1][:,i])) 
                 for i in range(len(self.m_W[n-1][0,:]))])
        
        return self.m_V[-1]
    
    def back_prop(self, a_ans, a_dfunc, a_del = 1e-9):
        """
        The function to calculate the back-propagation by using the answer
        
        @ Arguments
        - a_ans [#.node in final layer] (1-dim ndarray) : answer values
        - a_dfunc                            (function) : deriv. of the activation function
        - a_del                          (scalar float) : step value
        """
        if (len(self.m_V[-1]) != len(a_ans)):
            print("\nERROR: #.node of answer data is differ.\n")
            return None
        
        l_tmpV = np.append(self.m_V[-2], np.array([1.0]))
        l_Ninp = len(self.m_W[-1][:,0])
        l_Nout = len(self.m_W[-1][0,:])
        l_ary1 = (a_ans-self.m_V[-1])
        l_ary2 = np.array([a_dfunc(np.sum(l_tmpV*self.m_W[-1][:,j])) for j in range(l_Nout)])
        l_ary  = l_ary1 * l_ary2
        
        for i in     range(l_Ninp):
            for j in range(l_Nout):
                self.m_W[-1][i,j] += a_del * l_tmpV[i] * l_ary[j]
        
        for n in range(-2, -self.m_L-2, -1):
            l_tmpV = np.append(self.m_V[n-1], np.array([1.0]))
            l_Ninp = len(self.m_W[n][:,0])
            l_Nout = len(self.m_W[n][0,:])
            l_ary1 = np.array([        np.sum(l_ary1 *
                                              self.m_W[n+1][j,:])  for j in range(l_Nout)])
            l_ary2 = np.array([a_dfunc(np.sum(l_tmpV *
                                              self.m_W[n  ][:,j])) for j in range(l_Nout)])
            l_ary  = l_ary1 * l_ary2
            
            for i in     range(l_Ninp):
                for j in range(l_Nout):
                    self.m_W[n][i,j] -= a_del * l_tmpV[i] * l_ary[j]
    
    def save_w(self, a_ofname, verbose_flg = True):
        """
        The function to output the current weights
        
        @ File format
        - (1) Magic number         (32bit-integer)           : 19900518
        - (2) #.intermediate layer (32bit-integer)           : m_L
        - (3) m_N[#.in-layer + 2]                                      (1-dim ndarray)
        - (4) m_W[#.in-layer + 1][#.inode + 1, #.onode] (1-dim list and 2-dim ndarray)
        
        Note: Save as Little endian
        """
        from struct import pack
        
        MagicNum_BinData_Miya = 19900518
        
        with open(a_ofname, 'wb') as ofile:
            ofile.write(pack('<i', MagicNum_BinData_Miya))
            ofile.write(pack('<i', self.m_L))
            ofile.write(pack('<%di' % (self.m_L+2), *self.m_N))
            
            for i in range(self.m_L+1):
                for j in range(self.m_N[i]+1):
                    ofile.write(pack('<%dd' % self.m_N[i+1], *self.m_W[i][j,:]))
        
        if (verbose_flg):
            print("# Successful to output bin data.")
    
    def load_w(self, a_ifname, verbose_flg = True):
        """
        The function to input the current weights
        
        @ File format
        - (1) Magic number         (32bit-integer)           : 19900518
        - (2) #.intermediate layer (32bit-integer)           : m_L
        - (3) m_N[#.in-layer + 2]                                      (1-dim ndarray)
        - (4) m_W[#.in-layer + 1][#.inode + 1, #.onode] (1-dim list and 2-dim ndarray)
        
        Note: Save as Little endian
        """
        from struct import unpack
        
        MagicNum_BinData_Miya = 19900518
        
        with open(a_ifname, 'rb') as ifile:
            l_MagicNum = unpack('<i', ifile.read(4))[0]
            
            if (l_MagicNum != MagicNum_BinData_Miya):
                print("\nERROR: Invalid magic number '%d'.\n" % l_MagicNum)
                return None
            
            self.m_L = unpack('<i', ifile.read(4))[0]
            self.m_N = np.array([unpack('<i', ifile.read(4))[0] for i in range(self.m_L+2)])
            self.m_W = [np.array([[unpack('<d', ifile.read(8))[0] 
                                   for k in range(self.m_N[i+1])]
                                  for j in range(self.m_N[i]+1)])
                        for i in range(self.m_L+1)]
        
        if (verbose_flg):
            print("# Successful to input bin data.")
            print("# N.intermediate  layer =", self.m_L)
            print("# N.nodes in each layer =", self.m_N)
