#!/usr/bin/python
#
# Author: Takaya Miyamoto
# Date  : Thu Sep  8 23:56:36 JST 2016
# Brief : Take difference for the header of two NBSdata
#

if __name__ == "__main__":
   import sys
   import math
   import os.path
   import struct
   
   argv = sys.argv
   argc = len(argv)
   
   if (argc != 3):
      print("usage: %s [ifile1] [ifile2]" % argv[0])
      quit(1)
   
   ifile1 = open(argv[1].strip(), 'rb')
   ifile2 = open(argv[2].strip(), 'rb')

   for i in range(5):
      d1 = struct.unpack('<i', ifile1.read(4))[0]
      d2 = struct.unpack('<i', ifile2.read(4))[0]
      if (d1 != d2):
         print("Two file are differ.")
         quit(1)
   
   length = d1
   
   for i in range(length + 11):
      d1 = struct.unpack('<i', ifile1.read(4))[0]
      d2 = struct.unpack('<i', ifile2.read(4))[0]
      if (d1 != d2):
         print("Two file are differ.")
         quit(1)
   
   ifile1.close()
   ifile2.close()
