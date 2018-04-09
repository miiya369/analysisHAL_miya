#!/usr/bin/python
#
# Author: Takaya Miyamoto
# Date  : Thu Sep  8 23:56:36 JST 2016
# Brief : Take difference for two data and show results
#

fname1 = "UNKNOWN"
fname2 = "UNKNOWN"

data_type = -1
DataSize  = -1

ratio_flg = 0
max_flg   = 0
spin_flg  = 2
criterion = 999

ArraySize = 512

Count_Zero_r = [0] * 1
Count_p_r    = [0] * ArraySize
Count_m_r    = [0] * ArraySize
Count_Zero_i = [0] * 1
Count_p_i    = [0] * ArraySize
Count_m_i    = [0] * ArraySize

Count_Zero_r_s1 = [0] * 1 # For BBwave (Total spin=1)
Count_p_r_s1    = [0] * ArraySize
Count_m_r_s1    = [0] * ArraySize
Count_Zero_i_s1 = [0] * 1
Count_p_i_s1    = [0] * ArraySize
Count_m_i_s1    = [0] * ArraySize

### Print Usage ###
def print_usage(Fname):
   print("usage: %s [ifile1] [ifile2] [0/1/2] {option}" % Fname)
   print
   print("args3: 0 = Ishii-san's 1/48 conpressed NBS data")
   print("     : 1 = Not compressed NBS data")
   print("     : 2 = text 2pt-correlator data")
   print
   print("option: -ratio, -max_diff, -spin0, -spin1, -diff_rate [rate]");
   quit(1)

### Arguments setting ###
def set_args(Argc, Argv):
   global fname1, fname2, data_type, ratio_flg, max_flg, spin_flg, criterion
   
   if (Argc < 4):
      print_usage(Argv[0])
   
   fname1    = Argv[1].strip()
   fname2    = Argv[2].strip()
   data_type = int(Argv[3].strip())

   iarg = 4
   while iarg < Argc:
      if   (argv[iarg].strip() == "-ratio"):
         ratio_flg = 1
         iarg += 1
      elif (argv[iarg].strip() == "-max_diff"):
         max_flg   = 1
         criterion = 999
         iarg += 1
      elif (argv[iarg].strip() == "-spin0"):
         spin_flg = 0
         iarg += 1
      elif (argv[iarg].strip() == "-spin1"):
         spin_flg = 1
         iarg += 1
      elif (argv[iarg].strip() == "-diff_rate"):
         if (iarg == Argc-1):
            print("One argument is needed after '-diff_rate' option.")
            quit(1)
         criterion = float(Argv[iarg+1].strip())
         max_flg   = 0
         iarg += 2
      else:
         print("ERROR: Invalid option, %s" % argv[iarg].strip())
         print
         print_usage(Argv[0])
         quit(1)

   fsize1 = os.path.getsize(fname1)
   fsize2 = os.path.getsize(fname2)
   if (fsize1 != fsize2 and data_type != 2):
      print("ERROR: The size of each file are different.")
      quit(1)

### Culculate the difference between two files ###
def diff_data(Data1, Data2, Count_Zero, Count_p, Count_m):
   if (Data1 == 0.0 and Data2 == 0.0):
      Count_Zero[0] += 1
   elif	(Data1 + Data2 == 0.0):
      Count_Zero[0] += 1
   else:
      if (ratio_flg == 0):
         diff = 2 * abs((Data1 - Data2) / (Data1 + Data2))
      else:
         diff = abs(Data1 / Data2)
      
      if (diff == 0.0):
         Count_Zero[0] += 1
      else:
         tmp = math.log10(diff)
         
         if (tmp > 0):
            Count_p[ int(tmp)] += 1
         else:
            Count_m[-int(tmp)] += 1

### Read & Difference calculation for 1/48 compressed NBSwave files ###
def read_bin_compressed(Fname1, Fname2):
   global DataSize
   
   ifile1 = open(Fname1, 'rb')
   ifile2 = open(Fname2, 'rb')
   
   dummy = struct.unpack('<i', ifile1.read(4))[0]
   if (dummy != 16):
      print("ERROR: The file '%s' may be not Ishii-san's compressed NBS data." % Fname1)
      quit(1)
   dummy = struct.unpack('<i', ifile2.read(4))[0]
   if (dummy != 16):
      print("ERROR: The file '%s' may be not Ishii-san's compressed NBS data." % Fname2)
      quit(1)
   dummy = ifile1.read(4*3)[0]
   dummy = ifile2.read(4*3)[0]
   length1 = struct.unpack('<i', ifile1.read(4))[0]
   length2 = struct.unpack('<i', ifile2.read(4))[0]
   if (length1 != length2):
      print("ERROR: The data length of each files is different.")
      quit(1)
   dummy = ifile1.read(4*(11+length1))[0]
   dummy = ifile2.read(4*(11+length2))[0]

   DataSize = length1

   tmp1_r = [0.0] * 16
   tmp1_i = [0.0] * 16
   tmp2_r = [0.0] * 16
   tmp2_i = [0.0] * 16
   for i in range(DataSize):
      for j in range(16):
         tmp1_r[j] = struct.unpack('<d', ifile1.read(8))[0]
         tmp1_i[j] = struct.unpack('<d', ifile1.read(8))[0]
         tmp2_r[j] = struct.unpack('<d', ifile2.read(8))[0]
         tmp2_i[j] = struct.unpack('<d', ifile2.read(8))[0]
      
      TmpData1 = tmp1_r[5] + tmp1_r[10] - tmp1_r[6] - tmp1_r[9] # Spin=0 projection (Not-normalize)
      TmpData2 = tmp2_r[5] + tmp2_r[10] - tmp2_r[6] - tmp2_r[9]
      diff_data(TmpData1, TmpData2, Count_Zero_r, Count_p_r, Count_m_r)

      TmpData1 = tmp1_i[5] + tmp1_i[10] - tmp1_i[6] - tmp1_i[9]
      TmpData2 = tmp2_i[5] + tmp2_i[10] - tmp2_i[6] - tmp2_i[9]
      diff_data(TmpData1, TmpData2, Count_Zero_i, Count_p_i, Count_m_i)
      
      TmpData1 = tmp1_r[5] + tmp1_r[10] + tmp1_r[6] + tmp1_r[9] # Spin=1 projection (Not-normalize)
      TmpData2 = tmp2_r[5] + tmp2_r[10] + tmp2_r[6] + tmp2_r[9]
      diff_data(TmpData1, TmpData2, Count_Zero_r_s1, Count_p_r_s1, Count_m_r_s1)
      
      TmpData1 = tmp1_i[5] + tmp1_i[10] + tmp1_i[6] + tmp1_i[9]
      TmpData2 = tmp2_i[5] + tmp2_i[10] + tmp2_i[6] + tmp2_i[9]
      diff_data(TmpData1, TmpData2, Count_Zero_i_s1, Count_p_i_s1, Count_m_i_s1)
   
   ifile1.close()
   ifile2.close()

### Read & Difference calculation for Not-compressed NBSwave files ###
def read_bin_Not_compressed(Fname1, Fname2):
   global DataSize
   
   Fsize = os.path.getsize(Fname1)
   if ((Fsize % 16) != 0):
      print("ERROR: Invalid file size, fsize % sizeof(complex double) != 0")
      quit(1)

   ifile1 = open(Fname1, 'rb')
   ifile2 = open(Fname2, 'rb')
   
   DataSize = Fsize / 16

   for i in range(DataSize):
      TmpData1 = struct.unpack('>d', ifile1.read(8))[0]
      TmpData2 = struct.unpack('>d', ifile2.read(8))[0]
      diff_data(TmpData1, TmpData2, Count_Zero_r, Count_p_r, Count_m_r)
      
      TmpData1 = struct.unpack('>d', ifile1.read(8))[0]
      TmpData2 = struct.unpack('>d', ifile2.read(8))[0]
      diff_data(TmpData1, TmpData2, Count_Zero_i, Count_p_i, Count_m_i)

   ifile1.close()
   ifile2.close()

### Read & Difference calculation for Text 2pt-correlator files ###
def read_text_2pt_corr(Fname1, Fname2):
   global DataSize
   
   ifile1 = open(Fname1, 'r')
   ifile2 = open(Fname2, 'r')
   
   lines1 = ifile1.readlines()
   lines2 = ifile2.readlines()
   
   ifile1.close()
   ifile2.close()
   
   if (len(lines1) != len(lines2)):
      print("ERROR: #.time of each files is different.")
      quit(1)
   if (len(lines1[0].split()) != len(lines2[0].split())):
      print("ERROR: #.mom of each files is different.")
      quit(1)
   
   Ntime = len(lines1)
   Nmom  = (len(lines1[0].split()) - 1) / 2
   
   DataSize = Ntime * Nmom
   
   for l in range(Ntime):
      line1 = lines1[l].split()
      line2 = lines2[l].split()
      for m in range(Nmom):
         diff_data(float(line1[1+2*m  ]), float(line2[1+2*m  ]), Count_Zero_r, Count_p_r, Count_m_r)
         diff_data(float(line1[1+2*m+1]), float(line2[1+2*m+1]), Count_Zero_i, Count_p_i, Count_m_i)

def show_results(Count_Zero, Count_p, Count_m):
   if (Count_Zero[0] != 0):
      print("      Zero  : %12d" % Count_Zero[0])

   for i in range(ArraySize-1, 0, -1):
      if (Count_p[i] != 0):
         print("      1e%-+4d: %12d" % (i, Count_p[i]))

   if (Count_p[0] + Count_m[0] != 0):
      print("      1e+0  : %12d" % (Count_p[0] + Count_m[0]))

   for i in range(1, ArraySize, 1):
      if (Count_m[i] != 0):
         print("      1e%-4d: %12d" % (-i, Count_m[i]))

def count_check(Count_Zero, Count_p, Count_m, Comments):
   count_sum = Count_Zero[0] + sum(Count_p) + sum(Count_m)
   if (count_sum != DataSize):
      print("ERROR: Data counts is not equarl to Data size")
      print("     : %d != %d (%s)" % (count_sum, DataSize, Comments))
      quit(1)

def max_results(Count_p, Count_m):
   for i in range(ArraySize-1, 0, -1):
      if (Count_p[i] != 0):
         return ("1e%-+d" % i)

   if (Count_p[0] + Count_m[0] != 0):
      return "1e+0"
   
   for i in range(1, ArraySize, 1):
      if (Count_m[i] != 0):
         return ("1e%-d" % -i)

   return "0"

def count_above_criterion(Count, Count_p, Count_m):
   for i in range(ArraySize-1, 0, -1):
      if (Count_p[i] != 0 and criterion < pow(10.0, i)):
         Count[0] += Count_p[i]

   if (Count_p[0] + Count_m[0] != 0 and criterion < 1.0):
      Count[0] += (Count_p[0] + Count_m[0])

   for i in range(1, ArraySize, 1):
      if (Count_m[i] != 0 and criterion < pow(10.0, -i)):
         Count[0] += Count_m[i]

###################################
############ Main part ############
if __name__ == "__main__":
   import sys
   import math
   import os.path
   import struct
   
   ### Set arguments
   argv = sys.argv
   argc = len(argv)
   
   set_args(argc, argv)

   ### Read data
   if   (data_type == 0):
      read_bin_compressed(fname1, fname2)
   elif (data_type == 1):
      read_bin_Not_compressed(fname1, fname2)
   elif (data_type == 2):
      read_text_2pt_corr(fname1, fname2)
   else:
      print_usage(argv[0])

   ### Check count
   count_check(Count_Zero_r, Count_p_r, Count_m_r, "Real")
   count_check(Count_Zero_i, Count_p_i, Count_m_i, "Imag")

   if (data_type == 0):
      count_check(Count_Zero_r_s1, Count_p_r_s1, Count_m_r_s1, "Real, s=1")
      count_check(Count_Zero_i_s1, Count_p_i_s1, Count_m_i_s1, "Imag, s=1")

   ### Show resutlts
   if (max_flg == 1):
      tmp_max_r = max_results(Count_p_r, Count_m_r)
      tmp_max_i = max_results(Count_p_i, Count_m_i)
      MaxValue1 = tmp_max_r if float(tmp_max_r) > float(tmp_max_i) else tmp_max_i

      if (data_type == 0 and spin_flg != 0):
         tmp_max_r = max_results(Count_p_r_s1, Count_m_r_s1)
         tmp_max_i = max_results(Count_p_i_s1, Count_m_i_s1)
         MaxValue2 = tmp_max_r if float(tmp_max_r) > float(tmp_max_i) else tmp_max_i
         
      if   (data_type == 0 and spin_flg == 2):
         print("%s %s" % (MaxValue1, MaxValue2))
      elif (data_type == 0 and spin_flg == 1):
         print("%s" % MaxValue2)
      else:
         print("%s" % MaxValue1)
   elif (criterion != 999):
      tmp_count = [0] * 1
      count_above_criterion(tmp_count, Count_p_r, Count_m_r)
      count_above_criterion(tmp_count, Count_p_i, Count_m_i)

      if (data_type == 0 and spin_flg != 0):
         tmp_count_s1 = [0] * 1
         count_above_criterion(tmp_count_s1, Count_p_r_s1, Count_m_r_s1)
         count_above_criterion(tmp_count_s1, Count_p_i_s1, Count_m_i_s1)

      if   (data_type == 0 and spin_flg == 2):
         print("%lf %lf" % (float(tmp_count[0])/float(2*DataSize), float(tmp_count_s1[0])/float(2*DataSize)))
      elif (data_type == 0 and spin_flg == 1):
         print("%lf" % (float(tmp_count_s1[0])/float(2*DataSize)))
      else:
         print("%lf" % (float(tmp_count[0])/float(2*DataSize)))
   else:
      if (data_type == 0):
         if (Count_Zero_r   [0] == DataSize and Count_Zero_i   [0] == DataSize and
             Count_Zero_r_s1[0] == DataSize and Count_Zero_i_s1[0] == DataSize):
            print("Totally Same")
            quit(0)
      else:
         if (Count_Zero_r[0] == DataSize and Count_Zero_i[0] == DataSize):
            print("Totally Same")
            quit(0)

      if (ratio_flg == 0):
         print("\n========= #.value of Relative error =========")
      else:
         print("\n============= #.value of ratio ==============")
      
      print
      
      if (data_type == 0 and spin_flg != 1):
         print("   <S=0, Sz=0>")

      if (data_type != 0 or (data_type == 0 and spin_flg != 1)):
         print("         (Real)")
         show_results(Count_Zero_r, Count_p_r, Count_m_r)
         print
         print("         (Imag)")
         show_results(Count_Zero_i, Count_p_i, Count_m_i)
         print
      
      if (data_type == 0 and spin_flg != 0):
         print("   <S=1, Sz=0>")
         print("         (Real)")
         show_results(Count_Zero_r_s1, Count_p_r_s1, Count_m_r_s1)
         print
         print("         (Imag)")
         show_results(Count_Zero_i_s1, Count_p_i_s1, Count_m_i_s1)
         print
         
      print "=============================================\n"
