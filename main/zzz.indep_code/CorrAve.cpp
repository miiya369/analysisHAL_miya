//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup none
 * @brief   Take average for Text 2pt-correlator data
 * @author  Takaya Miyamoto
 * @since   Tue Jun 21 21:39:44 JST 2016
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

static string *ipath = NULL;
static string  opath;
static int     Naverage;

static int Tsize;

inline int set_args(int, char**);

inline int input (double*, int);
inline int output(const double*);

//========================================================================//
int main(int argc, char **argv)
{
   if (set_args(argc, argv) == 1) return 1;
      
   double *Corr = new double[Tsize * 10]; // 10 -> Momentum=5 * Real,Imag
   for (int i=0; i<Tsize * 10; i++) Corr[i] = 0.0;
   
   for (int i=0; i<Naverage; i++) if (input(Corr, i) == 1) return 1;
   
   if (output(Corr) == 1) return 1;
   
   delete [] Corr;
   if (ipath  != NULL) delete [] ipath;
   return 0;
}
//========================================================================//

inline int set_args(int argc, char** argv)
{
   if (argc < 3) {
      printf("usage: %s [ofile] [ifile1] [ifile2] ...\n", argv[0]);
      return 1;
   }
   
   Naverage = argc - 2;
   ipath    = new string [Naverage];
   
   for (int i=0; i<Naverage; i++) {
      ipath[i] = argv[i+2];
   }
   opath = argv[1];
   
   FILE* fp;
   if ((fp = fopen(ipath[0].c_str(), "r")) == NULL) {
      printf("The file '%s' is not found.\n", ipath[0].c_str());
      return 1;
   }
   int count = 0;
   char tmp[1000];
   while(fgets(tmp, 1000, fp) != NULL) count++;
   fclose(fp);

   Tsize = count;
   
//   printf("Nave = %d\n", Naverage); // For Debug
//   return 1;
   
   return 0;
}

inline int input(double* corr, int data_num)
{ 
   FILE* fp;
   if ((fp = fopen(ipath[data_num].c_str(), "rb")) == NULL) {
      printf("The file '%s' is not found.\n", ipath[data_num].c_str());
      return 1;
   }

   int    t_tmp;
   double tmp[10];
   for (int it=0; it<Tsize; it++) {
      if (fscanf(fp, " %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &t_tmp
		 , &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4]
		 , &tmp[5], &tmp[6], &tmp[7], &tmp[8], &tmp[9] ) != 11) {
	 printf("ERROR: Invalid File '%s'\n", ipath[data_num].c_str());
	 return 1;
      }
      if (it != t_tmp) {
	 printf("ERROR: Invalid File '%s'\n", ipath[data_num].c_str());
	 return 1;
      }
      for (int loop=0; loop<10; loop++) corr[loop+10*it] += tmp[loop] / double(Naverage);
   }
   fclose(fp);

   return 0;
}

inline int output(const double* CORR)
{
   FILE* fp;
   if ((fp = fopen(opath.c_str(), "wb")) == NULL) {
      printf("The file '%s' can not be open.\n", opath.c_str());
      return 1;
   }
   
   for (int it=0; it<Tsize; it++)
     if ( fprintf(  fp, "%d %1.16e %1.16e %1.16e %1.16e "
		    "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", it
		    , CORR[0+10*it], CORR[1+10*it], CORR[2+10*it]
		    , CORR[3+10*it], CORR[4+10*it], CORR[5+10*it]
		    , CORR[6+10*it], CORR[7+10*it]
		    , CORR[8+10*it], CORR[9+10*it] ) < 0 ) {
       printf("invalid writing for file '%s'\n", opath.c_str());
       return 1;
     }
   fclose(fp);
   
   return 0;
}
