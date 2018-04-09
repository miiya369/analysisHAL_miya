extern "C" 
{
#include <lime_config.h>
#include <lime.h>
#include <lime_fixed_types.h>
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gauge_link.h>

int ssize;
int tsize;

bool machine_is_little_endian();
void change_endian(double*, size_t);
void read_conf(const char*,double*);

int main(int argc, char** argv)
{
  if (argc != 5)
    {
      printf("usage : %s [ssize] [tsize] [limed-ildg file1] [limed-ildg file2]\n",argv[0]);
      return 1;
    }

  ssize = atoi(argv[1]);
  tsize = atoi(argv[2]);
  
  size_t csize = ssize * ssize * ssize * tsize * 4 * 3*3 * 2;

  double *conf1 = new double[csize];
  double *conf2 = new double[csize];

  read_conf(argv[3], conf1);
  read_conf(argv[4], conf2);

  double plaq1 = plaquette(conf1, ssize, tsize);
  double plaq2 = plaquette(conf2, ssize, tsize);

  printf(" plaquette difference = %e\n", fabs((plaq1-plaq2)/plaq1));

  delete [] conf1;
  delete [] conf2;

  return 0;
}

bool machine_is_little_endian() {
   
  int endianTEST = 1;
  if (*(char*)&endianTEST) return true;
  else                     return false;
}

void change_endian(double* DATA, size_t DATA_size)
{
//  printf("endian changing ... ");

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t k=0; k<DATA_size; k++)
    {
      char dummy[8];
      for (int j=0; j<8; j++) dummy[j] = ((char*)&DATA[k])[j];
      for (int j=0; j<8; j++) ((char*)&DATA[k])[j] = dummy[7-j];
    }
//  printf("DONE\n");
}

void read_conf(const char* fname, double* conf)
{
  LimeReader *reader;
  int        status;
  n_uint64_t nbytes, read_bytes; 
  char*      lime_type;

  FILE* fp;
  if((fp = fopen(fname, "r"))==NULL){
    fprintf(stderr,
            "cannot fopen '%s' for reading\n", fname);
    exit(1);
  }

  if((reader = limeCreateReader(fp)) == NULL){
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(1);
  }

  while((status = limeReaderNextRecord(reader)) != LIME_EOF){
    if(LIME_SUCCESS != status){
      fprintf(stderr, "limeReaderNextRecord returns status = %d\n", status);
      exit(1);
    }

    nbytes    = limeReaderBytes   (reader);
    lime_type = limeReaderType    (reader);

    if(strcmp(lime_type,"ildg-binary-data")==0)
      {
	if(2 *3*3*4 *ssize *ssize *ssize *tsize *sizeof(double) != nbytes){
	  fprintf(stderr,
		  "ERROR(misc::read_ildg_config): inconsistent data size ! %ld != %ld\n",
		  nbytes,
		  2*3*3*4 *ssize *ssize *ssize *tsize *sizeof(double));
	  exit(1);
	}

	read_bytes = nbytes;
	status = limeReaderReadData((void*)conf, &read_bytes, reader);

	if(status < 0){
	  if(status != LIME_EOR){
	    fprintf(stderr,
		    "LIME read error ccurred in reading ildg-binary-data: "
		    "status = %d %llu bytes wanted, %llu read\n",
		    status,
		    (unsigned long long)nbytes,
		    (unsigned long long)read_bytes);
	    exit(1);
	  }
	}

	if(machine_is_little_endian()){
	  change_endian((double*)conf, read_bytes/sizeof(double));
	}
      }
  }
  limeDestroyReader(reader);  
  fclose(fp);
}
