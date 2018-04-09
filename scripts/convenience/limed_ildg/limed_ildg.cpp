extern "C" 
{
#include <lime_config.h>
#include <lime.h>
#include <lime_fixed_types.h>
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      printf("usage : %s [limed-ildg file]\n",argv[0]);
      return 1;
    }

  LimeReader *reader;
  int        status;
  n_uint64_t nbytes, read_bytes;
  int        msg, rec, first;
  char*      lime_type;
  size_t     bytes_pad;
  int        MB_flag, ME_flag;
  char       lfn[64000];

  FILE* fp;
  if((fp = fopen(argv[1], "r"))==NULL){
    fprintf(stderr,
            "cannot fopen '%s' for reading\n", argv[1]);
    return 1;
  }

  if((reader = limeCreateReader(fp)) == NULL){
    fprintf(stderr, "Unable to open LimeReader\n");
    return 1;
  }

  msg = 0; first = 1; rec = 0;
  while((status = limeReaderNextRecord(reader)) != LIME_EOF){
    if(LIME_SUCCESS != status){
      fprintf(stderr, "limeReaderNextRecord returns status = %d\n", status);
      return 1;
    }

    nbytes    = limeReaderBytes   (reader);
    lime_type = limeReaderType    (reader);
    bytes_pad = limeReaderPadBytes(reader);
    MB_flag   = limeReaderMBFlag  (reader);
    ME_flag   = limeReaderMEFlag  (reader);

    if (MB_flag == 1 || first)
      {
        first = 0;
        rec = 0;
        msg++;
      }

    rec++;

    if(true){
      printf("Message:        %d\n",   msg);
      printf("Record:         %d\n",   rec);
      printf("Type:           %s\n",   lime_type);
      printf("Data Length:    %llu\n", (unsigned long long)nbytes);
      printf("Padding Length: %lu\n",  bytes_pad);
      printf("MB flag:        %d\n",   MB_flag);
      printf("ME flag:        %d\n",   ME_flag);
      printf("\n");
    }

    //---------------------------------------------------------------------
    // reading ildg-binary-data
    //---------------------------------------------------------------------
    if(strcmp(lime_type,"ildg-binary-data")==0)
      {
	printf("@ Read ildg binary data ... SKIP\n\n");
      }
    //---------------------------------------------------------------------
    // reading ildg-format
    //---------------------------------------------------------------------
    if(strcmp(lime_type,"ildg-format")==0)
      {
	char* buf = new char[nbytes+1];
	read_bytes = nbytes;
	status = limeReaderReadData((void*)buf, &read_bytes, reader);

	if(status < 0)
	  if(status != LIME_EOR)
	    {
	      fprintf(stderr,
		      "LIME read error ccurred in reading ildg-format: "
		      "status = %d %llu bytes wanted, %llu read\n",
		      status,
		      (unsigned long long)nbytes,
		      (unsigned long long)read_bytes);
	      return 1;
	    }
	
	printf("@ Read ildg format ...\n\n%s\n\n", buf);
      }
    //---------------------------------------------------------------------
    // reading ildg-data-lfn
    //---------------------------------------------------------------------
    if(strcmp(lime_type,"ildg-data-lfn")==0)
      {
	char* buf = new char[nbytes+1];
	read_bytes = nbytes;
	status = limeReaderReadData((void*)buf, &read_bytes, reader);
	
	if(status < 0){
	  if(status != LIME_EOR){
	    fprintf(stderr,
		    "LIME read error ccurred in reading ildg-data-lfn: "
		    "status = %d %llu bytes wanted, %llu read\n",
		    status,
		    (unsigned long long)nbytes,
		    (unsigned long long)read_bytes);
	    return 1;
	  }
	}
	snprintf(lfn, sizeof(lfn), "%s", buf);
	printf("@ Read ildg data lfn ...\n\n%s\n\n", lfn);
      }
  }
  
  limeDestroyReader(reader);
  
  fclose(fp);

  return 0;
}
