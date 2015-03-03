#include "gmock/gmock.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "lime_reader.h"

#ifdef __cplusplus
}
#endif

TEST( ALimeReader, ReadContent )
{


	LimeReader* reader;
	FILE *fp;

	fp = fopen( "qcdsf.721.00100.lime" , "r");

	reader = limeCreateReader(fp);

	int status;
	n_uint64_t nbytes, read_bytes;
	int msg = 0;
	int rec = 0;
	int first = 0;
	char *lime_type;
	size_t bytes_pad;
	int MB_flag, ME_flag;

	while( (status = limeReaderNextRecord(reader)) != LIME_EOF ){

	    if( status != LIME_SUCCESS ) {
	      fprintf(stderr, "limeReaderNextRecord returned status = %d\n",
		      status);
//	      return EXIT_FAILURE;
	    }

	    nbytes    = limeReaderBytes(reader);
	        lime_type = limeReaderType(reader);
	        bytes_pad = limeReaderPadBytes(reader);
	        MB_flag   = limeReaderMBFlag(reader);
	        ME_flag   = limeReaderMEFlag(reader);

	        if (MB_flag == 1 || first)
	          {
	    	first = 0;
	    	rec = 0;
	    	msg++;
	          }

	        rec++;

	        printf("\n\n");
	        printf("Message:        %d\n", msg);
	        printf("Record:         %d\n", rec);
	        printf("Type:           %s\n", lime_type);
	        printf("Data Length:    %llu\n", (unsigned long long)nbytes);
	        printf("Padding Length: %lu\n", (unsigned long)bytes_pad);
	        printf("MB flag:        %d\n", MB_flag);
	        printf("ME flag:        %d\n", ME_flag);
	}
}
