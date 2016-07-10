/*----------------------------------------------------------------------
    This file is part of aaphoto.

    aaphoto is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    aaphoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------*/



/* ----------------------------------------------------- */
/* ----------- STRING CONVERT TO LOWER CASE ------------ */
/* ----------------------------------------------------- */
void STRING_CONVERT_TO_LCASE(char *strin, char *strout)
{
	int i;
	for (i=0; strin[i]!=0; i++){
		if (i < max_char-1){
			strout[i]=tolower(strin[i]); }
	}
	strout[i]=0;
}




/* -------------------------------------------------- */
/* ----------- STRING CONVERT TO INTEGER ------------ */
/* -------------------------------------------------- */
int STRING_CONVERT_TO_INTEGER(char *str, int *number)
{
	int i, c, d;
	int xx = 0;
	int yy;
	int num;
	int num_max = 8;
	i = 0;
	while ((str[i]) && (i < num_max)){ i++; }
	if (i >= num_max) return 1;
	c = i;
	num = 0;
	if (c > 0){
		d = 1;
		for (i=1; i<=c; i++){

			yy = 1;
			if (str[c-i] == '0') { xx = 0; yy = 0; }
			if (str[c-i] == '1') { xx = 1; yy = 0; }
			if (str[c-i] == '2') { xx = 2; yy = 0; }
			if (str[c-i] == '3') { xx = 3; yy = 0; }
			if (str[c-i] == '4') { xx = 4; yy = 0; }
			if (str[c-i] == '5') { xx = 5; yy = 0; }
			if (str[c-i] == '6') { xx = 6; yy = 0; }
			if (str[c-i] == '7') { xx = 7; yy = 0; }
			if (str[c-i] == '8') { xx = 8; yy = 0; }
			if (str[c-i] == '9') { xx = 9; yy = 0; }
			if (yy) return 1;

			num += xx * d;
			d *= 10;
		}
	}
	else{ return 1; }

	*number = num;
	return 0;
}




/* ------------------------------------- */
/* ----------- STRING PRINT ------------ */
/* ------------------------------------- */
/* print text to stdout */
void STRING_PRINT(const char str[])
{
	if (!(opt_quiet)){
		fprintf(stdout, "%s", str);
		fflush(stdout);
	}
}




/* ------------------------------------- */
/* -------- STRING PRINT ERROR --------- */
/* ------------------------------------- */
/* print text to stderr */
void STRING_PRINTE(const char str[])
{
/*	if (!(opt_quiet)){ */
		fprintf(stderr, "%s", str);
		fflush(stderr);
/*	} */
}




/* ------------------------------------- */
/* ------- STRING PRINT VERBOSE -------- */
/* ------------------------------------- */
/* print time and text to stdout if in verbose mode */
void STRING_PRINTV(const char str[])
{
    if(opt_verbose){
	if (!(opt_quiet)){
		/* print time since start in seconds if bigger than 1 */
		/* and decrement the time by 1 because the program may start in a middle of a second, that i don't know */
		int t = time(NULL) - mytime;
		if (t > 1){ printf("(%d) ", t-1); }
		/* print text */
		fprintf(stdout, "%s", str);
		fflush(stdout);
	}
    }
}

void STRING_PRINTV2(const char str[])
{
    if(opt_verbose){
	if (!(opt_quiet)){
		/* print text */
		fprintf(stdout, "%s", str);
		fflush(stdout);
	}
    }
}




/* -------------------------------------- */
/* ----------- STRING PRINTVD ----------- */
/* -------------------------------------- */
/* print decimal number to stdout */
void STRING_PRINTVD(int num)
{
	if(opt_verbose){
		if (!(opt_quiet)){
       			fprintf(stdout, "%d", num);
			fflush(stdout);
		}
	}
}




/* -------------------------------------- */
/* ----------- STRING PRINTED ----------- */
/* -------------------------------------- */
/* print decimal number to stderr */
void STRING_PRINTED(int num)
{
	if (!(opt_quiet)){
		fprintf(stderr, "%d", num);
		fflush(stderr);
	}
}




/* -------------------------------------- */
/* ----------- STRING PRINTF ------------ */
/* -------------------------------------- */
/* print floating number to stdout */
void STRING_PRINTF(double num)
{
	if (!(opt_quiet)){
       		fprintf(stdout, "%.2f", num);
		fflush(stdout);
	}
}




/* --------------------------------------- */
/* ----------- STRING COMPARE ------------ */
/* --------------------------------------- */
/* compare two string values */
/* result is 0 if true */
int STRING_COMPARE(char *str1, char *str2)
{
	int i = 0;
	int result = 0;
	while ((str1[i]) && (str2[i])) {
		if (str1[i] != str2[i]) result = 1;
		i++;
	}
	if (str1[i] != str2[i]) result = 1;
	return result;
}




/* --------------------------------------------------------- */
/* ----------- STRING COMPARE WITH FIXED LENGTH ------------ */
/* --------------------------------------------------------- */
int STRING_COMPARE_FIX(char *str1, char *str2, int count)
{
	int result = 0;
	int i;
	for (i=0; i<count; i++)
		if (str1[i] != str2[i]) result = 1;
	return result;
}




/* -------------------------------------- */
/* ----------- GET FILE NAME ------------ */
/* -------------------------------------- */
/* get file name from path */
void GET_FILE_NAME(char *strin, char *strout)
{
	/* check length of string */
	long len;
	long i, c;
	for (i=0; strin[i]!='\0'; i++);
	len = i;
	/* check only path */
	c = 0;
	i = len-1;
	while ((i >= 0) && (strin[i] != slsh)) { c++; i--; }
	/*if (i >= 0) { */
		for (i=0; i<c; i++) {
			if (i < max_char) {
				strout[i] = strin[i+len-c]; }
		}
		strout[i] = 0;
	/*} */
}




/* ------------------------------------------- */
/* ----------- GET FILE NAME ONLY ------------ */
/* ------------------------------------------- */
/* get only the file name of the path */
void GET_FILE_NAME_ONLY(char *strin, char *strout)
{
	/* check length of string */
	long len;
	long i, c, st;
	for (i=0; strin[i]!='\0'; i++);
	len = i;
	/* check only path */
	i = len-1;
	while ((i >= 0) && (strin[i] != slsh)) { i--; }
    st = i + 1;
	/* check extension of file */
	c = 0;
	i = len-1;
	while ((i >= 0) && (strin[i] != '.')) { i--; }
	c = i - 1;
	if (c >= st) {
		for (i=0; i<c-st+1; i++) {
			if (i < max_char-1) {
				strout[i] = strin[i+st]; }
			}
		strout[i] = 0;
	}
	else {
		strout[0] = 0;
	}
}




/* ------------------------------------------- */
/* ----------- GET FILE EXTENSION ------------ */
/* ------------------------------------------- */
/* get extenstion of the file in path */
void GET_FILE_EXTENSION(char *strin, char *strout)
{
	/* check length of string */
	long len;
	long i, c;
	for (i=0; strin[i]!='\0'; i++);
	len = i;
	/* check extension of file */
	c = 0;
	i = len-1;
	while ((i >= 0) && (strin[i] != '.')) { c++; i--; }
	c++;
	if (c > 0) {
		for (i=0; i<c; i++) {
			if (i < max_char-1) {
				strout[i] = strin[i+len-c]; } }
		strout[i] = 0;
	}
	else {
		strout[0] = 0;
	}
}




/* ------------------------------------------ */
/* ------------- GET FILE PATH -------------- */
/* ------------------------------------------ */
/* get file path */
void GET_FILE_PATH(char *strin, char *strout)
{
	/* check length of string */
	long len;
	long i, c;
	for (i=0; strin[i]!='\0'; i++);
	len = i;
	/* check only path */
	c = 0;
	i = len-1;
	while ((i >= 0) && (strin[i] != slsh)) { c++; i--; }
	if (i >= 0) {
		for (i=0; i<len-c; i++) {
			if (i < max_char-1) {
				strout[i] = strin[i]; }
		}
		strout[i] = 0;
	}
	else {
		strout[0] = 0;
	}
}




/* ---------------------------------------- */
/* ------------ FILE EXIST? --------------- */
/* ---------------------------------------- */
int FILE_EXIST(char *file_name)
{
	FILE *fhandle;
	fhandle = fopen(file_name, "rb");
	if (fhandle == 0) return 0;
	fclose(fhandle);
	return 1;
}




/* ---------------------------------------------- */
/* ------------- GET FILE NAME NEW -------------- */
/* ---------------------------------------------- */
/* put together the new output file name */
int GET_FILE_NAME_NEW(char *strin, char *strout)
{
	char fpath [max_char];
	char fname [max_char];
	char fext  [max_char];
	char fextj1 [] = ".jpg\0";
	char fextj3 [] = ".png\0";
	char fextj4 [] = ".bmp\0";
	char fnew   [] = "_new";
	long i, c;

	GET_FILE_PATH(strin, fpath);
	GET_FILE_NAME_ONLY(strin, fname);
	GET_FILE_EXTENSION(strin, fext);

	c = 0;

	if (!(opt_output))
	{
		i = 0;
		while (fpath[i] != '\0')
		{
			if (c >= max_char) return 255;
			strout[c] = fpath[i];
			i++; c++;
		}
	}
	else
	{
		i = 0;
		while (opt_output_path[i] != '\0')
		{
			if (c >= max_char) return 255;
			strout[c] = opt_output_path[i];
			i++; c++;
		}
	}

	i = 0;
	while (fname[i] != '\0')
	{
		if (c >= max_char) return 255;
		strout[c] = fname[i];
		i++; c++;
	}

	if (!(opt_overwrite))
	{
		i = 0;
		while (fnew[i] != '\0')
		{
			if (c >= max_char) return 255;
			strout[c] = fnew[i];
			i++; c++;
		}
	}

	i = 0;
	while (fext[i] != '\0')
	{
		if (c >= max_char) return 255;
		if ((!opt_jpg) && (!opt_png) && (!opt_bmp)) strout[c] = fext[i];
		if (opt_jpg) strout[c] = fextj1[i];
		if (opt_png) strout[c] = fextj3[i];
		if (opt_bmp) strout[c] = fextj4[i];
		i++; c++;
	}
	if (c >= max_char) return 255;
	strout[c] = fext[i];

	if (opt_jpg) bitmap_format_jpg_file_type = 1;

	if (!(opt_overwrite) && (FILE_EXIST(strout))) return 1;

	return 0;

}




/* -------------------------------------------- */
/* ------------- GET FILE FORMAT -------------- */
/* -------------------------------------------- */
int GET_FILE_FORMAT(char *file_name)
{
	int res;
	char fext[max_char];
	char fextl[max_char];

	/* format codes */
	/* --------------------------- */
	/* 0 - bmp */
	/* 1 - jpg */
	/* 2 - png */

	GET_FILE_EXTENSION(file_name, fext);

	STRING_CONVERT_TO_LCASE(fext, fextl);

	res = -1;

	if (!STRING_COMPARE(fextl, ".bmp"))  res = 0;
	if (!STRING_COMPARE(fextl, ".jpg"))  res = 1;
	if (!STRING_COMPARE(fextl, ".jpeg")) res = 1;
	if (!STRING_COMPARE(fextl, ".jpe"))  res = 1;
	if (!STRING_COMPARE(fextl, ".png"))  res = 2;

	return res;
}




/* ----------------------------------------- */
/* ----------- FILE LIST ADD --------------- */
/* ----------------------------------------- */
int FILE_LIST_ADD(char *file_name)
{
	/* isn't the file name buffer full yet? */
	if (file_name_buffer_pointer + max_char < max_file_name_buffer){

/*
        char fpath [max_char];
		char fname [max_char];
		char fext  [max_char];
		GET_FILE_PATH(file_name, fpath);
		GET_FILE_NAME_ONLY(file_name, fname);
		GET_FILE_EXTENSION(file_name, fext);
*/

        DIR *                   dp;
        DIR *                   dp2;
        const struct dirent *   ent;
        int                     cnt;

	/* if the name points to a directory */
        dp = opendir(file_name);
        if (dp != NULL){

            /* list files in directory */
            cnt = 0;
            while (ent = readdir(dp), ent != NULL)
            {
		int i, i2, res, flag;
            	char file_name_new [max_char];

                res = 0;
	            if (!STRING_COMPARE((char*)(ent->d_name), "."))  res = 1;
	            if (!STRING_COMPARE((char*)(ent->d_name), "..")) res = 1;
                if (res == 0){

                cnt++;

		/* create new file name with path */
       		i = 0;
       		while (file_name[i]){
		    if (i >= max_char) return 255;
                    file_name_new[i] = file_name[i]; i++; }

		/* remove '/' characters from the end of directory names when more than 1 */
                flag = 0;
                while ((i > 1) && (flag == 0)) {
                    if (file_name_new[i-1] == slsh) { i--; }
                    else { flag = 1; } }
		/* add '/' character to directory path */
                file_name_new[i] = slsh;

		/* add found file name to directory path */
                i++;
          	i2 = 0;
           	while (ent->d_name[i2]){
                    file_name_new[i] = ent->d_name[i2]; i++; i2++; }
                file_name_new[i] = 0;


		/* if new name is not a dir, then store file name */
                dp2 = opendir(file_name_new);
                if (dp2 != NULL){ closedir(dp2); }
                else{
			int i, len;

            		i = 0;
            		while (file_name_new[i]){ i++; }
        	    	len = i;
        			for (i=0; i<len; i++){
        				file_name_buffer[file_name_buffer_pointer] = file_name_new[i];
        				file_name_buffer_pointer++;
        			}
        			file_name_buffer[file_name_buffer_pointer] = '\0';
        			file_name_buffer_pointer++;
        			file_name_counter++;

                    /* printf("%s\n", file_name_new); */
                }

                }

            }
            closedir(dp);
            dp = NULL;
        }

	/* the name is a file, so i store it */
        else{
		int i, len;

    		i = 0;
    		while (file_name[i]){ i++; }
	    	len = i;
			for (i=0; i<len; i++){
				file_name_buffer[file_name_buffer_pointer] = file_name[i];
				file_name_buffer_pointer++;
			}
			file_name_buffer[file_name_buffer_pointer] = '\0';
			file_name_buffer_pointer++;
			file_name_counter++;
        }

    }
    else { return 1; }
	return 0;

}




/* -----------------------------------
             --- EXIF ---
   -----------------------------------
    char *exif_buffer;
    long  exif_buffer_length;
    long  exif_file_length;
    int   exif_flag;

    http://en.wikipedia.org/wiki/JPEG
    http://www.media.mit.edu/pia/Research/deepview/exif.html
  ------------------------------------
  ------------------------------------
*/


/* -------------------------------------- */
/* ----------- EXIF CLEAR --------------- */
/* -------------------------------------- */
/* clear exif information and deallocate memory */
int EXIF_CLEAR()
{
    if (exif_flag) {
	/* print info */ STRING_PRINTV("freeing exif buffer\n");
        free(exif_buffer);
        exif_buffer_length = 0;
        exif_flag = 0;
    }
    return 0;
}




/* ------------------------------------ */
/* ----------- EXIF GET --------------- */
/* ------------------------------------ */
/* read exif information from file and store it in memory */
int EXIF_GET(char *file_name)
{
    int exif_ok, exif_bad;
    long exif_start, exif_offset;
    int ch1, ch2;
    FILE *fhandle;

    exif_flag = 0;

    /* if --noexif flag was specified, then i don't store exif from image */
    if (opt_noexif) { return 1; }

    /* print info */ STRING_PRINTV("checking exif info\n");

    /* open file for reading */
    fhandle = fopen(file_name, "rb");
    if (fhandle == 0) return 1;
    /* determine the length of file */
    fseek(fhandle, 0, SEEK_END);
    exif_file_length = ftell(fhandle);
    fseek(fhandle, 0, SEEK_SET);


    /* examine exif information and check its length */
    exif_start = 0;
    exif_offset = 0;
    ch1 = 0;
    ch2 = 0;

    /* FFD8 JPEG indicator */
    fseek(fhandle, 0, SEEK_SET);
    ch1 = fgetc(fhandle);
    ch2 = fgetc(fhandle);
    if ((ch1 != 0xff) || (ch2 != 0xd8)) {
    	fclose(fhandle);
        return 1;
    }

    /* seek for beginning of exif info (start offset = 2) */
    exif_start = 2;
    exif_ok = 0;
    exif_bad = 0;
    while ((exif_ok == 0) && (exif_bad == 0)) {

        /* FFE1 Exif indicator */
        fseek(fhandle, exif_start + 0, SEEK_SET);
        ch1 = fgetc(fhandle);
        ch2 = fgetc(fhandle);
	/* if indicator does not start with FF */
	/* (the id tag of the next block within JPEG format) */
	/* then exit, because this is an error */
        if (ch1 != 0xff) { exif_bad = 1; }
        else {
	    /* search for FFE1 exif array indicator */
            if ((ch1 != 0xff) || (ch2 != 0xe1)) {

		/* check length of exif array */
                fseek(fhandle, exif_start + 2, SEEK_SET);
                ch1 = fgetc(fhandle);
                ch2 = fgetc(fhandle);
                exif_offset = (long)(ch1) * 256 + (long)(ch2);
                exif_start += 2 + exif_offset;
		/* exit if pointer reaches the end of file */
		if (exif_start >= exif_file_length) { exif_bad = 1; }

            }
            else {

                exif_ok = 1;

		/* check 'Exif00' pattern 45 78 69 66 00 00 */
                fseek(fhandle, exif_start + 4, SEEK_SET);
                ch1 = fgetc(fhandle);
                ch2 = fgetc(fhandle);
                if ((ch1 != 0x45) || (ch2 != 0x78)) {
                    exif_ok = 0;
                }
                fseek(fhandle, exif_start + 6, SEEK_SET);
                ch1 = fgetc(fhandle);
                ch2 = fgetc(fhandle);
                if ((ch1 != 0x69) || (ch2 != 0x66)) {
                    exif_ok = 0;
                }
                fseek(fhandle, exif_start + 8, SEEK_SET);
                ch1 = fgetc(fhandle);
                ch2 = fgetc(fhandle);
                if ((ch1 != 0x00) || (ch2 != 0x00)) {
                    exif_ok = 0;
                }

                /* check length of exif array */
                if (exif_ok == 1) {
                    fseek(fhandle, exif_start + 2, SEEK_SET);
                    ch1 = fgetc(fhandle);
                    ch2 = fgetc(fhandle);
		    /* exif length = exif pointer + 2 */
		    /* with the extra 2 bytes we take the FFE1 exif marker too into the buffer */
                    exif_buffer_length = (long)(ch1) * 256 + (long)(ch2) + 2;
                }
		else {
			exif_bad = 1;
		}

            }
        }

    }

    if (exif_ok == 0) {
    	fclose(fhandle);
        return 1;
    }

    /* allocate memory for file load */
    /* print info */ STRING_PRINTV("allocating memory for exif buffer\n");
    exif_buffer = malloc(exif_buffer_length * sizeof (*exif_buffer));
    if (exif_buffer == 0) {
        fclose(fhandle);
        return 1;
    }

    /* load exif part of file into memory */
    /* print info */ STRING_PRINTV("reading exif info\n");
    fseek(fhandle, exif_start, SEEK_SET);
    if (fread(exif_buffer, 1, exif_buffer_length, fhandle) == 0) {
    	fclose(fhandle);
        return 1;
    }

    /* close file */
    fclose(fhandle);
    exif_flag = 1;
    return 0;
}




/* ------------------------------------ */
/* ----------- EXIF PUT --------------- */
/* ------------------------------------ */
/* put exif information back to the file from memory (if there was any) */
int EXIF_PUT(char *file_name)
{
    if (exif_flag) {
    	FILE *fhandle;

	/* print info */ STRING_PRINTV("writing exif info ");
	/* print info */ STRING_PRINTVD((int)(exif_buffer_length));
	/* print info */ STRING_PRINTV2(" bytes\n");

	/* open file for writing */
    	fhandle = fopen(file_name, "wb+");
    	if (fhandle == 0) return 1;
	/* write FFD8 JPEG indicator into the first 2 bytes */
        if (fputc(0xff, fhandle) == 0) return 1;
        if (fputc(0xd8, fhandle) == 0) return 1;
	/* write out the rest of the exif info */
	/* here the exif info has to be written with 2 bytes less from the end */
	/* because the JPEG writer wants to add the FFD8 marker himself too */
	/* so this prevents FFD8 to be 2 times wrongly */
        if (fwrite(exif_buffer, 1, exif_buffer_length - 2, fhandle) == 0) return 1;
        fclose(fhandle);

    }
    return 0;
}




/* ------------------------------------ */
/* ----------- EXIF CORRECT ----------- */
/* ------------------------------------ */
/* correct the last 2 bytes of the exif data back from FFD8 to FFD9 */
/* cause when we put the exif back to the file, we put 2 bytes less */
/* and the jpeg writer can append the jpeg file starting with FFD8 */
/* so these 2 bytes needs to be changed back to FFD9 */
int EXIF_CORRECT(char *file_name)
{
    if (exif_flag) {

	/* open file for writing */
    	FILE *fhandle;
    	fhandle = fopen(file_name, "rb+");
    	if (fhandle == 0) return 1;
	/* write FFD9 EXIF END indicator to the end of exif data */
	fseek(fhandle, exif_buffer_length + 0, SEEK_SET);
        if (fputc(0xff, fhandle) == 0) return 1;
        if (fputc(0xd9, fhandle) == 0) return 1;
        fclose(fhandle);

    }
    return 0;
}




/* -------------------------------------------------- */
/* ----------- BITMAP READ IN BMP FORMAT ------------ */
/* -------------------------------------------------- */
int BITMAP_READ_BMP(char *file_name)
{
	unsigned long f_bm;
	unsigned long f_bitcount;
	unsigned long f_compressed;
	unsigned long f_headersize;
	unsigned long f_offs;
	unsigned long f_width;
	unsigned long f_height;
	unsigned long f_xpixelpermeter;
	unsigned long f_ypixelpermeter;
    	unsigned long addr, addr2;
    	unsigned long addr_offset;
	unsigned long bw, bh;
	unsigned long file_length;
	unsigned char *file_buffer;
	unsigned long i, x, y;
	FILE *fhandle;

	/* print info */ STRING_PRINTV("bmp initializations\n");

	/* open file for reading */
	/* print info */ STRING_PRINTV("opening file for reading\n");
	fhandle = fopen(file_name, "rb");
	if (fhandle == 0) return 1;
	/* check length of file */
	fseek(fhandle, 0, SEEK_END);
	file_length = ftell(fhandle);
	fseek(fhandle, 0, SEEK_SET);
	/* allocate memory for file load */
	/* print info */ STRING_PRINTV("allocating memory for bmp object\n");
	file_buffer = malloc(file_length * sizeof (*file_buffer));
	if (file_buffer == 0) return 1;
	/* load file into memory */
	/* print info */ STRING_PRINTV("reading image file\n");
	if (fread(file_buffer, 1, file_length, fhandle) == 0) {
		/* print info */ STRING_PRINTV("freeing bmp object\n");
		free(file_buffer); return 1; }
	/* close file */
	fclose(fhandle);


	/* read BMP indicator */
	f_bm = 0;
	f_bm += file_buffer[0] << 0;
	f_bm += file_buffer[1] << 8;
	f_bitcount = 0;
	f_bitcount += file_buffer[28] << 0;
	f_bitcount += file_buffer[29] << 8;
	f_compressed = 0;
	f_compressed += file_buffer[30] << 0;
	f_compressed += file_buffer[31] << 8;
	f_compressed += file_buffer[32] << 16;
	f_compressed += file_buffer[33] << 24;

	/* check BMP format (BMP header + uncompressed + 24 bit colors) */
	if ((f_bm == 0x00004d42) && ((f_bitcount == 24) || (f_bitcount == 8)) && (f_compressed == 0)) {
        bitmap_format_bmp_clrspc_type = f_bitcount;

        if (f_bitcount == 8){
		int gray_flag;

		/* offset pointing to color palette */
	    	f_headersize = 14;
	    	f_headersize += file_buffer[14] << 0;
	    	f_headersize += file_buffer[15] << 8;
	    	f_headersize += file_buffer[16] << 16;
	    	f_headersize += file_buffer[17] << 24;
		/* check whether the 8 bit image contains only gray colors? */
		gray_flag = 1;

		#ifdef __OPENMP__
		#pragma omp parallel for num_threads(max_threads)
		#endif
		for (i=0; i<256; i++){
			if ((file_buffer[i*4+f_headersize+0] != file_buffer[i*4+f_headersize+1]) ||
				(file_buffer[i*4+f_headersize+1] != file_buffer[i*4+f_headersize+2])){ gray_flag = 0; }
		}
		/* if not gray then exit, because minimim bitdepth of colors to correct is 24 (or 8 bit gray) */
		if (gray_flag == 0) {
			/* print info */ STRING_PRINTV("freeing bmp object\n");
			free(file_buffer); return 1; }
		}

		/* offset pointing to RGB datas */
		f_offs = 0;
		f_offs += file_buffer[10] << 0;
		f_offs += file_buffer[11] << 8;
		f_offs += file_buffer[12] << 16;
		f_offs += file_buffer[13] << 24;
		/* width of image in pixels */
		f_width = 0;
		f_width += file_buffer[18] << 0;
		f_width += file_buffer[19] << 8;
		f_width += file_buffer[20] << 16;
		f_width += file_buffer[21] << 24;
		/* height of image in pixels */
		f_height = 0;
		f_height += file_buffer[22] << 0;
		f_height += file_buffer[23] << 8;
		f_height += file_buffer[24] << 16;
		f_height += file_buffer[25] << 24;
		/* x pixel per meter value */
		f_xpixelpermeter = 0;
		f_xpixelpermeter += file_buffer[38] << 0;
		f_xpixelpermeter += file_buffer[39] << 8;
		f_xpixelpermeter += file_buffer[40] << 16;
		f_xpixelpermeter += file_buffer[41] << 24;
		/* y pixel per meter value */
		f_ypixelpermeter = 0;
		f_ypixelpermeter += file_buffer[42] << 0;
		f_ypixelpermeter += file_buffer[43] << 8;
		f_ypixelpermeter += file_buffer[44] << 16;
		f_ypixelpermeter += file_buffer[45] << 24;

		bw = f_width;
		bh = f_height;
		bitmap_width = f_width;
		bitmap_height = f_height;
		xdpi = f_xpixelpermeter;
		ydpi = f_ypixelpermeter;

		/* print info */ STRING_PRINTV("dimension is "); STRING_PRINTVD(bitmap_width);
		/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(bitmap_height); STRING_PRINTV2(" pixels\n");
		/* print info */ STRING_PRINTV("resolution is "); STRING_PRINTVD(xdpi);
		/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(ydpi); STRING_PRINTV2(" pixels per meter\n");
		/* print info */ STRING_PRINTV("colors have "); STRING_PRINTVD(f_bitcount);
		/* print info */ STRING_PRINTV2(" bit depth\n");

#ifndef __BMP_ONLY__
        if (f_bitcount == 8) { bitmap_format_jpg_clrspc_type = 8; }
        if (f_bitcount == 24){ bitmap_format_jpg_clrspc_type = 24; }
        if (f_bitcount == 8) { bitmap_format_png_clrspc_type = 0; }
        if (f_bitcount == 24){ bitmap_format_png_clrspc_type = 2; }
#endif

	/* allocate memory for the unpacked RGB colors */
	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
    	bitmap_buffer = malloc(bw * bh * 3 * sizeof(*bitmap_buffer));
    	if (bitmap_buffer == 0) return 1;

	    /* print info */ STRING_PRINTV("copying colors from bmp object to bitmap buffer\n");
            if (f_bitcount == 8){
       		addr_offset = f_width % 4;
               	if (addr_offset) addr_offset = 4 - addr_offset;
		#ifdef __OPENMP__
		#pragma omp parallel for private(x, y, addr, addr2) num_threads(max_threads)
		#endif
    	    	for (y=0; y<=bh-1; y++){
    	    		addr2 = y * bw * 1 + y * addr_offset;
    	    		for (x=0; x<=bw-1; x++){
    	    			addr = f_offs + addr2 + x * 1;
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 0] = file_buffer[addr + 0];
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 1] = file_buffer[addr + 0];
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 2] = file_buffer[addr + 0];
    	    		}
    	    	}
            }
            else{
       		addr_offset = (f_width) * 3 % 4;
               	if (addr_offset) addr_offset = 4 - addr_offset;
		#ifdef __OPENMP__
		#pragma omp parallel for private(x, y, addr, addr2) num_threads(max_threads)
		#endif
    	    	for (y=0; y<=bh-1; y++){
    	    		addr2 = y * bw * 3 + y * addr_offset;
    	    		for (x=0; x<=bw-1; x++){
    	    			addr = f_offs + addr2 + x * 3;
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 0] = file_buffer[addr + 2];
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 1] = file_buffer[addr + 1];
    	    			bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 2] = file_buffer[addr + 0];
    	    		}
    	    	}
            }

	    /* print info */ STRING_PRINTV("freeing bmp object\n");
            free(file_buffer);
	    return 0;
	}

	/* free memory */
	/* print info */ STRING_PRINTV("freeing bmp object\n");
	free(file_buffer);

	return 1;
}




/* -------------------------------------------------- */
/* ---------- BITMAP WRITE IN BMP FORMAT ------------ */
/* -------------------------------------------------- */
int BITMAP_WRITE_BMP(char *file_name)
{
	unsigned long f_offs;
	unsigned long addr, addr2;
	unsigned long addr_offset;
	unsigned long bw, bh;
	unsigned long f_xpixelpermeter;
	unsigned long f_ypixelpermeter;
	unsigned long temp;
	unsigned long file_length;
	unsigned char *file_buffer;
	unsigned long col;
	unsigned long i, x, y;
	FILE *fhandle;

	/* print info */ STRING_PRINTV("bmp initializations\n");

	bw = bitmap_width;
	bh = bitmap_height;

	f_xpixelpermeter = xdpi;
	f_ypixelpermeter = ydpi;

	/* allocate memory for unpacked RGB colors (3 plus bytes more in every width because of the BMP's 4 byte align adjust) */
	/* print info */ STRING_PRINTV("allocating memory for bmp object\n");
	if (bitmap_format_bmp_clrspc_type == 8) {
	    	addr_offset = bw % 4;
	        if (addr_offset) addr_offset = 4 - addr_offset;
		file_length = 54 + 4 * 256 + (bw + addr_offset) * bh;
	}
	else {
	    	addr_offset = (bw * 3) % 4;
       		if (addr_offset) addr_offset = 4 - addr_offset;
		file_length = 54 + (bw * 3 + addr_offset) * bh;
	}

	file_buffer = calloc(file_length, sizeof (*file_buffer));
	if (file_buffer == 0) return 1;

	/* write BMP indicator */
	file_buffer[0] = 0x42;
	file_buffer[1] = 0x4d;

	/* file length marker */
	file_buffer[2] = (file_length & 0x000000ff) >> 0;
	file_buffer[3] = (file_length & 0x0000ff00) >> 8;
	file_buffer[4] = (file_length & 0x00ff0000) >> 16;
	file_buffer[5] = (file_length & 0xff000000) >> 24;

	/* zero */
	file_buffer[6] = 0;
	file_buffer[7] = 0;
	file_buffer[8] = 0;
	file_buffer[9] = 0;

	if (bitmap_format_bmp_clrspc_type == 8) {
	    	/* bitmap offset (standard = 54) */
    		file_buffer[10] = 54;
	    	file_buffer[11] = 4; /* + 4 * 256 pieces of RGB gray colors = 0x0400 */
    		file_buffer[12] = 0;
	    	file_buffer[13] = 0;
    		/* number of bits per pixel */
	    	file_buffer[28] = 8;
    		file_buffer[29] = 0;
	}
	else {
	    	/* bitmap offset (standard = 54) */
    		file_buffer[10] = 54;
		file_buffer[11] = 0;
		file_buffer[12] = 0;
		file_buffer[13] = 0;
    		/* number of bits per pixel */
		file_buffer[28] = 24;
		file_buffer[29] = 0;
	}

	/* bitmap info header (standard = 40) */
	file_buffer[14] = 40;
	file_buffer[15] = 0;
	file_buffer[16] = 0;
	file_buffer[17] = 0;
	/* width of image in pixels */
	file_buffer[18] = (bw & 0x000000ff) >> 0;
	file_buffer[19] = (bw & 0x0000ff00) >> 8;
	file_buffer[20] = (bw & 0x00ff0000) >> 16;
	file_buffer[21] = (bw & 0xff000000) >> 24;
	/* height of image in pixels */
	file_buffer[22] = (bh & 0x000000ff) >> 0;
	file_buffer[23] = (bh & 0x0000ff00) >> 8;
	file_buffer[24] = (bh & 0x00ff0000) >> 16;
	file_buffer[25] = (bh & 0xff000000) >> 24;
	/* number of planes */
	file_buffer[26] = 1;
	file_buffer[27] = 0;
	/* compression (standard = 0) */
	file_buffer[30] = 0;
	file_buffer[31] = 0;
	file_buffer[32] = 0;
	file_buffer[33] = 0;
	/* size of bitmap image in bytes */
	if (bitmap_format_bmp_clrspc_type == 8) { temp = (bw + addr_offset) * bh; }
	else                                    { temp = (bw * 3 + addr_offset) * bh; }
	file_buffer[34] = (temp & 0x000000ff) >> 0;
	file_buffer[35] = (temp & 0x0000ff00) >> 8;
	file_buffer[36] = (temp & 0x00ff0000) >> 16;
	file_buffer[37] = (temp & 0xff000000) >> 24;
	/* x pixel per meter value */
	file_buffer[38] = (f_xpixelpermeter & 0x000000ff) >> 0;
	file_buffer[39] = (f_xpixelpermeter & 0x0000ff00) >> 8;
	file_buffer[40] = (f_xpixelpermeter & 0x00ff0000) >> 16;
	file_buffer[41] = (f_xpixelpermeter & 0xff000000) >> 24;
	/* y pixel per meter value */
	file_buffer[42] = (f_ypixelpermeter & 0x000000ff) >> 0;
	file_buffer[43] = (f_ypixelpermeter & 0x0000ff00) >> 8;
	file_buffer[44] = (f_ypixelpermeter & 0x00ff0000) >> 16;
	file_buffer[45] = (f_ypixelpermeter & 0xff000000) >> 24;

	for (i=46; i<54; i++) file_buffer[i] = 0;

	/* print info */ STRING_PRINTV("copying colors from bitmap buffer to bmp object\n");
	if (bitmap_format_bmp_clrspc_type == 8) {
		#ifdef __OPENMP__
		#pragma omp parallel for num_threads(max_threads)
		#endif
		for (i=0; i<256; i++){
			file_buffer[54 + i*4 + 0] = i;
			file_buffer[54 + i*4 + 1] = i;
			file_buffer[54 + i*4 + 2] = i;
			file_buffer[54 + i*4 + 3] = 0;
        	}

	    	f_offs = 54 + 4 * 256;
    		addr = 0;
		#ifdef __OPENMP__
		#pragma omp parallel for private(x, y, addr, addr2, col) num_threads(max_threads)
		#endif
    		for (y=0; y<=bh-1; y++){
    			addr2 = y * (bw + addr_offset);
	    		for (x=0; x<=bw-1; x++){
    				addr = f_offs + addr2 + x;
				col  = bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 0];
				col += bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 1];
				col += bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 2];
    				file_buffer[addr + 0] = col / 3;
    			}
    		}
	}

	else {
	    	f_offs = 54;
    		addr = 0;
		#ifdef __OPENMP__
		#pragma omp parallel for private(x, y, addr, addr2) num_threads(max_threads)
		#endif
    		for (y=0; y<=bh-1; y++){
    			addr2 = y * (bw * 3 + addr_offset);
	    		for (x=0; x<=bw-1; x++){
    				addr = f_offs + addr2 + x * 3;
    				file_buffer[addr + 0] = bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 2];
    				file_buffer[addr + 1] = bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 1];
	    			file_buffer[addr + 2] = bitmap_buffer[x * 3 + bw * 3 * (bh-1-y) + 0];
    			}
	    	}
	}


	/* open file for writing */
	/* print info */ STRING_PRINTV("opening file for writing\n");
	fhandle = fopen(file_name, "wb");
	if (fhandle == 0) return 1;
	/* write and flush file from memory */
	/* print info */ STRING_PRINTV("writing image file\n");
	if (fwrite(file_buffer, 1, file_length, fhandle) == 0) return 1;
	/* close file */
	fclose(fhandle);

	/* print info */ STRING_PRINTV("freeing bmp object\n");
	free(file_buffer);

	return 0;
}


#ifndef __BMP_ONLY__

/* ----------------------------------------------- */
/* --------- BITMAP READ IN PNG FORMAT ----------- */
/* ----------------------------------------------- */
/* source code examples for PNG read taken from: */
/* http://www.libpng.org/pub/png/book/chapter13.html */
int BITMAP_READ_PNG(char *file_name)
{
	FILE *fhandle;
	unsigned long x, y;
	long offs1, offs2;
	unsigned long   rowbytes;
	unsigned char **row_pointers;
	unsigned long	i;
	int alpha_flag;
	png_uint_32 res_x, res_y;
	int unit_type;
	png_uint_32 png_width, png_height;
	int bit_depth;
	png_structp png_ptr;
	png_infop info_ptr;
	unsigned char sig[8];


        /* open file for reading */
	/* print info */ STRING_PRINTV("opening file for reading\n");
        fhandle = fopen(file_name, "rb");
        if (fhandle == 0) return 1;

	/* check PNG signature */
	if (fread(sig, 1, 8, fhandle) != 8) return 1;

/*	libpng version change to 1.4.0 */
/*	if (!png_check_sig(sig, 8)) return 1; */
	if (png_sig_cmp(sig, 0, 8)) return 1;


	/* print info */ STRING_PRINTV("png initializations\n");
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	/* out of memory */
	if (!png_ptr) return 1;

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
        	png_destroy_read_struct(&png_ptr, NULL, NULL);
        	return 1;
	}

/*	libpng version change to 1.4.0 */
/*	if (setjmp(png_ptr->jmpbuf)) { */
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 1;
	}

	png_init_io(png_ptr, fhandle);
	png_set_sig_bytes(png_ptr, 8);
	png_read_info(png_ptr, info_ptr);

	/* get image header info */
	png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height, &bit_depth, \
		&bitmap_format_png_clrspc_type, &bitmap_format_png_interlace_type, \
		&bitmap_format_png_compression_type, &bitmap_format_png_filter_type);

	bitmap_width = png_width;
	bitmap_height = png_height;

	/* get image resolution */
	png_get_pHYs(png_ptr, info_ptr, &res_x, &res_y, &unit_type);

	if (unit_type == 0) { res_x = 0; res_y = 0; }

	xdpi = res_x;
	ydpi = res_y;
	udpi = unit_type;

	/* print info */ STRING_PRINTV("dimension is "); STRING_PRINTVD(bitmap_width);
	/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(bitmap_height); STRING_PRINTV2(" pixels\n");
	/* print info */ if (udpi == 1)	{
	/* print info */ STRING_PRINTV("resolution is "); STRING_PRINTVD(xdpi);
	/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(ydpi);
	/* print info */ STRING_PRINTV2(" pixels per meter\n"); }
	/* print info */ else {
	/* print info */ STRING_PRINTV2("resolution is of unknown type\n"); }
	/* print info */ STRING_PRINTV("colors have "); STRING_PRINTVD(bit_depth);
	/* print info */ STRING_PRINTV2(" bit depth\n");

	/* bitmap depth must be 8 bit color, RGB or Grayscale */
	if (bit_depth != 8) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 1;
	}

	/* check color types */
	/* 0 - Gray */
	/* 2 - RGB */
	/* 4 - Gray + Alpha */
	/* 6 - RGB + Alpha */
        if (bitmap_format_png_clrspc_type != 0 &&
	    bitmap_format_png_clrspc_type != 2 &&
	    bitmap_format_png_clrspc_type != 4 &&
	    bitmap_format_png_clrspc_type != 6) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 1;
	}

	/* set color types */
        if (bitmap_format_png_clrspc_type == 0 || bitmap_format_png_clrspc_type == 4) {
		bitmap_format_jpg_clrspc_type = 8;
		bitmap_format_bmp_clrspc_type = 8; }
        if (bitmap_format_png_clrspc_type == 2 || bitmap_format_png_clrspc_type == 6) {
		bitmap_format_jpg_clrspc_type = 24;
		bitmap_format_bmp_clrspc_type = 24; }

	/* if there is Alpha channel, then allocate bigger memory for it */
	/* RGB + Alpha needs 5 times of the pixels because the bytes need to be rearranged */
	alpha_flag = 0;
        if (bitmap_format_png_clrspc_type == 4) { alpha_flag = 1; }
        if (bitmap_format_png_clrspc_type == 6) { alpha_flag = 2; }

	if (alpha_flag){ /* print info */ STRING_PRINTV("alpha channel exists\n"); }
	else           { /* print info */ STRING_PRINTV("no alpha channel\n");     }

	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
	bitmap_buffer = malloc(bitmap_width * bitmap_height * (3 + alpha_flag) * sizeof (*bitmap_buffer));
	if (bitmap_buffer == 0) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		return 1;
	}


	png_read_update_info(png_ptr, info_ptr);
	rowbytes = png_get_rowbytes(png_ptr, info_ptr);

	row_pointers = malloc(bitmap_height * sizeof (long));
	if (row_pointers == 0) { return 1; }

	/* print info */ STRING_PRINTV("creating row pointers for bitmap buffer\n");
	#ifdef __OPENMP__
	#pragma omp parallel for num_threads(max_threads)
	#endif
	for (i = 0; i < bitmap_height; i++) {
		row_pointers[i] = (unsigned char*) (bitmap_buffer + i*rowbytes);
	}

	/* print info */ STRING_PRINTV("reading image file and copy colors to bitmap buffer\n");
	png_read_image(png_ptr, row_pointers);
	png_read_end(png_ptr, NULL);

	/* destroy png structure and free memory */
	/* print info */ STRING_PRINTV("freeing png object\n");
	if (png_ptr && info_ptr) {
		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		png_ptr = NULL;
		info_ptr = NULL;
	}



	/* if it's Grayscale, then we pack the gray bytes to same RGB bytes */
	/* because the aaRGB function expects RGB bytes */
	/* if it's RGB picture with Alpha channel, then we pack the */
	/* RGB bytes to the begining of the allocated memory */
	/* and the Alpha bytes to the end */
	/* so it will be compatible in case of JPG output too */
	/* only the alpha channel will be lost */
	/* cause JPG does not support alpha */


	/* print info */ STRING_PRINTV("converting colors to different format\n");
	/* Gray */
        if (bitmap_format_png_clrspc_type == 0){
		/* copy Gray bytes to be RGB bytes */
		offs1 = bitmap_width * bitmap_height * 1 - 1;
		offs2 = bitmap_width * bitmap_height * 3 - 3;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2 + 0] = bitmap_buffer[offs1];
				bitmap_buffer[offs2 + 1] = bitmap_buffer[offs1];
				bitmap_buffer[offs2 + 2] = bitmap_buffer[offs1];
				offs1--;
				offs2 = offs2 - 3;
			}
		}
	}
	/* RGB */
        if (bitmap_format_png_clrspc_type == 2){
		/* do nothing */
	}
	/* Gray + Alpha */
        if (bitmap_format_png_clrspc_type == 4){
		/* move Alpha bytes to the end */
		offs1 = bitmap_width * bitmap_height * 2 - 1;
		offs2 = bitmap_width * bitmap_height * 4 - 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1 = offs1 - 2;
				offs2--;
			}
		}
		/* copy Gray bytes to be RGB bytes */
		offs1 = bitmap_width * bitmap_height * 2 - 2;
		offs2 = offs2 - 2;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2 + 0] = bitmap_buffer[offs1];
				bitmap_buffer[offs2 + 1] = bitmap_buffer[offs1];
				bitmap_buffer[offs2 + 2] = bitmap_buffer[offs1];
				offs1 = offs1 - 2;
				offs2 = offs2 - 3;
			}
		}
	}
	/* RGB + Alpha */
        if (bitmap_format_png_clrspc_type == 6){
		/* move Alpha bytes to the 5. endpart */
		offs1 = bitmap_width * bitmap_height * 4 - 1;
		offs2 = bitmap_width * bitmap_height * 5 - 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1 = offs1 - 4;
				offs2--;
			}
		}
		/* transfer RGBA bytes into RGB bytes */
		offs1 = 0;
		offs2 = 0;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2 + 0] = bitmap_buffer[offs1 + 0];
				bitmap_buffer[offs2 + 1] = bitmap_buffer[offs1 + 1];
				bitmap_buffer[offs2 + 2] = bitmap_buffer[offs1 + 2];
				offs1 = offs1 + 4;
				offs2 = offs2 + 3;
			}
		}
		/* move Alpha bytes back to the 4. part */
		offs1 = bitmap_width * bitmap_height * 5 - 1;
		offs2 = bitmap_width * bitmap_height * 4 - 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1--;
				offs2--;
			}
		}
	}

	free (row_pointers);
	return 0;
}




/* ------------------------------------------------ */
/* --------- BITMAP WRITE IN PNG FORMAT ----------- */
/* ------------------------------------------------ */
/* source code examples for PNG write taken from: */
/* http://www.libpng.org/pub/png/book/chapter15.html */
int BITMAP_WRITE_PNG(char *file_name)
{
        unsigned long   rowbytes;
	unsigned long	i;
	unsigned char **row_pointers;
        FILE *fhandle;
	png_structp png_ptr;
	png_infop info_ptr;

	unsigned long x, y;
	long offs1, offs2;
	long col;

	/* repack the Gray and RGB bytes with the Alpha bytes */
	/* print info */ STRING_PRINTV("converting colors to different format\n");

	/* Gray */
        if (bitmap_format_png_clrspc_type == 0){
		/* copy RGB bytes back as Gray bytes */
		offs1 = 0;
		offs2 = 0;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				col = 0;
				col += bitmap_buffer[offs2 + 0];
				col += bitmap_buffer[offs2 + 1];
				col += bitmap_buffer[offs2 + 2];
				col /= 3;
				bitmap_buffer[offs1] = col;
				offs1++;
				offs2 = offs2 + 3;
			}
		}
	}
	/* RGB */
        if (bitmap_format_png_clrspc_type == 2){
		/* do nothing */
	}
	/* Gray + Alpha */
        if (bitmap_format_png_clrspc_type == 4){
		/* copy RGB bytes back as Gray bytes */
		offs1 = 0;
		offs2 = 0;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				col = 0;
				col += bitmap_buffer[offs2 + 0];
				col += bitmap_buffer[offs2 + 1];
				col += bitmap_buffer[offs2 + 2];
				col /= 3;
				bitmap_buffer[offs1] = col;
				offs1 = offs1 + 2;
				offs2 = offs2 + 3;
			}
		}
		/* copy Alpha bytes back */
		offs1 = bitmap_width * bitmap_height * 3 + 0;
		offs2 = 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1++;
				offs2 = offs2 + 2;
			}
		}
	}
	/* RGB + Alpha */
        if (bitmap_format_png_clrspc_type == 6){
		/* move Alpha bytes back from the 4. to the 5. endpart */
		offs1 = bitmap_width * bitmap_height * 4 - 1;
		offs2 = bitmap_width * bitmap_height * 5 - 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1--;
				offs2--;
			}
		}
		/* transfer RGB bytes into RGBA bytes */
		offs1 = bitmap_width * bitmap_height * 3 - 3;
		offs2 = bitmap_width * bitmap_height * 4 - 4;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2 + 2] = bitmap_buffer[offs1 + 2];
				bitmap_buffer[offs2 + 1] = bitmap_buffer[offs1 + 1];
				bitmap_buffer[offs2 + 0] = bitmap_buffer[offs1 + 0];
				offs1 = offs1 - 3;
				offs2 = offs2 - 4;
			}
		}
		/* move Alpha bytes back from the 5. endpart to RGBA */
		offs1 = bitmap_width * bitmap_height * 5 - 1;
		offs2 = bitmap_width * bitmap_height * 4 - 1;
		for (y=0; y < bitmap_height; y++){
			for (x=0; x < bitmap_width; x++){
				bitmap_buffer[offs2] = bitmap_buffer[offs1];
				offs1--;
				offs2 = offs2 - 4;
			}
		}
	}


	/* print info */ STRING_PRINTV("png initializations\n");
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	/* out of memory */
	if (!png_ptr) return 1;

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr, NULL);
		return 1;
	}

/*	libpng version change to 1.4.0 */
/*	if (setjmp(png_ptr->jmpbuf)) { */
	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_write_struct(&png_ptr, &info_ptr);
		return 1;
	}

        /* open file for writing */
	/* print info */ STRING_PRINTV("opening file for writing\n");
        fhandle = fopen(file_name, "wb");
        if (fhandle == 0) return 1;

	png_init_io(png_ptr, fhandle);
/*	png_set_compression_level(png_ptr, Z_BEST_COMPRESSION);*/
/*	set the compression level manually */
	png_set_compression_level(png_ptr, 5);

	/* set image header */
	png_set_IHDR(png_ptr, info_ptr, bitmap_width, bitmap_height,
/*		8, bitmap_format_png_clrspc_type, bitmap_format_png_interlace_type, */
/*		bitmap_format_png_compression_type, bitmap_format_png_filter_type); */
		8, bitmap_format_png_clrspc_type, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

        /* set image resolution */
        png_set_pHYs(png_ptr, info_ptr, xdpi, ydpi, udpi);


	png_write_info(png_ptr, info_ptr);
/*	png_set_packing(png_ptr); */


        png_read_update_info(png_ptr, info_ptr);
        rowbytes = png_get_rowbytes(png_ptr, info_ptr);

	row_pointers = malloc(bitmap_height * sizeof (long));
	if (row_pointers == 0) { return 1; }

	/* print info */ STRING_PRINTV("creating row pointers for bitmap buffer\n");
	#ifdef __OPENMP__
	#pragma omp parallel for num_threads(max_threads)
	#endif
        for (i = 0; i < bitmap_height; i++) {
                row_pointers[i] = (unsigned char*) (bitmap_buffer + i*rowbytes);
        }

	/* print info */ STRING_PRINTV("copying colors from bitmap buffer and write image file\n");
	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, NULL);

	/* print info */ STRING_PRINTV("freeing png object\n");
	if (png_ptr && info_ptr)
		png_destroy_write_struct(&png_ptr, &info_ptr);

	free (row_pointers);
	return 0;
}




/* ------------------------------------------------ */
/* --------- BITMAP READ IN JPEG FORMAT ----------- */
/* ------------------------------------------------ */
/* source code examples for reading JPEG taken from libjpeg's example.c */
int BITMAP_READ_JPEG(char *file_name)
{
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *infile;
	JSAMPARRAY buffer;
	int row_stride;

	EXIF_GET(file_name);

	/* print info */ STRING_PRINTV("jpeg initializations\n");

	/* print info */ STRING_PRINTV("opening file for reading\n");
	if ((infile = fopen(file_name, "rb")) == NULL) { return 1; }

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, infile);
	(void) jpeg_read_header(&cinfo, TRUE);

	 /* print info */ STRING_PRINTV("decoding jpeg image\n");
	(void) jpeg_start_decompress(&cinfo);

	row_stride = cinfo.output_width * cinfo.output_components;
	buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
	bitmap_width = cinfo.output_width;
	bitmap_height = cinfo.output_height;
	bitmap_buffer = malloc(bitmap_width * bitmap_height * 3 * sizeof (*bitmap_buffer));
	if (bitmap_buffer == 0) return 1;

	xdpi = cinfo.X_density;
	ydpi = cinfo.Y_density;
	udpi = 0;
	/* unit = 0, no units */
	/* unit = 1, unit is pixels per inch */
	/* unit = 2, unit is pixels per centimetre */
	/* convert the unit to meter, cause BMP and PNG stores it as meter too with unit being 1 */
	if (cinfo.density_unit == 1){ xdpi = xdpi * 1000 / 25.4; ydpi = ydpi * 1000 / 25.4; udpi = 1; }
	if (cinfo.density_unit == 2){ xdpi = xdpi * 100; ydpi = ydpi * 100; udpi = 1; }

	/* determine colorspace type */
	if (cinfo.output_components == 1) { bitmap_format_jpg_clrspc_type = 8; }
	if (cinfo.output_components == 3) { bitmap_format_jpg_clrspc_type = 24; }
        if (bitmap_format_jpg_clrspc_type == 8)  { bitmap_format_png_clrspc_type = 0;  }
        if (bitmap_format_jpg_clrspc_type == 24) { bitmap_format_png_clrspc_type = 2;  }
        if (bitmap_format_jpg_clrspc_type == 8)  { bitmap_format_bmp_clrspc_type = 8;  }
        if (bitmap_format_jpg_clrspc_type == 24) { bitmap_format_bmp_clrspc_type = 24; }

	/* print info */ STRING_PRINTV("dimension is "); STRING_PRINTVD(bitmap_width);
	/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(bitmap_height); STRING_PRINTV2(" pixels\n");
	/* print info */ STRING_PRINTV("resolution is "); STRING_PRINTVD(xdpi);
	/* print info */ STRING_PRINTV2(" x "); STRING_PRINTVD(ydpi);
	/* print info */ if (udpi == 1)	{ STRING_PRINTV2(" pixels per meter\n"); }
	/* print info */ else		{ STRING_PRINTV2(" with unknown type\n"); }
	/* print info */ STRING_PRINTV("colors have "); STRING_PRINTVD(bitmap_format_jpg_clrspc_type);
	/* print info */ STRING_PRINTV2(" bit depth\n");

	/* print info */ STRING_PRINTV("copying colors from jpeg object to bitmap buffer\n");
	while (cinfo.output_scanline < cinfo.output_height) {
		unsigned int i;

		(void) jpeg_read_scanlines(&cinfo, buffer, 1);
/*		put_scanline_someplace(buffer[0], row_stride); */

		if (cinfo.output_components == 3){
			#ifdef __OPENMP__
			#pragma omp parallel for num_threads(max_threads)
			#endif
			for (i=0; i<cinfo.output_width; i++){
				bitmap_buffer[0 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 0 + i*3);
				bitmap_buffer[1 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 1 + i*3);
				bitmap_buffer[2 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 2 + i*3);
			}
		}
		else{
			#ifdef __OPENMP__
			#pragma omp parallel for num_threads(max_threads)
			#endif
			for (i=0; i<cinfo.output_width; i++){
				bitmap_buffer[0 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 0 + i);
				bitmap_buffer[1 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 0 + i);
				bitmap_buffer[2 + i*3 + (cinfo.output_scanline-1) * cinfo.output_width*3] = *(*buffer + 0 + i);
			}
		}
	}

	/* print info */ STRING_PRINTV("freeing jpeg object\n");
	(void) jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(infile);

	return 0;
}




/* ------------------------------------------------- */
/* --------- BITMAP WRITE IN JPEG FORMAT ----------- */
/* ------------------------------------------------- */
/* source code examples for writing JPEG taken from libjpeg's example.c */
int BITMAP_WRITE_JPEG(char *file_name)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *outfile;
	JSAMPROW row_pointer[1];
	int row_stride;

	/* write back exif info if there exists any */
	EXIF_PUT(file_name);

	/* print info */ STRING_PRINTV("jpeg initializations\n");

	cinfo.err = jpeg_std_error(&jerr);
        jpeg_create_compress(&cinfo);

	/* print info */ STRING_PRINTV("opening file for writing\n");
	/* if there was exif info, then open file in append mode */
	/* otherwise open it in write mode */
	if (exif_flag){
		if ((outfile = fopen(file_name, "a+b")) == NULL) { return 1; } }
	else{
		if ((outfile = fopen(file_name, "wb")) == NULL) { return 1; } }

	jpeg_stdio_dest(&cinfo, outfile);

	cinfo.image_width = bitmap_width;
	cinfo.image_height = bitmap_height;
	if (bitmap_format_jpg_clrspc_type == 8){
		cinfo.input_components = 1;
		cinfo.in_color_space = JCS_GRAYSCALE;
	}
	if (bitmap_format_jpg_clrspc_type == 24){
		cinfo.input_components = 3;
		cinfo.in_color_space = JCS_RGB;
	}

	jpeg_set_defaults(&cinfo);

	/* set image resolution */
	if (udpi == 0){
		/* unknown unit */
		cinfo.X_density = xdpi;
		cinfo.Y_density = ydpi;
		cinfo.density_unit = 0;
	}
	if (udpi == 1){
		/* convert from pixel per meter to inch with rounding up to the nearest decimal */
		int xnew = (int)(xdpi * 25.4 / 100);
		int ynew = (int)(ydpi * 25.4 / 100);
		if ((xnew % 10) > 5){ xnew += 5; }
		if ((ynew % 10) > 5){ ynew += 5; }
		xnew /= 10;
		ynew /= 10;
		cinfo.X_density = xnew;
		cinfo.Y_density = ynew;
		cinfo.density_unit = 1;
	}

	/* set default quality to 95 if unset */
	if (opt_quality < 0)  { opt_quality = 0;   }
	if (opt_quality > 100){ opt_quality = 100; }
	if (!(opt_quality))   { opt_quality = 95;  }
	jpeg_set_quality(&cinfo, opt_quality, TRUE);

	/* print info */ STRING_PRINTV("encoding jpeg image\n");
	jpeg_start_compress(&cinfo, TRUE);

	row_stride = bitmap_width * cinfo.input_components;

	/* if image is 1 component Gray, then copy RGB bytes back as Gray bytes */
	if (cinfo.input_components == 1){
                int offs1 = 0;
                int offs2 = 0;
		unsigned int col, x, y;
		/* this cycle cannot be paralleled because of writing back to the same memory */
                for (y=0; y < bitmap_height; y++){
                        for (x=0; x < bitmap_width; x++){
                                col = 0;
                                col += bitmap_buffer[offs2 + 0];
                                col += bitmap_buffer[offs2 + 1];
                                col += bitmap_buffer[offs2 + 2];
                                col /= 3;
                                bitmap_buffer[offs1] = col;
                                offs1++;
                                offs2 += 3;
                        }
                }
	}

	/* print info */ STRING_PRINTV("copying colors from bitmap buffer to jpeg object\n");
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &bitmap_buffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	/* print info */ STRING_PRINTV("freeing jpeg object\n");
	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);


	/* correct the end of exif data in file and free exif buffer if there was any */
	EXIF_CORRECT(file_name);
	EXIF_CLEAR();

	return 0;
}

#endif



/* --------------------------------------- */
/* ----------- BITMAP LOAD --------------- */
/* --------------------------------------- */
/* load image and unpack it into memory (memory will be allocated) */
int BITMAP_LOAD(char *file_name)
{
    FILE *fhandle;

    /* clear exif flag at new file */
    exif_flag = 0;

	/* does file exist? */
	fhandle = fopen(file_name, "rb");
	if (fhandle == 0) {
		return 2;
	}
	fclose(fhandle);

    /* print info */ STRING_PRINTV("file exists\n");
    /* if it's a BMP format, then load BMP with custom procedure */
    if (GET_FILE_FORMAT(file_name) == 0){
        if (BITMAP_READ_BMP(file_name)) return 1;
    }
    else{

#ifdef __BMP_ONLY__
    STRING_PRINT("\n"); STRING_PRINTE("BMP_ONLY version. Only BMP format is supported here.n");
    return 1;
#endif
#ifndef __BMP_ONLY__


	/* if it's a JPEG, then save exif info and load image with libjpeg */
        if (GET_FILE_FORMAT(file_name) == 1){
            if (BITMAP_READ_JPEG(file_name)) return 1;
            return 0;
        }

	/* if it's a PNG format, then load that one */
        if (GET_FILE_FORMAT(file_name) == 2){
            if (BITMAP_READ_PNG(file_name)) return 1;
            return 0;
        }

        return 0;


#endif
    }

	return 0;
}




/* --------------------------------------- */
/* ----------- BITMAP SAVE --------------- */
/* --------------------------------------- */
/* pack the image file and save it from memory (memory will be freed) */
int BITMAP_SAVE(char *file_name)
{
    /* if the format is not JPEG, then clear exif info (memory is freed) */
    if (GET_FILE_FORMAT(file_name) != 1){ EXIF_CLEAR(); }

    /* if it's a BMP format, then save BMP with custom procedure */
    if (GET_FILE_FORMAT(file_name) == 0){
        if (BITMAP_WRITE_BMP(file_name)) { free(bitmap_buffer); return 1; }

	/* free memory on success */
	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(bitmap_buffer);
	return 0;
    }
    else{

#ifdef __BMP_ONLY__
    STRING_PRINT("\n"); STRING_PRINTE("BMP_ONLY version. Only BMP format is supported here.n");
    /* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
    free(bitmap_buffer);
    return 1;
#endif
#ifndef __BMP_ONLY__


    /* if it's a JPEG, then save it */
    if (GET_FILE_FORMAT(file_name) == 1){
        if (BITMAP_WRITE_JPEG(file_name)) {
	    /* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	    free(bitmap_buffer); return 1; }

	/* free memory on success */
	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(bitmap_buffer);
	return 0;
    }

    /* if it's a PNG, then save PNG */
    if (GET_FILE_FORMAT(file_name) == 2){
        if (BITMAP_WRITE_PNG(file_name)) {
   	    /* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	    free(bitmap_buffer); return 1; }

	/* free memory on success */
	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(bitmap_buffer);
	return 0;
    }


    /* free memory on success */
    /* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
    free(bitmap_buffer);
    return 0;

#endif
    }

}

