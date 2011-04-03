/*-----------------------------------------------------------------------
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


/* #define __BMP_ONLY__ */
/* #define __OPENMP__ */



/* --------------------------------------------------- */
/* ----------- Auto Adjust Photo --------------------- */
/* ----------- András Horváth (C) 2006-2011 ---------- */
/* ----------- Hungary, http://log69.com ------------- */
/* --------------------------------------------------- */

/*
aaphoto Changelog:
--------------------
2011/01/26 - aaphoto v0.41 - new aaRGB v0.64 version update (see aaRGB changelog)
                           - add -- switch to mark the end of option list for posix compatibility
                           - fix some warnings given by -Wextra compile option
2010/12/18 - aaphoto v0.40 - new aaRGB v0.63 version update
                           - fix some warning messages during build
                           - some changes in documentation
                           - error message when --verbose switch is used with no other ones to avoid misunderstanding
                             that other switches are still needed (aaphoto -V file)
2010/09/14 - aaphoto v0.39 - new aaRGB v0.62 version update
                           - bugfix: an ugly misconception in my paralleled code caused weird behavior when using more threads
                           - bugfix: BMP image writer function didn't zero out the BMP align bytes in file buffer
                             therefore the same output differed in some cases
                           - bugfix: PNG image reader function stored the unreliable values of pixel resolution
                             when the type was unknown
                           - rewrite the code to suffice the ISO C90 ANSI standard C form (GCC -pedantic option)
                           - new windows platform patch of libjasper for static build to replace mkstemp() function call
                             because it is not yet supported under MinGW
                           - some changes in documentation
                           - new default LDFLAG -lgomp for OpenMP
                           - update all build scripts for Debian 6 platform
2010/07/18 - aaphoto v0.38 - add verbose message showing presence of alpha channel in PNG images
                           - bugfix: fix compile with only BMP support
                           - fix some warning messages during build
                           - remove unnecessary __WIN32__ macro
2010/05/10 - aaphoto v0.37 - add OpenMP support for multi processing, all possible time intensive codes paralleled
                             __OPENMP__ directive and -fopenmp option are needed at compile time (supported since GCC v4.2)
                             (not available on windows platform yet)
                           - new -t, --threads switch added to manually set the number of working threads
                           - new aaRGB v0.61 version update
                           - bugfix: exif info was left out after last version's change in JPEG handling
                           - improve exif handling and verbose messages
                           - refine additional verbose messages
                           - remove pgx file type support because it is outdated
                           - some changes in documentation
                           - some minor code cleanup and bugfixes
2010/03/16 - aaphoto v0.36 - bugfix: fix for tmpfile() patches of libjasper and libjpeg (windows platform only)
                             when running more than one instances of aaphoto at a time
                             they all used the same temporary files and therefore the images became corrupt
                           - bugfix: it doesn't ask for administrative privileges anymore while running in an admin account
                             under vista and windows 7 (windows platform only)
                           - bugfix: the --rotate180 switch didn't turn the middle line in images with odd heights
                           - rewrite JPEG format handling entirely to be able to handle extra parameters in this format
                             separately, now libjpeg is used directly instead of libjasper for reading / writing JPEG images
                             so libjpeg is a new dependency from now, formerly only libjasper was depended on it
                           - restore original DPI values of images in BMP, JPEG and PNG formats during conversion
                           - refine program messages
                           - print the time elapsed in seconds since program start in verbose mode if not zero
                           - print extra infos of bitmap dimension, resolution and color depth in verbose mode
2010/02/25 - aaphoto v0.35 - bugfix: possible buffer overflows fixed
2010/02/19 - aaphoto v0.34 - __UNIX__ macro removed from Makefile, not needed anymore
                           - bugfix: static binary update with patches of libjasper and libjpeg
                             1) sleep() function missing on mingw32 platform
                             2) bad implementation of tmpfile() function on windows platform
                                it tries to create temporary files in the root of current directory
                                instead of the system temporary path so that causes failure for unprivileged users
                                who don't have permissions to write there
                           - update: changes in new version of libpng 1.4.0, aaio.c updated as necessary
                                png_check_sig() function replaced with png_sig_cmp()
                                setjmp(png_ptr->jmpbuf) has been deprecated, changed to setjmp(png_jmpbuf(png_ptr))
                                see more at http://www.libpng.org/pub/png/src/libpng-1.2.x-to-1.4.x-summary.txt
2010/01/10 - aaphoto v0.33 - some changes in documentation
                           - bugfix: unfreed space caused memory leak
                           - bugfix: uninitialized variable caused --resize to misbehave
                           - fix: change of return values in procedures to reflect standard exit codes
                             now it has a sense to run something like "aaphoto image.jpg && echo OK"
                             formerly return codes meant opposite
                           - fix a warning message during compile time, an include was missing
                           - boundary check of fixed size arrays added for safety reasons
                           - the --speed switch removed, it made the code less platform independent and was fussy anyway
                           - error messages printed to stderr instead of stdout from now
                           - more verbose error messages on failure of image load
2009/10/18 - aaphoto v0.32 - new aaRGB v0.60 version update
                           - new --noexif switch added to save new image without exif info
                           - new --bmp switch added for BMP format output
                           - new -o, --output switch added for alternate directory output
2009/08/23 - aaphoto v0.31 - bugfix: __BMP_ONLY__ directive is fixed in source code
                           - bugfix: writing of BMP images could result in corrupt BMP structure
                           - code cleanup in BMP write function
                           - parameters of switches also work with spaces between them
                           - new aaRGB v0.59 version update
2009/02/22 - aaphoto v0.30 - implement PNG format (RGB and Gray images read / write with alpha channel support)
                           - bugfix: reading corrupt exif info in JPEG files could get into an infinite loop
                           - bugfix: length of exif info was determined wrongly
                           - rework of the parameters and switches parsing part
                           - lots of code cleanup
                           - most of the comments in code translated to english
                           - print messages get flushed out with fflush now during process
                           - the --info switch removed
                           - the -o switch removed for safety reasons, --overwrite still available
                           - the -h switch added for unix compatibility
                           - the -j1 and -j2 switch removed, --jpg and --jp2 still available
                           - the --png switch is now new for PNG output
                           - the -s switch changed from --speed to --silent for compatibility
                           - the --mute switch changed to -s, --silent and --quiet for unix / posix compatibility
                           - the -V, --verbose switch is now new for more detailed output during image process
                           - the --test switch now turns the --autoadjust switch on by default
			   - the -d, --description and -l, --license switches removed
			   - thanks to Rezső Páder (rezso.net) for the suggestions for the option switches
                           - new aaRGB v0.58 version update
2008/02/02 - aaphoto v0.29 - bugfix: color space variable was not defined during the load of BMP format
                           - bugfix: BMP format handling fixed for JPEG conversion
                           - grayscale images can also be used as an input
                           - increase file name buffer for processing files in folders
                           - Exif meta data information is now restored during conversion in JPEG images
2007/08/11 - aaphoto v0.28 - new aaRGB v0.57 version update
                           - bugfix: remove extra slashes from the end of folders
2007/07/04 - aaphoto v0.27 - new aaRGB v0.56 version update v0.56 with "Apply only on selection" function
2007/05/26 - aaphoto v0.26 - bugfix: Win32 version crashed during JPEG-2000 conversion
2007/05/19 - aaphoto v0.25 - expand functions: rotate 90, 180, 270, flip x, flip y
2007/05/01 - aaphoto v0.24 - improve timing values
                           - input parameter can be folders beside files too
                           - bitmap info parameter now works with bmp too
                           - simplify parameter input: no --autoadjust parameter needed from now, default is on
2007/04/03 - aaphoto v0.23 - new aaRGB v0.55 version update
2007/04/01 - aaphoto v0.22 - new aaRGB v0.54 version update
2007/03/30 - aaphoto v0.21 - create bmp_only macro in source code for other platforms
                           - bug fixes
2007/03/29 - aaphoto v0.20 - new aaRGB v0.53 version update
2007/02/25 - aaphoto v0.19 - extra information output within image for testing purposes
2007/02/22 - aaphoto v0.18 - custom code for BMP input/output
2007/01/24 - aaphoto v0.17 - implement JasPer encoder for further image formats
2007/01/04 - aaphoto v0.16 - stable working command-line version for linux environment with BMP format support

aaphoto end of Changelog.

*/




/* global constants and variables */

#define max_char 1024
#define max_file_name_buffer 4096 * 1024

int   file_counter;
char *file_name;

unsigned char *file_name_buffer;
long  file_name_counter;
long  file_name_buffer_pointer;

unsigned char *bitmap_buffer;
unsigned long  bitmap_width;
unsigned long  bitmap_height;

double  xdpi;
double  ydpi;
int     udpi;

int   bitmap_format_bmp_clrspc_type;
int   bitmap_format_jpg_clrspc_type;
int   bitmap_format_jpg_file_type;
int   bitmap_format_png_clrspc_type;
int   bitmap_format_png_interlace_type;
int   bitmap_format_png_compression_type;
int   bitmap_format_png_filter_type;


char *exif_buffer;
long  exif_buffer_length;
long  exif_file_length;
int   exif_flag;

int opt_help;
int opt_version;
int opt_autoadjust;
int opt_overwrite;
int opt_jpg;
int opt_jp2;
int opt_png;
int opt_bmp;
int opt_resize;
int opt_resize_percent;
int opt_rotate90;
int opt_rotate180;
int opt_rotate270;
int opt_flipx;
int opt_flipy;
int opt_output;
char opt_output_path [max_char];
int opt_quality;
int opt_threads;
int opt_verbose;
int opt_recursive;
int opt_quiet;
int opt_test;
int opt_noexif;
int opt_doubledash;

int char_temp_x;
int char_temp_y;

int mytime;


/* OpenMP */
int max_cpus;
int max_threads;
int max_num_threads;
int opt_openmp;
#ifdef __OPENMP__
#include <omp.h>
#endif


/* file name separation character */
char slsh='/';


#ifndef __BMP_ONLY__
#include <jasper/jasper.h>
#include <png.h>
#include <jpeglib.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <time.h>
/*#include <stddir.h> */

#include "aaio.c"
#include "aargb.c"
#include "aaresize.c"




void PRINT_VERSION(void){
STRING_PRINT("Auto Adjust Photo\n");
STRING_PRINT("Copyright (C) 2006-2011 Andras Horvath\n");
STRING_PRINT("E-mail: mail@log69.com - suggestions & feedbacks are welcome\n");
STRING_PRINT("URL: http://log69.com - the official site\n");
STRING_PRINT("aaphoto (command-line) version - v0.41\n");
STRING_PRINT("aaRGB (color-correction engine) version - v0.64\n");
STRING_PRINT("last update = 20/03/2011\n");
STRING_PRINT("\n");
#ifndef __BMP_ONLY__
STRING_PRINT("The following libraries are used by this program:\n");
#ifdef __OPENMP__
STRING_PRINT("libgomp - OpenMP for parallel programming, http://gcc.gnu.org/onlinedocs/libgomp/\n");
#endif
STRING_PRINT("libjasper - JasPer software, http://www.ece.uvic.ca/~mdadams/jasper/\n");
STRING_PRINT("libjpeg - IJG JPEG software, http://www.ijg.org/\n");
STRING_PRINT("libpng - PNG software, http://www.libpng.org/\n");
STRING_PRINT("libz - Compression library, http://www.zlib.net/\n");
#endif
#ifdef __BMP_ONLY__
#ifdef __OPENMP__
STRING_PRINT("The following libraries are used by this program:\n");
STRING_PRINT("libgomp - OpenMP for parallel programming, http://gcc.gnu.org/onlinedocs/libgomp/\n");
STRING_PRINT("\n");
#endif
STRING_PRINT("BMP_ONLY version. Only BMP format is supported here.\n");
#endif
STRING_PRINT("\n");
}



void PRINT_DESCRIPTION(void){
STRING_PRINT("[DESCRIPTION]\n");
STRING_PRINT("Auto Adjust Photo is a tiny command-line image manipulation tool ");
STRING_PRINT("for automatic color correction of photos. ");
STRING_PRINT("It tries to make the picture look better. ");
STRING_PRINT("The program does this by analyzing the input image and then sets the ");
STRING_PRINT("most optimal contrast, gamma, color balance and saturation for it.\n");
STRING_PRINT("\n");
}



void PRINT_LICENSE(void){
STRING_PRINT("[LICENSE]\n");
STRING_PRINT("This program is free software; you can redistribute it and/or modify ");
STRING_PRINT("it under the terms of the GNU General Public License as published by ");
STRING_PRINT("the Free Software Foundation; either version 3 of the License, or ");
STRING_PRINT("(at your option) any later version.\n");
STRING_PRINT("\n");
STRING_PRINT("This program is distributed in the hope that it will be useful, ");
STRING_PRINT("but WITHOUT ANY WARRANTY; without even the implied warranty of ");
STRING_PRINT("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ");
STRING_PRINT("GNU General Public License for more details.\n");
STRING_PRINT("\n");
STRING_PRINT("You should have received a copy of the GNU General Public License ");
STRING_PRINT("along with this program.  If not, see <http://www.gnu.org/licenses/>.\n");
STRING_PRINT("\n");
}



void PRINT_HELP(void){
STRING_PRINT("[HELP]\n");
STRING_PRINT("USAGE: aaphoto [options] [source files]\n");
STRING_PRINT("\n");
STRING_PRINT("The following image types are supported (thanks to JasPer, JPEG and PNG):\n");
STRING_PRINT("mif, pnm / pgm / ppm, bmp, ras, jp2, jpc, jpg, png\n");
STRING_PRINT("\n");
STRING_PRINT("Quality settings can be applied only to jp2, jpc, jpg formats\n");
STRING_PRINT("\n");
STRING_PRINT("The following options are supported:\n");
STRING_PRINT("    -h   --help          Print this help\n");
STRING_PRINT("    -v   --version       Print version information\n");
STRING_PRINT("    -a   --autoadjust    Auto adjust the colors of the image\n");
STRING_PRINT("    -o   --output        Set output directory\n");
STRING_PRINT("         --overwrite     Overwrite mode, the original source file is replaced\n");
STRING_PRINT("         --jpg           JPEG image output\n");
STRING_PRINT("         --jp2           JPEG 2000 image output\n");
STRING_PRINT("         --png           PNG image output with alpha channel support\n");
STRING_PRINT("         --bmp           BMP image output\n");
STRING_PRINT("    -r   --resize        Resize image taking the longer side in % or pixels\n");
STRING_PRINT("         --rotate90      Rotate image with 90 degrees clockwise\n");
STRING_PRINT("         --rotate180     Rotate image with 180 degrees\n");
STRING_PRINT("         --rotate270     Rotate image with 90 degrees counter-clockwise\n");
STRING_PRINT("         --flipx         Mirror image horizontally\n");
STRING_PRINT("         --flipy         Mirror image vertically\n");
STRING_PRINT("         --noexif        Save image without EXIF info\n");
STRING_PRINT("    -q   --quality       Set image quality from 1 to 100\n");
STRING_PRINT("    -t   --threads       Set number of working threads (default: autodetect)\n");
STRING_PRINT("    -s   --silent        Silent mode, no information printed during operation\n");
STRING_PRINT("         --quiet         ...same as above\n");
STRING_PRINT("    -V   --verbose       Print verbose information about processing\n");
/* STRING_PRINT("    -R   --recursive     Recursive directory search\n"); */
STRING_PRINT("         --test          Print detailed test information into image\n");
STRING_PRINT("\n");
STRING_PRINT("EXAMPLES:\n");
STRING_PRINT("     aaphoto image.jpg\n");
STRING_PRINT("     aaphoto -a -r600 -q85 *.jpg\n");
STRING_PRINT("     aaphoto mydir\n");
STRING_PRINT("     aaphoto -V --resize70% image.png\n");
STRING_PRINT("     aaphoto --quality60 image.jp2\n");
STRING_PRINT("\n");
STRING_PRINT("REMARKS:\n");
STRING_PRINT("- auto adjust parameter is set by default without any other parameters\n");
STRING_PRINT("- _new file name will be generated without --overwrite parameter\n");
STRING_PRINT("- every file is processed on directory input but not recursively\n");
STRING_PRINT("- order of parameters does not matter\n");
STRING_PRINT("- resize value can be set in percentage too\n");
STRING_PRINT("- resize parameter should be less or equal than original\n");
STRING_PRINT("- resize uses the best (and slowest) resampling method\n");
STRING_PRINT("- default jpeg compression quality is 95%\n");
STRING_PRINT("- EXIF information is restored in jpeg images by default\n");
STRING_PRINT("- number of threads is set to the number of online processors by default\n");
STRING_PRINT("\n");
}



int MAIN_RESIZE(void)
{
	/* if there is an Alpha channel beside the RGB too (clrspc_type says) */
	/* then i use the alpha flag to notify the resize function */
	int alpha_flag = 0;
	if ((bitmap_format_png_clrspc_type == 4) || (bitmap_format_png_clrspc_type == 6)){
		alpha_flag = 1; }

	if (opt_resize){

		unsigned long opr = opt_resize;
		/* calculate resize in percentage */
		if (opt_resize_percent){
			if (bitmap_width > bitmap_height) { opr = bitmap_width  * opr / 100; }
			else                              { opr = bitmap_height * opr / 100; }
		}

		/* if the new size is bigger or equal than the original size, then no process is done */
		if (!((opr > bitmap_width) && (opr > bitmap_height))){

			unsigned long new_width, new_height;
			if (bitmap_width > bitmap_height){
				new_width = opr;
				new_height = bitmap_height * opr / bitmap_width;
			}
			else{
				new_width = bitmap_width * opr / bitmap_height;
				new_height = opr;
			}

			if (new_width < 1) new_width = 1;
			if (new_height < 1) new_height = 1;

			if (RESIZE(&bitmap_buffer, &bitmap_width, &bitmap_height, new_width, new_height, alpha_flag)){
				STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; }

		}

/*		else { */
/*			STRING_PRINT("\n"); STRING_PRINTE("error: resize parameter too big, " */
/*				"target size should be less or equal than original\n"); return 1; } */

	}

	/* Rotate and Flip calls */
	if (opt_rotate90) { if (ROTATE90(&bitmap_buffer, &bitmap_width, &bitmap_height, alpha_flag))
			{ STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; } }
	if (opt_rotate180) { if (ROTATE180(&bitmap_buffer, &bitmap_width, &bitmap_height, alpha_flag))
			{ STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; } }
	if (opt_rotate270) { if (ROTATE270(&bitmap_buffer, &bitmap_width, &bitmap_height, alpha_flag))
			{ STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; } }
	if (opt_flipx) { if (FLIPX(&bitmap_buffer, &bitmap_width, &bitmap_height, alpha_flag))
			{ STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; } }
	if (opt_flipy) { if (FLIPY(&bitmap_buffer, &bitmap_width, &bitmap_height, alpha_flag))
			{ STRING_PRINT("\n"); STRING_PRINTE("error: memory allocation error\n"); return 1; } }


	return 0;
}




void MAIN_AARGB(void)
{
/*	AARGB_NORMAL(bitmap_buffer, bitmap_width, bitmap_height); */
	AARGB_MAIN(bitmap_buffer, bitmap_width, bitmap_height, 0, 0, bitmap_width-1, bitmap_height-1, 0, 0, opt_test);
}




int MAIN_RUN(char *file_name)
{
	int result, res;
	char file_name_new [max_char];

	/* store the unix time here in seconds */
	mytime = time(NULL);

	/* print info */ STRING_PRINTV("verbose mode on\n");

	/* print info */ STRING_PRINTV("available processors ");
	/* print info */ STRING_PRINTVD(max_cpus);
	/* print info */ STRING_PRINTV2("\n");
	if (opt_openmp){
		/* print info */ STRING_PRINTV("requested threads ");
		/* print info */ STRING_PRINTVD(max_threads);
		/* print info */ STRING_PRINTV2("\n"); }
	else{
		/* print info */ STRING_PRINTV("no multi processing support\n"); }


	/* print info */ if (!(opt_verbose)) STRING_PRINT(file_name);
	/* print info */ /* STRING_PRINTV(file_name); STRING_PRINTV2("\n"); */
	/* print info */ STRING_PRINTV(file_name); if (opt_verbose) STRING_PRINT("\n");
	/* print info */ STRING_PRINTV("loading\n");

	/* try to load image */
	result = 0;
	result = BITMAP_LOAD(file_name);
	switch (result) {
		case  2:
			STRING_PRINT("\n");
			STRING_PRINTE("error: file ");
			STRING_PRINTE(file_name);
			STRING_PRINTE(" does not exist\n");
			return 1;
			break;
		case  1:
			STRING_PRINT("\n");
			STRING_PRINTE("error: file ");
			STRING_PRINTE(file_name);
			STRING_PRINTE(" cannot be loaded\n");
			return 1;
			break;
	}

	/* determine new file name */
	res = GET_FILE_NAME_NEW(file_name, file_name_new);
	if (res > 0){
		if (res == 1){
			/* error code 1 for existing file name */
			STRING_PRINT("\n");
			STRING_PRINTE("error: file already exists with the new name ");
			STRING_PRINTE(file_name_new); STRING_PRINTE("\n"); return 1; }
		else {
			/* error code 255 for array boundary error */
			STRING_PRINT("\n");
			STRING_PRINTE("error: file name error at file ");
			STRING_PRINTE(file_name);
			STRING_PRINTE("\n");
			return 1; }
	}

	/* print info */ if (!(opt_verbose)) STRING_PRINT(" .");
	if (opt_autoadjust){
		/* print info */ STRING_PRINTV("auto adjusting colors\n");
		MAIN_AARGB();
	}

	/* print info */ if (!(opt_verbose)) STRING_PRINT(" .");
	if ((opt_resize) || (opt_rotate90) || (opt_rotate180) || (opt_rotate270) || (opt_flipx) || (opt_flipy))
		{ if (MAIN_RESIZE()) return 1; }

	/* print info */ if (!(opt_verbose)) STRING_PRINT(" .");
	/* print info */ STRING_PRINTV("saving\n");
	if (BITMAP_SAVE(file_name_new)){ STRING_PRINT("\n");
		STRING_PRINTE("error: file ");
		STRING_PRINTE(file_name);
		STRING_PRINTE(" cannot be saved\n");
		return 1; }

	/* print info */ if (!(opt_verbose)) STRING_PRINT(" done\n");
	/* print info */ STRING_PRINTV("done\n");

	return 0;
}




int MAIN_ARGUMENTS_READ(int argc, char **argv)
{
	int result, opt_wrong, opt_wrong_total, opt_counter;
	int opt_quality_space_flag;
	int opt_threads_space_flag;
	int opt_resize_space_flag;
	int opt_output_space_flag;

	opt_help = 0;
	opt_version = 0;
	opt_autoadjust = 0;
	opt_overwrite = 0;
	opt_noexif = 0;
	opt_jpg = 0;
	opt_jp2 = 0;
	opt_png = 0;
	opt_bmp = 0;
	opt_resize = 0;
	opt_rotate90 = 0;
	opt_rotate180 = 0;
	opt_rotate270 = 0;
	opt_flipx = 0;
	opt_flipy = 0;
	opt_resize_percent = 0;
	opt_output = 0;
	opt_quality = 0;
	opt_threads = 0;
	opt_verbose = 0;
	opt_recursive = 0;
	opt_quiet = 0;
	opt_test = 0;
	opt_doubledash = 0;

	result = 0; /* store the result value of the MAIN_RUN function */

	opt_wrong = 1;
	opt_wrong_total = 0;
	opt_counter = 0;

	opt_quality_space_flag = 0;
	opt_threads_space_flag = 0;
	opt_resize_space_flag = 0;
	opt_output_space_flag = 0;

	if (argc >= 2) {
		char *myarg;

		argc--;
		while (argc--){
			/*printf("%s\n", *argv); */

			opt_wrong = 1;
			myarg = *argv++;
			myarg = *argv;

			if ((opt_doubledash == 0) && ((myarg[0] == '-') || (opt_quality_space_flag) || (opt_threads_space_flag)
				|| (opt_resize_space_flag) || (opt_output_space_flag))) {

				/* check switches */
				if (!STRING_COMPARE(myarg, "--\0")){
					opt_wrong = 0; opt_doubledash = 1; }
				if ((!STRING_COMPARE(myarg, "-h\0")) || (!STRING_COMPARE(myarg, "--help\0"))){
					opt_wrong = 0; opt_help = 1; }
				if ((!STRING_COMPARE(myarg, "-v\0")) || (!STRING_COMPARE(myarg, "--version\0"))){
					opt_wrong = 0; opt_version = 1; }
				if ((!STRING_COMPARE(myarg, "-a\0")) || (!STRING_COMPARE(myarg, "--autoadjust\0"))){
					opt_wrong = 0; opt_autoadjust = 1; }
				if (!STRING_COMPARE(myarg, "--overwrite\0")){
					opt_wrong = 0; opt_overwrite = 1; }
				if (!STRING_COMPARE(myarg, "--noexif\0")){
					opt_wrong = 0; opt_noexif = 1; }
				if (!STRING_COMPARE(myarg, "--jpg\0")){
					opt_wrong = 0; opt_jpg = 1; }
				if (!STRING_COMPARE(myarg, "--jp2\0")){
					opt_wrong = 0; opt_jp2 = 1; }
				if (!STRING_COMPARE(myarg, "--png\0")){
					opt_wrong = 0; opt_png = 1; }
				if (!STRING_COMPARE(myarg, "--bmp\0")){
					opt_wrong = 0; opt_bmp = 1; }
/*				if ((!STRING_COMPARE(myarg, "-R\0")) || (!STRING_COMPARE(myarg, "--recursive\0"))){ */
/*					opt_wrong = 0; opt_recursive = 1; } */
				if ((!STRING_COMPARE(myarg, "-s\0")) || (!STRING_COMPARE(myarg, "--silent\0"))
					|| (!STRING_COMPARE(myarg, "--quiet\0"))){
					opt_wrong = 0; opt_quiet = 1; }
				if ((!STRING_COMPARE(myarg, "-V\0")) || (!STRING_COMPARE(myarg, "--verbose\0"))){
					opt_wrong = 0; opt_verbose = 1; }
				if (!STRING_COMPARE(myarg, "--test\0")){
					opt_wrong = 0; opt_test = 1; opt_autoadjust = 1; }


				/* check --output switch */

				/* check if there was space between the --output switch and its parameter */
				/* and if the parameter needs to be parsed separately */
				if (opt_output_space_flag){
					DIR * dp;
					int char_offs, i;
					char output_ending;

					opt_output_space_flag = 0;

					char_offs = 0;
					i = 0;
					output_ending = '\0';
					while ((myarg[char_offs+i]) && (i < max_char)){
						opt_output_path[i] = myarg[char_offs+i];
						output_ending = opt_output_path[i];
						i++;
					}
					/* check if output path ends with slash character */
					/* if not, then append it */
					if (output_ending != slsh){
						opt_output_path [i] = slsh;
						i++;
					}
					opt_output_path [i] = 0;

					/* check if specified output directory is a dir and if it exists */
					dp = opendir(opt_output_path);
					if (dp == NULL){
						STRING_PRINTE("error: bad directory parameter\n"); return 1; }

					opt_output = 1;
					opt_wrong = 0;
				}

				if ((!STRING_COMPARE_FIX(myarg, "-o", 2)) || (!STRING_COMPARE_FIX(myarg, "--output", 8))){
					int char_offs, i;

					if (!STRING_COMPARE_FIX(myarg, "--output", 8)){
						char_offs = 8;
					}
					else{
						char_offs = 2;
					}

					i = 0;
					/* check if there is space between the switch and its parameter */
					if (myarg[char_offs+i] == 0){
						opt_output_space_flag = 1;
						opt_wrong = 0;
					}
					else{
						DIR * dp;

						char output_ending = '\0';
						while ((myarg[char_offs+i]) && (i < max_char)){
							opt_output_path[i] = myarg[char_offs+i];
							output_ending = opt_output_path[i];
							i++;
						}
						/* check if output path ends with slash character */
						/* if not, then append it */
						if (output_ending != slsh){
							opt_output_path [i] = slsh;
							i++;
						}
						opt_output_path [i] = 0;

						/* check if specified output directory is a dir and if it exists */
						dp = opendir(opt_output_path);
						if (dp == NULL){
							STRING_PRINTE("error: bad directory parameter\n"); return 1; }

						opt_output = 1;
						opt_wrong = 0;
					}

				}

				/* check --quality switch */

				/* check if there was space between the --quality switch and its parameter */
				/* and if the parameter needs to be parsed separately */
				if (opt_quality_space_flag){
					char num [8+1];
					int num_max, num_offs, i;

					opt_quality_space_flag = 0;
					num_max = 8;
					num_offs = 0;
					i = 0;
					while ((myarg[num_offs+i]) && (i < num_max)){
						num[i] = myarg[num_offs+i];
						i++;
					}
					num [i] = 0;
					if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }
					if (STRING_CONVERT_TO_INTEGER(num, &opt_quality)){
						STRING_PRINTE("error: bad parameters\n"); return 1; }

					if ((opt_quality >= 1) && (opt_quality <= 100)){ opt_wrong = 0; }
					else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
				}

				if ((!STRING_COMPARE_FIX(myarg, "-q", 2)) || (!STRING_COMPARE_FIX(myarg, "--quality", 9))){
					char num [8+1];
					int num_max, num_offs, i;

					num_max = 8;
					if (!STRING_COMPARE_FIX(myarg, "--quality", 9)){
						num_offs = 9;
					}
					else{
						num_offs = 2;
					}

					i = 0;
					/* check if there is space between the switch and its parameter */
					if (myarg[num_offs+i] == 0){
						opt_quality_space_flag = 1;
						opt_wrong = 0;
					}
					else{
						while ((myarg[num_offs+i]) && (i < num_max)){
							num[i] = myarg[num_offs+i];
							i++;
						}
						num [i] = 0;
						if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }

						if (STRING_CONVERT_TO_INTEGER(num, &opt_quality)){
							STRING_PRINTE("error: bad parameters\n"); return 1; }
						if ((opt_quality >= 1) && (opt_quality <= 100)){ opt_wrong = 0; }
						else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
					}

				}

				/* check --threads switch */

				/* check if there was space between the --threads switch and its parameter */
				/* and if the parameter needs to be parsed separately */
				if (opt_threads_space_flag){
					char num [8+1];
					int num_max, num_offs, i;

					opt_threads_space_flag = 0;
					num_max = 8;
					num_offs = 0;
					i = 0;
					while ((myarg[num_offs+i]) && (i < num_max)){
						num[i] = myarg[num_offs+i];
						i++;
					}
					num [i] = 0;
					if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }
					if (STRING_CONVERT_TO_INTEGER(num, &opt_threads)){
						STRING_PRINTE("error: bad parameters\n"); return 1; }

					if ((opt_threads >= 1) && (opt_threads <= max_num_threads)){ opt_wrong = 0; }
					else{
						STRING_PRINTE("error: bad parameters, threads value must be between 1 and ");
						STRING_PRINTED(max_num_threads);
						STRING_PRINTE("\n");
						return 1; }

					max_threads = opt_threads;
				}

				if ((!STRING_COMPARE_FIX(myarg, "-t", 2)) || (!STRING_COMPARE_FIX(myarg, "--threads", 9))){
					char num [8+1];
					int num_max, num_offs, i;

					num_max = 8;
					if (!STRING_COMPARE_FIX(myarg, "--threads", 9)){
						num_offs = 9;
					}
					else{
						num_offs = 2;
					}

					i = 0;
					/* check if there is space between the switch and its parameter */
					if (myarg[num_offs+i] == 0){
						opt_threads_space_flag = 1;
						opt_wrong = 0;
					}
					else{
						while ((myarg[num_offs+i]) && (i < num_max)){
							num[i] = myarg[num_offs+i];
							i++;
						}
						num [i] = 0;
						if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }

						if (STRING_CONVERT_TO_INTEGER(num, &opt_threads)){
							STRING_PRINTE("error: bad parameters\n"); return 1; }
						if ((opt_threads >= 1) && (opt_threads <= max_num_threads)){ opt_wrong = 0; }
						else{
							STRING_PRINTE("error: bad parameters, threads value must be between 1 and ");
							STRING_PRINTED(max_num_threads);
							STRING_PRINTE("\n");
							return 1; }

						max_threads = opt_threads;
					}

				}

				/* check --resize switch */

				/* check if there was space between the --resize switch and its parameter */
				/* and if the parameter needs to be parsed separately */
				if (opt_resize_space_flag){
					char num [16+1];
					int num_max, num_offs, i;

					opt_resize_space_flag = 0;

					num_max = 16;
					num_offs = 0;
					i = 0;
					while ((myarg[num_offs+i]) && (myarg[num_offs+i] != '%') && (i < num_max)){
						num[i] = myarg[num_offs+i];
						i++;
					}
					num [i] = 0;
					if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }
					/* is the resize value entered in percentage? */
					if (myarg[num_offs+i] == '%') opt_resize_percent = 1;

					if (STRING_CONVERT_TO_INTEGER(num, &opt_resize)){
						STRING_PRINTE("error: bad parameters\n"); return 1; }

					if (opt_resize_percent){
						if ((opt_resize >= 1) && (opt_resize <= 100)){ opt_wrong = 0; }
						else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
					}
					else{
						if ((opt_resize >= 1) && (opt_resize <= 100000)){ opt_wrong = 0; }
						else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
					}
				}

				if ((!STRING_COMPARE_FIX(myarg, "-r", 2)) || (!STRING_COMPARE_FIX(myarg, "--resize", 8))){
					int i;

					char num [16+1] = "";
					int num_max = 16;
					int num_offs;
					if (!STRING_COMPARE_FIX(myarg, "--resize", 8)){
						num_offs = 8;
					}
					else{
						num_offs = 2;
					}

					i = 0;
					/* check if there is space between the switch and its parameter */

					if (myarg[num_offs+i] == 0){
						opt_resize_space_flag = 1;
						opt_wrong = 0;
					}
					else{
						while ((myarg[num_offs+i]) && (myarg[num_offs+i] != '%') && (i < num_max)){
							num[i] = myarg[num_offs+i];
							i++;
						}
						num [i] = 0;
						if (i >= num_max){ STRING_PRINTE("error: bad parameters\n"); return 1; }
						/* is the resize value entered in percentage? */
						if (myarg[num_offs+i] == '%') opt_resize_percent = 1;

						if (STRING_CONVERT_TO_INTEGER(num, &opt_resize)){
							STRING_PRINTE("error: bad parameters\n"); return 1; }

						if (opt_resize_percent){
							if ((opt_resize >= 1) && (opt_resize <= 100)){ opt_wrong = 0; }
							else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
						}
						else{
							if ((opt_resize >= 1) && (opt_resize <= 100000)){ opt_wrong = 0; }
							else{ STRING_PRINTE("error: bad parameters\n"); return 1; }
						}
					}

				}

				/* check --rotate and --flip switches */
				if (!STRING_COMPARE(myarg, "--rotate90\0")){ opt_wrong = 0; opt_rotate90 = 1; }
				if (!STRING_COMPARE(myarg, "--rotate180\0")){ opt_wrong = 0; opt_rotate180 = 1; }
				if (!STRING_COMPARE(myarg, "--rotate270\0")){ opt_wrong = 0; opt_rotate270 = 1; }
				if (!STRING_COMPARE(myarg, "--flipx\0")){ opt_wrong = 0; opt_flipx = 1; }
				if (!STRING_COMPARE(myarg, "--flipy\0")){ opt_wrong = 0; opt_flipy = 1; }

				/* remember if there were switches but entered wrongly */
				if (opt_wrong){ opt_wrong_total = 1; }
				else{ opt_counter++; }
			}
			else {
				/* expand the file list that contains the files to be processed one by one */
				FILE_LIST_ADD(myarg);
			}

		}

		/* error message if wrong parameters */
		/* error message on missing extra parameter too */
		if ((opt_wrong_total) || (opt_quality_space_flag) || (opt_threads_space_flag) || (opt_resize_space_flag) || (opt_output_space_flag)){
			STRING_PRINTE("error: bad parameters\n"); return 1; }
		else{

			/* check if only 1 type of output format is specified on input */
			int cnt = 0;
			if (opt_jpg) cnt++;
			if (opt_jp2) cnt++;
			if (opt_png) cnt++;
			if (opt_bmp) cnt++;
			if (cnt > 1){
				STRING_PRINTE("error: only one output format allowed\n");
				return 1; }

			/* manage info switches */
			if (opt_version) { PRINT_VERSION(); PRINT_LICENSE(); }
			if (opt_help)    { PRINT_DESCRIPTION(); PRINT_HELP(); }

			/* check if --verbose switch is used alone? */
			cnt = 0;
			cnt += opt_help;
			cnt += opt_version;
			cnt += opt_autoadjust;
			cnt += opt_output;
			cnt += opt_overwrite;
			cnt += opt_jpg;
			cnt += opt_jp2;
			cnt += opt_png;
			cnt += opt_bmp;
			cnt += opt_resize;
			cnt += opt_rotate90;
			cnt += opt_rotate180;
			cnt += opt_rotate270;
			cnt += opt_flipx;
			cnt += opt_flipy;
			cnt += opt_noexif;
			cnt += opt_quality;
			cnt += opt_test;

			if ((opt_verbose) && (cnt <= 0)) {
				STRING_PRINTE("error: verbose option flag needs more switches\n");
				return 1;
			}

			/* AUTOADJUST DEFAULT = 1 IF NO OTHER PARAMETER EXCEPT FILE NAME */
			/* if only file names specified, then --autoadjust parameter is turned on by default */
			if (((opt_counter <= 0) && (file_name_counter > 0)) ||
				((opt_counter == 1) && (opt_doubledash) && (file_name_counter > 0)))
				{ opt_autoadjust = 1; opt_counter++; }

			/* Main process of the switches (load -> process -> save) */
			/* ------------------------------------------------------------------ */
			/* check switches that need file name as input */
			if ((opt_autoadjust) ||
			    (opt_overwrite ) ||
			    (opt_noexif    ) ||
			    (opt_jpg       ) ||
			    (opt_jp2       ) ||
			    (opt_png       ) ||
			    (opt_bmp       ) ||
			    (opt_resize    ) ||
			    (opt_rotate90  ) ||
			    (opt_rotate180 ) ||
			    (opt_rotate270 ) ||
			    (opt_flipx     ) ||
			    (opt_flipy     ) ||
			    (opt_output    ) ||
			    (opt_quality   ) ||
			    (opt_recursive ) ||
			    (opt_test      )){

				/* were there any file names on input? */
				if (file_name_counter > 0){

					/* no error messages during the process with 1 or less file names */
					/*if (file_name_counter > 1) opt_quiet = 1; */

					/* process them one by one when at least 1 or more input files */
					char tt[max_char];
					int c = 0;
					int d;
					int i;
					for (i=0; i<file_name_counter; i++){
						d = 0;
						while (file_name_buffer[c]){
							if (d >= max_char) return 255;
							if (c >= max_file_name_buffer) return 255;
							tt[d] = file_name_buffer[c];
							d++;
							c++;
						}
						tt[d] = 0; d++; c++;
						/* -------------------------------------------------------- */
							if (MAIN_RUN(tt)) { result = 1; }
						/* -------------------------------------------------------- */
					}
				}

				/* error here, because no input file names were specified */
				else {
					STRING_PRINTE("error: file name missing\n");
					return 1;
				}
			}
		}
	}

	/* if no switches or file names or any parameters were specified, */
	/* (the argument counter = 1 which is the command name itself) */
	/* then print something */
	else {
		/*PRINT_VERSION(); */
		/*PRINT_DESCRIPTION(); */
		/*PRINT_LICENSE(); */
		/*PRINT_HELP(); */

		STRING_PRINTE("No parameters specified. Use -h (--help) for more information\n");
	}

	return result;
}




int MAIN_INIT(void)
{
	/* store the unix time here in seconds */
	mytime = time(NULL);

	/* OpenMP */
	/* define the maximum number of threads, this is important to set it in the beginning already */
	/* cause fixed arrays will have to be defined knowing that value */
	max_cpus = 1;
	max_threads = 1;
	max_num_threads = 1024;
	opt_openmp = 0;
	#ifdef __OPENMP__
	/* is there openmp support? */
	opt_openmp = 1;
	/* get the number of available processors to know maximum number of threads */
	max_cpus = omp_get_num_procs();
	max_threads = max_cpus;
	/* set it to minimum 1 */
	if (max_threads < 1){ max_threads = 1; }
	#endif

	/* number of files */
	file_counter = 0;
	/* number of files stored in file_name_buffer */
	file_name_counter = 0;
	/* pointer pointing to the end of the file names read from the files that are stored in file_name_buffer */
	/* it is the offset pointer pointing to the next file name */
	file_name_buffer_pointer = 0;
	/* allocate memory */
	file_name_buffer = malloc(max_file_name_buffer * sizeof (*file_name_buffer));
	if (file_name_buffer == 0){
		STRING_PRINTE("error: memory allocation failure\n"); return 1; }

#ifndef __BMP_ONLY__
	bitmap_format_png_interlace_type = PNG_INTERLACE_NONE;
	bitmap_format_png_compression_type = PNG_COMPRESSION_TYPE_DEFAULT;
	bitmap_format_png_filter_type = PNG_FILTER_TYPE_DEFAULT;
#endif

	/* set default image resolution */
	/* 2835 dots per meter is equal to 72 dots per inch */
	xdpi = 2835;
	ydpi = 2835;
	udpi = 1;

	return 0;
}



int main(int argc, char **argv)
{
	int result;

	if (MAIN_INIT()) return 1;

	result = 0;
	if (MAIN_ARGUMENTS_READ(argc, argv)){ result = 1; }

	free(file_name_buffer);

	return result;
}
