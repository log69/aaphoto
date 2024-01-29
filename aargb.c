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





/* --------------------------------------------------- */
/* ----------- Auto Adjust RGB ----------------------- */
/* ----------- András Horváth (C) 2006-2016 ---------- */
/* ----------- Hungary, http://log69.com ------------- */
/* --------------------------------------------------- */

/*

aaRGB Changelog:
------------------
2013/09/30 - aaRGB v0.65 - fix some compile time warnings and uninitalized variables
2011/01/26 - aaRGB v0.64 - the contrast seemed strong with the default initial value of the former algorithm,
                           so now constant optimized to the new contrast algorithm introduced in 0.63
                         - add more detailed explanations to the contrast algorithm (in hungarian)
2010/12/18 - aaRGB v0.63 - include OpenMP defs in aargb.c too for standalone usage
                         - fix some warning messages during build
                         - make embedded test info image look more readable by changing the colors (--test switch)
                         - change saturation algorithm from linear to exponential
                         - improve contrast algorithm to include a self balance mechanism to avoid overexposure
                           on images with large blank areas that have relatively small details
                         - improve color balance algorithm to make it a bit more aggressive by raising the value
                           of the color difference factor of the white and black point to the 3rd power
                         - speed up process by skipping saturation part if no change is needed
                         - some memory allocation check
                         - some changes in documentation
2010/09/14 - aaRGB v0.62 - bugfix: an ugly misconception in my paralleled code caused weird behavior
                           when using more threads
                         - rewrite the code to suffice the ISO C90 ANSI standard C form (GCC -pedantic option)
2010/05/02 - aaRGB v0.61 - add OpenMP support for multi processing, all computing cycles paralleled
                         - solve warning issues with some uninitialized variables
                         - some code cleanup
2009/10/18 - aaRGB v0.60 - remove gamma handling of the lighter colors from the two-pole gamma computing
                           by setting the gamma_interval_high from 0.9 to 1, it proved to be inefficient
2009/04/05 - aaRGB v0.59 - some more code cleanup
2009/02/22 - aaRGB v0.58 - code cleanup
2007/08/11 - aaRGB v0.57 - improve black and white point analyzing
                           from now they are not scaled to perfect black and white, but to their darkest
                           and brightest color that have maximum saturation to fix overexposure problem
                         - improve saturation algorithm with full floating point computing and HSL conversion
                           to fix over saturated colors
                         - expand image information display with color balance circle for testing (--test switch)
                         - remove text information from test display (--test switch)
2007/06/21 - aaRGB v0.56 - expand functionality with "apply only on selection" to process only the selected area
2007/04/03 - aaRGB v0.55 - maximize saturation limit with a predefined constant to avoid overexposure in saturation
                           when reconverting the same image
2007/04/01 - aaRGB v0.54 - new two-pole gamma computing
                         - new saturation compensation
2007/03/29 - aaRGB v0.53 - improve contrast computing to avoid underexposure
2007/02/25 - aaRGB v0.52 - improve image information display for testing (--test switch)
2007/02/16 - aaRGB v0.51 - improve average RGB color computing for more accurate color balance calibration
2007/01/04 - aaRGB v0.49 - stable working version with gamma handling and more clever image analyzing
2006/08/29 - aaRGB project begun...

aaRGB end of Changelog.

*/


/* ---------------------------------------------------------------- */
/* -- Automatically adjusts contrast, color balance, gamma level -- */
/* -- and saturation of an RGB image ------------------------------ */
/* ---------------------------------------------------------------- */

/* 
   pictures might need some contrast adjustment when their colors don't spread
   the entire spectrum, usually when their darkest colors aren't dark enough and 
   the lightest colors not bright enough, this causing a veiled effect on them.
   
   sometime color balance adjustment is also necessary, when the brightest and darkest
   colors of the image are not black and white enough, but they are shifted towards
   some other color which may be the effect of the technology used,
   or simply the environment changes the colors unrealistically.
   
   the procedure assumes that there must be some detail on the picture that's
   supposed to be black (dark enough) and also white (bright enough) in real.
   the automatic contrast and color balance process is built upon this.
   
   the overall (or average) brilliance of the image can be set by the gamma level.
   this procedure makes sure that those pictures that seem too dark will be
   corrected by this factor by changing the middle dark colors while not changing
   the darkest and brightest colors.
   
   as the last step, a saturation adjustment takes place to intensify the colors
   themselves if needed.
*/



#include <math.h>
#include <stdlib.h>



/* ------------------------- */
/* ----- SIGN FUNCTION ----- */
/* ------------------------- */

int sgn(double x)
{
    if (x <  0) { return -1; }
    if (x >  0) { return  1; }
/*  if x is 0 then return 0 */
    return 0;
}



/* -------------------------------- */
/* ----- RGB - HSL CONVERSION ----- */
/* -------------------------------- */

void RGB_TO_HSL(double R, double G, double B, double *H1, double *S1, double *L1)
{
    double H, S, L;
    double LN, LK, LNX, LKX;

    if (R < 0) { R = 0; } if (R > 1) { R = 1; }
    if (G < 0) { G = 0; } if (G > 1) { G = 1; }
    if (B < 0) { B = 0; } if (B > 1) { B = 1; }

    H = 0;
    S = 0;
    L = 0;

    L = (R + G + B) / 3;

    /* To calculate the 'S' saturation value, first i determine the
       maximum amount of how much i can stretch the RGB color elements
       to their limits while keeping the same 'L' lightness and 'H' hue values */

    LN = R;
    if (LN < G) { LN = G; }
    if (LN < B) { LN = B; }
    LK = R;
    if (LK > G) { LK = G; }
    if (LK > B) { LK = B; }
    if (LN == LK) { S = 0; H = 0; }
    else {
        double R2, G2, B2;

        if ((LN < 1) && (LK > 0)) {
            /* L cannot be either 0 or 1 here */
            LKX = (L - LK) / L;
            LNX = (LN - L) / (1 - L);
            S = LNX;
            if (LKX >= LNX) { S = LKX; }
        }
        else { S = 1; }

        /* To get the 'H' color value, i stretch the RGB elements to their maximum,
           so it'll be a 6 case outcome. */

        LN = LN - LK;
        R2 = (R - LK) / LN;
        G2 = (G - LK) / LN;
        B2 = (B - LK) / LN;

        if ((R2 == 1) && (G2 <  1) && (G2 >= 0) && (B2 == 0))  { H = (  0 +    G2  * 60) / 360; }
        if ((G2 == 1) && (R2 <= 1) && (R2 >  0) && (B2 == 0))  { H = ( 60 + (1-R2) * 60) / 360; }
        if ((G2 == 1) && (B2 <  1) && (B2 >= 0) && (R2 == 0))  { H = (120 +    B2  * 60) / 360; }
        if ((B2 == 1) && (G2 <= 1) && (G2 >  0) && (R2 == 0))  { H = (180 + (1-G2) * 60) / 360; }
        if ((B2 == 1) && (R2 <  1) && (R2 >= 0) && (G2 == 0))  { H = (240 +    R2  * 60) / 360; }
        if ((R2 == 1) && (B2 <= 1) && (B2 >  0) && (G2 == 0))  { H = (300 + (1-B2) * 60) / 360; }

        if (H == 1) { H = 0; }
    }

    *H1 = H; *S1 = S; *L1 = L;

}



void HSL_TO_RGB(double H, double S, double L, double *R1, double *G1, double *B1)
{
    double R = 0;
    double G = 0;
    double B = 0;

    if (H < 0) { H = 0; } if (H > 1) { H = 1; }
    if (S < 0) { S = 0; } if (S > 1) { S = 1; }
    if (L < 0) { L = 0; } if (L > 1) { L = 1; }

    if (L == 0) { R = 0; G = 0; B = 0; }
    else {
        if (L == 1) { R = 1; G = 1; B = 1; }
        else {
            double L2, templ;

            /* set the things here based on 'H' value */
            int deg;
            double mul;
            if (H == 1) { H = 0; }
            deg = (int)(H * 6);
            mul = H * 6 - deg;

            switch (deg) {
                case 0:
                    R = 1; G = mul;   B = 0; break;
                case 1:
                    G = 1; R = 1-mul; B = 0; break;
                case 2:
                    G = 1; B = mul;   R = 0; break;
                case 3:
                    B = 1; G = 1-mul; R = 0; break;
                case 4:
                    B = 1; R = mul;   G = 0; break;
                case 5:
                    R = 1; B = 1-mul; G = 0; break;
            }

            /* scale the things here in the RGB field by the value of 'L' */
            L2 = (R + G + B) / 3;
            if (L > L2) {
                templ = (1-L) / (1-L2);
                R = 1 - (1-R) * templ;
                G = 1 - (1-G) * templ;
                B = 1 - (1-B) * templ;
            }
            else {
                templ = L / L2;
                R = R * templ;
                G = G * templ;
                B = B * templ;
            }

            /* scale the things here in the RGB field by the value of 'S' */
            if (R > L) { R = L + (R-L) * S; } else { R = L - (L-R) * S; }
            if (G > L) { G = L + (G-L) * S; } else { G = L - (L-G) * S; }
            if (B > L) { B = L + (B-L) * S; } else { B = L - (L-B) * S; }
        }
    }

    *R1 = R; *G1 = G; *B1 = B;

}



/* ------------------------------------------------------------------------------ */

void AARGB_MAIN(
    unsigned char *image_buffer,
    int image_width,
    int image_height,
    int x1,
    int y1,
    int x2,
    int y2,
    int format_flag,
    int apply_on_selection,
    int test_flag)
{


/* ------------------------------------------------------------------------------ */
/* ----------  Global variables for the procedure  ------------------------------ */
/* ------------------------------------------------------------------------------ */

      int max_threads2;

      double cont_max;
      double gamma_max;
      double gamma_interval_low;
      double gamma_interval_high;
      double satur_max;

      double gamma_weight_low_all;
      double gamma_weight_high_all;
      double gamma_weight_low;
      double gamma_weight_high;
      double gamma_low;
      double gamma_high;

      long hist1[256];
      long hist2[256];
      long hist3[256];
      long *hist1n;
      long *hist2n;
      long *hist3n;
      long hist_min;
      long hist_max;
      long hist_cut_columns;
      long hist_cut_weight;
      double hist_cut_limit;
      double hist_avg;

/*
      long hist_sum;
      long hist_min_test;
      long hist_max_test;
      double hist_avg_test;
*/

      long hist_satur[256];
      long  *hist_saturn;
      double hist_satur_avg;
      double hist_satur_low;
      double hist_satur_ok;

      double temp1, temp2, temp3;
      long flag1;
      long bw, bh;
      long x, y;
      long i1, i2, i3;

      long col_r, col_g, col_b;
      double col_r2, col_g2, col_b2;
      unsigned long col_r3[256];
      unsigned long col_g3[256];
      unsigned long col_b3[256];
      unsigned long *col_r3n;
      unsigned long *col_g3n;
      unsigned long *col_b3n;

      double H, S, L;

      double wp, bp;
      double wp_end, bp_end;
      double wp_r, wp_g, wp_b;
      double bp_r, bp_g, bp_b;
      double wp_r_end = 0, wp_g_end = 0, wp_b_end = 0;
      double bp_r_end = 0, bp_g_end = 0, bp_b_end = 0;

      long cc;
      long addr, addr2;
      long addr_offset;
      int N;


      long col;
      long xm, ym;
      long xma1, xma2;
      long xmb1, xmb2;
      long color_black =       0x00000000;
      long color_green =       0x0000ff00;
      long color_brown =       0x00707000;
      long color_blue =        0x000080ff;
      long color_white =       0x00ffffff;
      long color_gray =        0x00606060;
      long color_red =         0x00a00000;
      long color_yellow =      0x00ffff00;
/*    long color_green_dark =  0x00008000; */


	/* is there openmp support? */
	max_threads2 = 1;
	#ifdef __OPENMP__
        /* get the number of available processors to know maximum number of threads */
        max_threads2 = omp_get_num_procs();
        /* set it to minimum 1 */
        if (max_threads2 < 1){ max_threads2 = 1; }
	#endif


	/* allocate memory for the arrays that help to parallel cycles */
	/* the allocated memory here is way too small, so i don't print any error message on failure
	   cause i don't want to have stdio.h as a dependency for this procedure,
	   just simply exit */

	hist_saturn = 0;

	hist1n = 0;
	hist2n = 0;
	hist3n = 0;

	col_r3n = 0;
	col_g3n = 0;
	col_b3n = 0;

	if ((hist_saturn = calloc (256 * max_threads2, sizeof (*hist_saturn))) == 0){ goto exit; }

	if ((hist1n = calloc (256 * max_threads2, sizeof (*hist1n))) == 0){ goto exit; }
	if ((hist2n = calloc (256 * max_threads2, sizeof (*hist2n))) == 0){ goto exit; }
	if ((hist3n = calloc (256 * max_threads2, sizeof (*hist3n))) == 0){ goto exit; }

	if ((col_r3n = calloc (256 * max_threads2, sizeof (*col_r3n))) == 0){ goto exit; }
	if ((col_g3n = calloc (256 * max_threads2, sizeof (*col_g3n))) == 0){ goto exit; }
	if ((col_b3n = calloc (256 * max_threads2, sizeof (*col_b3n))) == 0){ goto exit; }


/* ------------------------------------------------------------------------------ */
/* Initialization and constants for the contrast and gamma (0...1) */
/* ------------------------------------------------------------------------------ */
/* Initial values and constants definition for contrast and gamma operations */

/* Contrast constant definition: specifies the maximum value for automatic */
/* contrast adjustment (ranging from 0 to 1), */
/* default = 0.1 */
/* Gamma constant definition: specifies the maximum value for automatic */
/* gamma adjustment (recommended range 1 to 10), */
/* default = 1.5 */

/* For gamma adjustment, specifying the maximum occurrence value (occur_max), */
/* which indicates that in calculating gamma value, if the occurrence of */
/* predominant colors is less than this, they are omitted from */
/* the calculation, meaning that predominant colors are not considered */
/* because the focus is on making details visible, */
/* i.e., the brightness is calculated preferably for details */

/* ------------------------------ */
/* Setting constant values */
/* ------------------------------ */
    /* maximum contrast level (0..1) */
    cont_max = 0.066666;

    /* maximum gamma correction level (1..10) */
    gamma_max = 1.5;
    gamma_interval_low = 0.333;
    gamma_interval_high = 1;

    /* maximum saturation limit (0..1) */
    satur_max = 0.333;

    bw = image_width;
    bh = image_height;

    /* Checking and setting boundary values of coordinates for the selected area */
    /* The purpose of the selection is to make a specific part of the image */
    /* visible (not to proportionally adjust the entire image) */

    /* converting width to absolute coordinates, currently inactive */
    /* x2 = x1 + x2 - 1; */
    /* y2 = y1 + y2 - 1; */

    if (x1 < 0) { x1 = 0; }
    if (x2 < 0) { x2 = 0; }
    if (y1 < 0) { y1 = 0; }
    if (y2 < 0) { y2 = 0; }
    if (x1 > bw-1) { x1 = bw-1; }
    if (x2 > bw-1) { x2 = bw-1; }
    if (y1 > bh-1) { y1 = bh-1; }
    if (y2 > bh-1) { y2 = bh-1; }
    /* Calculating the offset value for the 4-byte alignment in DIB format */
    /* This is not necessary for a normal array, in which case the format_flag value = 0 */
    /* otherwise, in the DIB format, the rows are vertically inverted, */
    /* and the row ends are aligned with 4 bytes */
    addr_offset = 0;
    /* Format = 0 --> NORMAL 3 byte RGB data in array */
    /* Format = 1 --> DIB data format in array */
    /* Format = 2 --> BMP data format in array */
    if ((format_flag == 1) || (format_flag == 2)){
        y1 = bh - 1 - y1;
        y2 = bh - 1 - y2;
        addr_offset = bw * 3 - 4 * (bw * 3 / 4);
        if (addr_offset) addr_offset = 4 - addr_offset;
    }
    /* Swapping selection coordinates if they are not given in correct order, */
    /* i.e., the top left corner is x1 and y1, the bottom right is x2 and y2 */
    if (x1 > x2) { i1 = x1; x1 = x2; x2 = i1; }
    if (y1 > y2) { i1 = y1; y1 = y2; y2 = i1; }



/* ------------------------------------------------------------------------------ */
/* Create Histogram and average RGB colors for the Image */
/* ------------------------------------------------------------------------------ */
/* Generating a histogram from the colors of the image, which represents the distribution of brightness in the image */
/* plus storing the average RGB values of colors with the same brightness, in anticipation of the work of a future routine */

/* The average brightness (grayness) of each point in the image is included in a */
/* 256-element array, where the brightness value increases the value of the corresponding index in the array by 1 */

/* Thus, this histogram array gives an accurate description of the color scale of the image from black to white */

    /* array initialization with zeros */
    #ifdef __OPENMP__
    #pragma omp parallel for num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){
        hist1[i1] = 0;
        hist2[i1] = 0;
        hist3[i1] = 0;
        hist_satur[i1] = 0;
        col_r3[i1] = 0;
        col_g3[i1] = 0;
        col_b3[i1] = 0;
    }

    #ifdef __OPENMP__
    #pragma omp parallel for private(x, y, addr, addr2, col_r, col_g, col_b, cc, N) num_threads(max_threads2)
    #endif
    for (y=0; y<=bh-1; y++){
        addr2 = y * bw * 3 + y * addr_offset;
        for (x=0; x<=bw-1; x++){
            if (x >= x1 &&
                x <= x2 &&
                y >= y1 &&
                y <= y2) {

                addr = addr2 + x * 3;
                col_r = image_buffer[addr + 0];
                col_g = image_buffer[addr + 1];
                col_b = image_buffer[addr + 2];
                cc = (col_r + col_g + col_b) / 3;
		/* multi processing, thread id */
#ifdef __OPENMP__
		N = omp_get_thread_num();
#else
		N = 0;
#endif
                /* gray histogram */
                hist1n[cc + N*256]++;
                /* Storing average RGB values for determining the average RGB of black and white points later */
                col_r3n[cc + N*256] += col_r;
                col_g3n[cc + N*256] += col_g;
                col_b3n[cc + N*256] += col_b;
            }
        }
    }

    /* this is a replacement code for arrays for multi processing */
    #ifdef __OPENMP__
    #pragma omp parallel for private(i1, i2) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){
	for (i2=0; i2<max_threads2; i2++){
		hist1[i1]  += hist1n[i1 + i2*256];
		col_r3[i1] += col_r3n[i1 + i2*256];
		col_g3[i1] += col_g3n[i1 + i2*256];
		col_b3[i1] += col_b3n[i1 + i2*256];
	}
    }


/* ------------------------------------------------------------------------------ */
/* Start analyzing to find the White and Black points */
/* ------------------------------------------------------------------------------ */
/* Determining the Black and White points for automatic contrast adjustment */

/* I start reading the magnitude of gray values from the left and right sides of the histogram */
/* and keep reading until it exceeds a predetermined value, at which point I obtain */
/* the positions of the black and white points */

/* This predetermined value is the product of the contrast constant and the average */
/* of all elements in the histogram */
/* The average product creates a balance at the limit value, because otherwise */
/* if this value were, say, the maximum value of the histogram, then */
/* it would characterize the function with drastic contrast over-adjustment */

/* In other words, this value indicates what percentage of the average brightness of the image */
/* is the threshold beyond which blacks and whites found below this percentage will be clipped to the limit (cut off) */

    hist_min = bw * bh;
    hist_max = 0;
    hist_avg = 0;

    /* Calculating average brightness: I get this by summing up all columns of the histogram 
    and dividing by 256 (number of columns), that is the mathematical average.
    Then I further divide it by a constant (to 10% by default) by multiplying with a number less than zero,
    the resulting value here is called the limit.

    (calculation of maximum and minimum values is only for testing purposes)
    */
    for (i1=0; i1<256; i1++){
        temp1 = hist1[i1];
        if (hist_min > (long)temp1){ hist_min = (long)temp1; }
        if (hist_max < (long)temp1){ hist_max = (long)temp1; }
        hist_avg = hist_avg + temp1;
    }
    /* total sum of the histogram */
/*    hist_sum = hist_avg; */

    /* mathematical average value of the histogram */
    hist_avg = hist_avg / 256;

/*
    hist_min_test = hist_min;
    hist_max_test = hist_max;
    hist_avg_test = hist_avg;
*/

    /* this is the limit */
    temp1 = hist_avg * cont_max;
    hist_cut_limit = temp1;



    /* searching for white and black points up to 0 (zero) maximum */

    bp = 255;
    flag1 = 0;

    for (i1=0; i1<256; i1++){
        if (flag1 == 0){
            if (hist1[i1] > 0){
                flag1 = 1;
                bp = i1;
            }
        }
    }

    wp = 0;
    flag1 = 0;
    hist_cut_columns = 0;

    for (i1=255; i1>=0; i1--){
        if (flag1 == 0){
            if (hist1[i1] > 0){
                flag1 = 1;
                wp = i1;
            }
        }
    }

    if (bp > wp){
        i1 = (long)(wp);
        wp = bp;
        bp = i1;
    }
    if (bp == wp){
       bp = bp - 1;
       wp = wp + 1;
    }
    if (bp < 0){ bp = 0; }
    if (bp >= 255){ bp = 254; }
    if (wp > 255){ wp = 255; }
    if (wp <= 0){ wp = 1; }


    /* -----------------------------------------------------------------------------------
    The automatic contrast adjustment algorithm sometimes receives inputs,
    where the details are too sparse in the image, and large areas of a single color dominate - in these cases, the phenomenon was
    that the algorithm cuts too much from the edges of the histogram, and the important small details
    on these images are lost.
    
    Therefore, I have incorporated a self-regulating mechanism into the algorithm, which, in the event of the above situation,
    increases the likelihood that the details will be preserved after contrast adjustment.
    
    This is done by checking how many columns fall below the limit value, and depending on this,
    I further reduce this limit, because the more portions that fall below the average (number of columns in the histogram),
    it means that there are more interesting details on the image that would be cut off,
    hence the many and significant large empty areas should matter less.

    The limit value is further reduced as follows: I recalculate the average
    and this limit value in such a way that during the calculation, the column values above this
    count less if there are more columns below the previous limit value,
    and if the sum of the heights of these columns is also smaller,
    because this tells us that the valuable detail is located on an increasingly smaller area.
    
    This results in a self-braking process for the too high limit value caused by the average of large empty areas,
    and thereby to the too drastic contrast cut, where the detail is lost.

    After this, I recalculate this limit value in the previous manner,
    and this will be used for further calculation of the histogram edge cuts.
    
    To what extent the small image detail should count, I want to decide in an analogous manner,
    so that the change in the original limit value does not happen in a stepwise fashion,
    but since I believe I need a curve that rises very little until the midpoint of the 0-1 interval, and then more drastically,
    I found multiplying by the 5th power times 3 to be appropriate.
    
    f(x) = x^5*3
    Wolfram Alpha link for illustration:
    http://www.wolframalpha.com/input/?i=x^5+*+3+from+0+to+1
    
    According to this, if the balance between small and large details is upset, then
    from 0.5 upwards, the contrast calculation becomes progressively less drastic,
    and this is how I achieve that the smaller values influence as little as possible,
    while the larger values influence more and more.
    -------------------------------------------------------------------------------------- */


    /* I only examine the section between the black and white points
       (this is the remaining width of the histogram after the cuts from the left and right sides)
       because I want to examine the proportion of detail to the total extent
       of the intended final result (how much would remain of it).
    */

    hist_cut_columns = 0;
    hist_cut_weight = 0;
    for (i1=bp; i1<wp; i1++){
        if (hist1[i1] < temp1){
            hist_cut_columns++;
            hist_cut_weight += hist1[i1];
        }
    }


    /* temp1 shows the original limit value to be cut,
       hist_cut_columns shows the number of columns below the limit (I call this detail below the limit),
       hist_cut_weight shows the sum of the values of these columns (their weight),

       temp2 is set in such a way that the greater the weight spread over fewer column numbers
       (i.e., the fewer columns are below the limit and their weight is greater),
       the greater its value will be - so if temp2 has a greater value,
       it means that there are more valuable details below the limit.
       
       I call the columns below the limit more valuable details, because the average of the entire histogram
       is skewed upwards by the large single-color areas, which are few columns with high weight,
       i.e., these are obviously large extent empty areas - hence they are the "non" details,
       while I consider those below the average to be the detail.
       
       Since we are looking at the columns below the limit, dividing by the limit value itself gives a
       ratio value in the 0..1 interval.
    */
    if ((hist_cut_columns == 0) || (temp1 == 0)){ temp2 = 0; }
    else { temp2 = (double)(hist_cut_weight) / hist_cut_columns / temp1; }


    /* here I determine the value of temp3 by dividing the number of columns below the original limit
   by the remaining (post-cut) middle width of the histogram (wp-bp)
   (but only within the remaining section, therefore the result will be between 0..1),

   meaning the more parts that fall off on the left and right sides, the smaller the value we divide by,
   and if the number of columns below the limit increases, dividing by this smaller value
   will yield a larger number, hence we can consider the originally applied contrast to be more drastic,
   and thus raising this value to the 5th power and multiplying by 3 results in a value
   that gives back an increasingly larger value above 50%, and the larger this value,
   the more I reduce the original planned contrast (cut) measure.
    */
    temp3 = (double)(hist_cut_columns) / (wp-bp);

    /* here I curve the linear value by raising it to a power, so that smaller values have less,
    while larger values increasingly influence the result */
    temp3 = temp3 * temp3 * temp3 * temp3 * temp3 * 3;
    if (temp3 > 1){ temp3 = 1; }

    /* here I pull back the limit line created by temp2 (new limit value) under temp1 (original limit value)
    upwards towards the original temp1 based on the temp3 curve.
    
    meaning I reduced the drastic contrast, and then allow it back based on the curve
    (which curves drastically after about 50%).
    */
    temp3 = temp1 - ((temp1 - temp2) * temp3);
    /* here I prevent the new reduced limit value from going below 10% of the original,
    this is just a lower limit for the extent of contrast reduction */
    if (temp3 < temp1 * 0.1){ temp1 = temp1 * 0.1; }
    else{ temp1 = temp3; }

    hist_cut_limit = temp1;

/* ------------------------------------------------------------ */


    bp = 255;
    flag1 = 0;

    /* I step from the right side of the histogram and cut off up to the value,
    which is still smaller than the average * cont_max (10% of the average) */
    for (i1=0; i1<256; i1++){
        if (flag1 == 0){
            if (hist1[i1] >= temp1){
                flag1 = 1;
                bp = i1;
            }
        }
    }

    wp = 0;
    flag1 = 0;
    hist_cut_columns = 0;

    for (i1=255; i1>=0; i1--){
        if (flag1 == 0){
            if (hist1[i1] >= temp1){
                flag1 = 1;
                wp = i1;
            }
        }
    }

    /* Setting and correcting boundary values */
    if (bp > wp){
        i1 = (long)(wp);
        wp = bp;
        bp = i1;
    }
    if (bp == wp){
       bp = bp - 1;
       wp = wp + 1;
    }
    if (bp < 0){ bp = 0; }
    if (bp >= 255){ bp = 254; }
    if (wp > 255){ wp = 255; }
    if (wp <= 0){ wp = 1; }

    bp = bp / 255;
    wp = wp / 255;



/* ------------------------------------------------------------------------------ */
/* Get the average RGB values for the White and Black points */
/* ------------------------------------------------------------------------------ */
/* Calculating the average RGB values for color balance adjustment */
/* based on the black and white point values */

/* Here, an average RGB value is generated for both the black and white points */
/* This value indicates how each color skews downward (average RGB of the black point) */
/* and upward (average RGB of the white point) during automatic contrast adjustment */

/* The average RGB of all colors above the white point will be the reference point */
/* for distortion towards white, */
/* i.e., this will be extended to perfect white */

/* In reality, this sets the color balance of the image appropriately so that the */
/* color average of the whites to be cut off becomes the perfect white, */
/* therefore, if their average is not perfect white, then the deviation from it */
/* should be proportionally shifted towards perfect white for all colors, */
/* and the same applies to the black point */

    bp_r = 0;
    bp_g = 0;
    bp_b = 0;
    wp_r = 0;
    wp_g = 0;
    wp_b = 0;

    i3 = 0;
    /* Calculating the average RGB of all colors below the black point */
    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:bp_r, bp_g, bp_b, i3) num_threads(max_threads2)
    #endif
    for (i1=(long)(bp * 255); i1>=0; i1--){
        bp_r += col_r3[i1];
        bp_g += col_g3[i1];
        bp_b += col_b3[i1];
        i3 += hist1[i1];
    }

    if (i3 > 0){
        bp_r = bp_r / i3;
        bp_g = bp_g / i3;
        bp_b = bp_b / i3;
    }

    i3 = 0;
    /* Calculating the average RGB of all colors above the white point */
    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:wp_r, wp_g, wp_b, i3) num_threads(max_threads2)
    #endif
    for (i1=(long)(wp * 255); i1<256; i1++){
        wp_r += col_r3[i1];
        wp_g += col_g3[i1];
        wp_b += col_b3[i1];
        i3 += hist1[i1];
    }
    if (i3 > 0){
        wp_r = wp_r / i3;
        wp_g = wp_g / i3;
        wp_b = wp_b / i3;
    }

    /* Scaling from 255 to the [0..1] interval */
    bp_r = bp_r / 255;
    bp_g = bp_g / 255;
    bp_b = bp_b / 255;
    wp_r = wp_r / 255;
    wp_g = wp_g / 255;
    wp_b = wp_b / 255;

    /* Restoring the brightness of the obtained average RGB value to the white point level. */
    /* Since we did not only calculate the average color for colors with the same brightness as the white point */
    /* but for all colors brighter than that, the brightness of the obtained */
    /* average color will be greater or equal to the starting white point */
    /* therefore, for correction, we restore the brightness of the RGB value */
    /* while maintaining the ratio of the R, G, and B components */
    /* and the same applies to the black point */
    RGB_TO_HSL (bp_r, bp_g, bp_b, &H, &S, &L);
    L = bp;
    HSL_TO_RGB (H, S, L, &bp_r, &bp_g, &bp_b);
    RGB_TO_HSL (wp_r, wp_g, wp_b, &H, &S, &L);
    L = wp;
    HSL_TO_RGB (H, S, L, &wp_r, &wp_g, &wp_b);

    /* Calculation of the target points for the black and white points. */
    /* This shows where to move the average RGB of the black and white points */
    /* so that it is parallel to the 'gray' line connecting the two vertices in the RGB cube, */
    /* maintaining the distance from the 'gray' line and the direction of the color, */
    /* but as dark or as bright as possible */
    /* */
    /* In other words, we shift along the gray line until we hit the wall of the RGB cube (in both directions) */
    /* */
    /* This is a change from the previous methods in that the black point is */
    /* now not just pulled to perfect black, but to the corresponding */
    /* darkest point where the color saturation is maximal */
    /* and whose color matches that of the black point, */
    /* thereby eliminating poor color balance and inappropriate contrast. */

    temp3 = bp_r;
    if (temp3 > bp_g) { temp3 = bp_g; }
    if (temp3 > bp_b) { temp3 = bp_b; }
    /*bp_r_end = 1 - (1 - bp_r) / (1 - temp3); */
    /*bp_g_end = 1 - (1 - bp_g) / (1 - temp3); */
    /*bp_b_end = 1 - (1 - bp_b) / (1 - temp3); */
    bp_r_end = bp_r - temp3;
    bp_g_end = bp_g - temp3;
    bp_b_end = bp_b - temp3;

    temp3 = wp_r;
    if (temp3 < wp_g) { temp3 = wp_g; }
    if (temp3 < wp_b) { temp3 = wp_b; }
    if (temp3 > 0){
        wp_r_end = wp_r / temp3;
        wp_g_end = wp_g / temp3;
        wp_b_end = wp_b / temp3;
    }



/* ----------------- */
/* ---- RGB SPACE ---- */
/* ----------------- */
/* The entire RGB space is encompassed by a regular 3D cube, */
/* where each edge represents the R, G, and B axes */
/* and in one corner lies the perfect white color, */
/* while the opposite corner holds the perfect black */
/* and the line connecting these two vertices contains the */
/* complete gray scale from black to white. */
/* The measure of color shift is nothing but the distance of the given color */
/* perpendicular from the gray scale line */
/* (which can be characterized by a value between 0 and 1, where 0 */
/* means the color is gray, i.e., it is located on the line) */
/* */
/* ------------------------- */
/* ---- DIRECTION OF RGB COLORS ---- */
/* ------------------------- */
/* I compare the average RGB shift values of the black and white points, */
/* to determine if they are shifted in the same direction, */
/* because if the color balance of the picture is off, I assume that */
/* all colors of the picture are shifted in the same direction due to the factor */
/* causing the imbalance. If their shifts do not point in the same direction, */
/* I assume that it is not due to a disturbed color balance. In this case, during the contrast operation, I do not */
/* apply color balance compensation (i.e., correction of the picture's colors */
/* towards the direction of perfect white and perfect black). */
/* */
/* The direction of a given color in the RGB cube refers to the angle of rotation */
/* of the perpendicular segment starting from its point in the cube towards the gray line. */
/* This value should fall between -180 and 180 degrees. */
/* Thus, the direction of the average RGB colors of the black and white points gives two angles. */
/* The difference between them indicates the extent to which color balance */
/* compensation needs to be performed. The more they point in the same direction, the stronger */
/* the color compensation needs to be. */



/* ------------------------------------------------------------------------------------ */
/* Get RGB color directions of the White and Black points and change average RGB colors */
/* ------------------------------------------------------------------------------------ */
/* We determine the direction of the average RGB values of the black and white points */
/* This returns two angles, and we examine their difference, */
/* the less different they are, the more we replace the average RGB values */
/* with the perfect black and white values, thus more color balance correction occurs during contrast adjustment */
/* */
/* Since their directions encompass a 360° angle, and different colors are at every 60°, */
/* therefore, we consider a deviation greater than 60° as completely different. */
/* Meaning, only at a deviation less than 60° do we proportionally shift the target points of the black and */
/* white points towards the direction of perfect black and perfect white */
/* (i.e., in case of matching directions, full color correction occurs) */

    RGB_TO_HSL (bp_r_end, bp_g_end, bp_b_end, &H, &S, &L);
    temp1 = H;
    RGB_TO_HSL (wp_r_end, wp_g_end, wp_b_end, &H, &S, &L);
    temp2 = H;
    temp2 = temp2 - temp1;
    
    /* change value to positive */
    if (temp2 < 0)   { temp2 = 0 - temp2; }
    /* if angle of direction is larger than 180 degree, then take the smaller section of the circle,
       it means it'll always be less or equal than 180 deg */
    if (temp2 > 0.5) { temp2 = 1 - temp2; }
    /* if the angle is greater than 60 degrees, that means the colors are totally different,
       so I check the amount of difference only on this 1/6 interval from 0 to 60 degrees. */
    temp2 = temp2 * 6;
    if (temp2 > 1){ temp2 = 1; }

    /* raise the value (0..1) of angle difference to 3rd power to make color balance a bit more aggressive */
    temp2 = temp2 * temp2 * temp2;

    /* With this, we have the value of the direction difference in a [0..1] interval, */
    /* where 0 shows complete agreement */
    /* now I only examine a sixth of the 'circle' and */
    /* create a value from it in the [0..1] interval, */
    /* so that I can use it to multiply the brightness of the black point, */
    /* meaning if the directions match, it goes to perfect black */
    /* and the same applies to the white point */
    if (temp2 < 1) {
        RGB_TO_HSL (bp_r_end, bp_g_end, bp_b_end, &H, &S, &L);
        L = L * temp2;
        HSL_TO_RGB (H, S, L, &bp_r_end, &bp_g_end, &bp_b_end);
        RGB_TO_HSL (wp_r_end, wp_g_end, wp_b_end, &H, &S, &L);
        L = 1 - (1 - L) * temp2;
        HSL_TO_RGB (H, S, L, &wp_r_end, &wp_g_end, &wp_b_end);
    }
    wp_end = (wp_r_end + wp_g_end + wp_b_end) / 3;
    bp_end = (bp_r_end + bp_g_end + bp_b_end) / 3;



/* ------------------------------------------------------------------------------ */
/* Convert original Histogram using White and Black point values */
/* ------------------------------------------------------------------------------ */
/* Creating a modified histogram from the original histogram based on */
/* the black and white point values */

/* Not re-analyzing the entire image, but just the original histogram, */
/* because this way only 256 values need to be processed instead of the total number of points in the image. */
/* This histogram represents the state of the image after automatic contrast adjustment */

/* This histogram will be used to determine the center of gravity for the gamma correction */

    #ifdef __OPENMP__
    #pragma omp parallel for private(temp2, temp3, cc, N) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){

	temp2 = (double)(i1) / 255;

        /* Calculation of contrast by pulling towards zero relative to bp and wp is as follows */
        /*temp2 = bp + ((temp2 - bp) * (1 - bp) / (wp - bp)); */
        /*temp2 = 1 - (1 - temp2) * 1 / (1 - bp); */

        /* Calculation of contrast with full range pull is as follows */
        /*temp2 = temp2 * wp_end / wp; */
        /* ---> bp = bp * wp_end / wp; */
        /* ---> bp_end = bp_end * wp_end / wp; */
        /*temp2 = 1 - (1 - temp2) * (1 - bp_end * wp_end / wp) / (1 - bp * wp_end / wp); */

        /* Calculation of contrast by pulling from bp and wp towards bp_end and wp_end is as follows */
        if ((temp2 > bp_end) && (wp > bp_end)) {
		temp2 = bp_end + (temp2 - bp_end) * (wp_end - bp_end) / (wp - bp_end); }
        /* here bp black point also needs to be multiplied for the following calculation, */
        /* because we are stretching the scale to the right towards white relative to bp_end */
        /* and therefore bp shifts */
	temp3 = 0;
        if (wp > bp_end) {
		temp3 = bp_end + (bp - bp_end) * (wp_end - bp_end) / (wp - bp_end); }
        if ((temp2 < wp_end) && (wp_end != temp3)) {
		temp2 = wp_end - (wp_end - temp2) * (wp_end - bp_end) / (wp_end - temp3); }
        /* in 'if' statements everywhere I check to prevent division by zero */

        if (temp2 > 1){ temp2 = 1; }
        if (temp2 < 0){ temp2 = 0; }
        cc = (long)(temp2 * 255);
#ifdef __OPENMP__
	N = omp_get_thread_num();
#else
	N = 0;
#endif
        hist2n[cc + N*256] += hist1[i1];
    }

    /* this is a replacement code for arrays for multi processing */
    #ifdef __OPENMP__
    #pragma omp parallel for private(i1, i2) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){
        for (i2=0; i2<max_threads2; i2++){ hist2[i1] += hist2n[i1 + i2*256]; } }



/* ------------------------------------------------------------------------------ */
/* Gamma value calculating */
/* ------------------------------------------------------------------------------ */
/* Determining the gamma center of gravity based on the second histogram, which shows */
/* the situation after contrast adjustment */
/* The final result indicates how much the image needs to be lightened or darkened */
/* to balance the overall brightness of the image */

/* The gamma center of gravity is the value in the histogram that indicates */
/* an equal number of image points on both its left and right sides */
/* (i.e., half of all the points in the image) */

/* This is obtained by starting to read the histogram values from one side inward */
/* and summing up the values we get */
/* If this sum reaches or exceeds half of the total number of points in the image */
/* we stop and the current index of the array gives us the appropriate value */
/* for the gamma center of gravity */

/*
    gamma_weight_mid_all = 0;
    for (i1=0;   i1<256; i1++){ gamma_weight_mid_all  += hist2[i1]; }
    i3 = 0;
    flag1 = 0;
    gamma_weight_mid = 0;
    for (i1=0; i1<256; i1++){
        i3 = i3 + hist2[i1];
        if (flag1 == 0){
            if (i3 > gamma_weight_mid_all / 2){
                flag1 = 1;
                gamma_weight_mid = i1;
            }
        }
    }
    gamma_weight_mid = gamma_weight_low / 255;
    gamma_mid = 1;
    gamma_mid = log(gamma_interval_mid) / log(gamma_weight_mid);
    if (gamma_mid < (1/gamma_max)){ gamma_mid = (1/gamma_max); }
    if (gamma_mid > gamma_max){ gamma_mid = gamma_max; }
*/
    /* Convert CONTRAST Histogram using MID GAMMA VALUE */ /*
    for (i1=0; i1<256; i1++){
        temp2 = i1;
        temp2 = pow(temp2 / 255, gamma_mid) * 255;
        i2 = (long)temp2;
        if (i2 > 255){ i2 = 255; }
        if (i2 < 0){ i2 = 0; }
        hist2b[i2] = hist2b[i2] + hist2[i1];
    }
*/

    /* Determining the histogram weight */
    gamma_weight_low_all = 0;
    gamma_weight_high_all = 0;

    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:gamma_weight_low_all) num_threads(max_threads2)
    #endif
    for (i1=0;   i1<128; i1++){ gamma_weight_low_all  += hist2[i1]; }
    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:gamma_weight_high_all) num_threads(max_threads2)
    #endif
    for (i1=128; i1<256; i1++){ gamma_weight_high_all += hist2[i1]; }

    /* Determining the histogram center of gravity */
    i3 = 0;
    flag1 = 0;
    gamma_weight_low = 0;
    for (i1=0; i1<128; i1++){
        i3 = i3 + hist2[i1];
        if (flag1 == 0){
            if (i3 > gamma_weight_low_all / 2){
                flag1 = 1;
                gamma_weight_low = i1;
            }
        }
    }
    i3 = 0;
    flag1 = 0;
    gamma_weight_high = 0;
    for (i1=128; i1<256; i1++){
        i3 = i3 + hist2[i1];
        if (flag1 == 0){
            if (i3 > gamma_weight_high_all / 2){
                flag1 = 1;
                gamma_weight_high = i1;
            }
        }
    }
    gamma_weight_low = gamma_weight_low / 255;
    gamma_weight_high = gamma_weight_high / 255;

    /* Determining the necessity of shifting the center of gravity */
    gamma_low = 1;
    gamma_high = 1;

    /* we only shift gamma in one direction, */
    /* meaning we only lighten if necessary, but never darken */
    if (gamma_weight_low < gamma_interval_low){
        gamma_low = log(gamma_interval_low) / log(gamma_weight_low);
    }
    if (gamma_weight_high > gamma_interval_high){
        gamma_high = log(gamma_interval_high) / log(gamma_weight_high);
    }

    if (gamma_low < (1/gamma_max)){ gamma_low = (1/gamma_max); }
    if (gamma_low > gamma_max){ gamma_low = gamma_max; }
    if (gamma_high < (1/gamma_max)){ gamma_high = (1/gamma_max); }
    if (gamma_high > gamma_max){ gamma_high = gamma_max; }



/* ------------------------------------------------------------------------------ */
/* Recalculate the RGB values by setting the CONTRAST, COLOR BALANCE and GAMMA */
/* ------------------------------------------------------------------------------ */
/* Recalculating the colors of the image (contrast, color balance, and gamma correction) */
/* All points of the image are recalculated and rewritten into the buffer */

/* In gamma adjustment, the RGB color of the image is adjusted not separately for each color channel, */
/* but together according to their brightness */
/* The gamma is adjusted such that the center is shifted to 128 */
/* (so if the value is less, it increases, and if more, it decreases) */
/* and this proportionally stretches the mid-colors towards 128 */

    #ifdef __OPENMP__
    #pragma omp parallel for private(x, y, addr, addr2, col_r, col_g, col_b, col_r2, col_g2, col_b2, temp2, temp3, cc, H, S, L, N) num_threads(max_threads2)
    #endif
    for (y=0; y<=bh-1; y++){
        addr2 = y * bw * 3 + y * addr_offset;
        for (x=0; x<=bw-1; x++){
            addr = addr2 + x * 3;

            /* apply changes ONLY on selected area of the image */
            if ((apply_on_selection == 0) || ((apply_on_selection) &&
                (x >= x1 &&
                x <= x2 &&
                y >= y1 &&
                y <= y2))) {

            col_r2 = image_buffer[addr + 0];
            col_g2 = image_buffer[addr + 1];
            col_b2 = image_buffer[addr + 2];

            col_r2 /= 255;
            col_g2 /= 255;
            col_b2 /= 255;


            /* CONTRAST SHIFT AND COLOR BALANCE */
            /* Calculation of contrast by pulling from bp_end and wp_end towards bp and wp as follows */
            if ((col_r2 > bp_r_end) && (wp_r > bp_r_end)) {
                col_r2 = bp_r_end + (col_r2 - bp_r_end) * (wp_r_end - bp_r_end) / (wp_r - bp_r_end); }
            if ((col_g2 > bp_g_end) && (wp_g > bp_g_end)) {
                col_g2 = bp_g_end + (col_g2 - bp_g_end) * (wp_g_end - bp_g_end) / (wp_g - bp_g_end); }
            if ((col_b2 > bp_b_end) && (wp_b > bp_b_end)) {
                col_b2 = bp_b_end + (col_b2 - bp_b_end) * (wp_b_end - bp_b_end) / (wp_b - bp_b_end); }

            /* here the bp black point also needs to be multiplied for the next calculation, */
            /* because we are stretching the scale to the right towards white relative to bp_end */
            /* and therefore bp shifts */
	    temp3 = 0;
            if (wp_r > bp_r_end) {
                temp3 = bp_r_end + (bp_r - bp_r_end) * (wp_r_end - bp_r_end) / (wp_r - bp_r_end); }
            if ((col_r2 < wp_r_end) && (wp_r_end != temp3)) {
                col_r2 = wp_r_end - (wp_r_end - col_r2) * (wp_r_end - bp_r_end) / (wp_r_end - temp3); }

	    temp3 = 0;
            if (wp_g > bp_g_end) {
                temp3 = bp_g_end + (bp_g - bp_g_end) * (wp_g_end - bp_g_end) / (wp_g - bp_g_end); }
            if ((col_g2 < wp_g_end) && (wp_g_end != temp3)) {
                col_g2 = wp_g_end - (wp_g_end - col_g2) * (wp_g_end - bp_g_end) / (wp_g_end - temp3); }

	    temp3 = 0;
            if (wp_b > bp_b_end) {
                temp3 = bp_b_end + (bp_b - bp_b_end) * (wp_b_end - bp_b_end) / (wp_b - bp_b_end); }
            if ((col_b2 < wp_b_end) && (wp_b_end != temp3)) {
                col_b2 = wp_b_end - (wp_b_end - col_b2) * (wp_b_end - bp_b_end) / (wp_b_end - temp3); }

            /* boundary check and rewrite */
            if (col_r2 > 1){ col_r2 = 1; }
            if (col_g2 > 1){ col_g2 = 1; }
            if (col_b2 > 1){ col_b2 = 1; }
            if (col_r2 < 0){ col_r2 = 0; }
            if (col_g2 < 0){ col_g2 = 0; }
            if (col_b2 < 0){ col_b2 = 0; }


            /* GAMMA CORRECTION */
            /* note: this method of gamma lifting gives a more colorful result */
            /* than the one where colors are shifted individually, not proportionally */
            temp2 = (col_r2 + col_g2 + col_b2) / 3;
            cc = (long)(temp2 * 255);

            /*temp2 = pow(temp2 / 255, gamma_mid) * 255; */
            if (temp2 <= 0.5){ temp3 = pow(temp2 * 2, gamma_low) / 2; }
            else{ temp3 = pow((temp2 - 0.5) * 2, gamma_high) / 2 + 0.5; }

            /*temp3 = pow(temp2 / 255, gamma_exp) * 255; */
            if (temp2 < temp3){
                if (temp2 < 1){
                    if (temp2 > 0){
                        col_r2 = ((1 - (1 - col_r2) * (1 - temp3)
                            / (1 - temp2)) + (col_r2 * temp3 / temp2)) / 2;
                        col_g2 = ((1 - (1 - col_g2) * (1 - temp3)
                            / (1 - temp2)) + (col_g2 * temp3 / temp2)) / 2;
                        col_b2 = ((1 - (1 - col_b2) * (1 - temp3)
                            / (1 - temp2)) + (col_b2 * temp3 / temp2)) / 2;
                    }
                }
            }
            else{
                if (temp2 > 0){
                    col_r2 = col_r2 * temp3 / temp2;
                    col_g2 = col_g2 * temp3 / temp2;
                    col_b2 = col_b2 * temp3 / temp2;
                }
            }

            /* boundary check and rewrite */
            if (col_r2 > 1){ col_r2 = 1; }
            if (col_g2 > 1){ col_g2 = 1; }
            if (col_b2 > 1){ col_b2 = 1; }
            if (col_r2 < 0){ col_r2 = 0; }
            if (col_g2 < 0){ col_g2 = 0; }
            if (col_b2 < 0){ col_b2 = 0; }


            col_r = (long)(col_r2 * 255);
            col_g = (long)(col_g2 * 255);
            col_b = (long)(col_b2 * 255);
            image_buffer[addr + 0] = col_r & 0xff;
            image_buffer[addr + 1] = col_g & 0xff;
            image_buffer[addr + 2] = col_b & 0xff;

            /* Build up Histogram for SATURATION */
            if (x >= x1 &&
                x <= x2 &&
                y >= y1 &&
                y <= y2) {
                    double H, S, L;
                    RGB_TO_HSL (col_r2, col_g2, col_b2, &H, &S, &L);
                    /* Here, I consider the 'S' value smaller as we move away from the center of the gray line, */
                    /* because optically it appears less colorful, */
                    /* and here I am examining optically */
                    if (L > 0.5) { L = 1 - L; }
                    S = S * L * 2;
                    /* If S = 0, meaning the color is gray, then I don't add it */
                    /* to the saturation histogram for obvious reasons */
                    if (S > 0) {

                        cc = (long)(S * 255);
#ifdef __OPENMP__
			N = omp_get_thread_num();
#else
			N = 0;
#endif
                        hist_saturn[cc + N*256]++;
                    }
            }

            }
        }
    }

    /* this is a replacement code for arrays for multi processing */
    #ifdef __OPENMP__
    #pragma omp parallel for private(i1, i2) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){
        for (i2=0; i2<max_threads2; i2++){ hist_satur[i1] += hist_saturn[i1 + i2*256]; } }


/* ------------------------------------------------------------------------------ */
/* Recalculate the RGB values by setting the SATURATION */
/* ------------------------------------------------------------------------------ */

    /* Determining the average value of the histogram */
    hist_satur_avg = 0;

    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:hist_satur_avg) num_threads(max_threads2)
    #endif
    for (i1=255; i1>=0; i1--){
        hist_satur_avg += hist_satur[i1];
    }
    hist_satur_avg = hist_satur_avg / 255;
    /* Establishing a cut-off limit, which is 10% of the average, has proven effective in contrast adjustment,
    with this value the histogram self-adjusts properly,
    and the edges are cut off at the appropriate size */
    temp1 = hist_satur_avg * 0.1;

    /* Searching the edge of the histogram to lift the colors */
    i3 = 0;
    flag1 = 0;
    for (i1=255; i1>=0; i1--){
        if (flag1 == 0){
            if (hist_satur[i1] >= temp1){
                flag1 = 1;
                i3 = i1;
            }
        }
    }
    hist_satur_low = (double)(i3) / 255;

    /* Boundary check */
    hist_satur_ok = 1;
    if (hist_satur_low > satur_max){ hist_satur_low = satur_max; }
    if (hist_satur_low > 0){ hist_satur_ok = log(satur_max) / log(hist_satur_low); }

  /* Run saturation recalculation only if necessary */
  if (hist_satur_ok != 1){

    #ifdef __OPENMP__
    #pragma omp parallel for private(x, y, addr, addr2, col_r, col_g, col_b, col_r2, col_g2, col_b2, H, S, L, cc, N) num_threads(max_threads2)
    #endif
    for (y=0; y<=bh-1; y++){
        addr2 = y * bw * 3 + y * addr_offset;
        for (x=0; x<=bw-1; x++){
            double S_new;
            double S_diff;

            addr = addr2 + x * 3;

            /* Apply changes ONLY on selected area of the image */
            if ((apply_on_selection == 0) || ((apply_on_selection) &&
                (x >= x1 &&
                x <= x2 &&
                y >= y1 &&
                y <= y2))) {

            col_r2 = image_buffer[addr + 0];
            col_g2 = image_buffer[addr + 1];
            col_b2 = image_buffer[addr + 2];

            col_r2 = col_r2 / 255;
            col_g2 = col_g2 / 255;
            col_b2 = col_b2 / 255;

            RGB_TO_HSL (col_r2, col_g2, col_b2, &H, &S, &L);

            /* Power raising of color saturation similar to gamma,
               meaning I raise the original saturation value to a power,
               and multiply it with twice the distance from 0.5,
               so that those at the 0.5 point increase maximally, while those farther away
               increasingly less, and those at 0 and 1 not at all,
               
               in other words, I adjust it exponentially, but such that those closer to the center
               adjust more, while those farther away less -
               this is necessary because with simple power raising the curve is too drastic
               and the lower parts of the histogram jump too much,
               this way, however, the change will be appropriate, no matter how big it is.
            */

            if (hist_satur_ok != 1){
            	S_new = pow(S, hist_satur_ok);
            	S_diff = 1 - (fabs(0.5 - S) * 2);
            	S = (S_new - S) * S_diff + S;
            }

            HSL_TO_RGB (H, S, L, &col_r2, &col_g2, &col_b2);

            col_r = (long)(col_r2 * 255);
            col_g = (long)(col_g2 * 255);
            col_b = (long)(col_b2 * 255);
            image_buffer[addr + 0] = col_r & 0xff;
            image_buffer[addr + 1] = col_g & 0xff;
            image_buffer[addr + 2] = col_b & 0xff;

            /* Additionally, creating the final Histogram */
            cc = (long)((double)(col_r + col_g + col_b) / 3);

            /*hist3[cc]++; */
        /* this is a replacement code for arrays for multi processing */
#ifdef __OPENMP__
	    N = omp_get_thread_num();
#else
	    N = 0;
#endif
            hist3n[cc + N*256]++;

            }
        }
    }

    /* this is a replacement code for arrays for multi processing */
    #ifdef __OPENMP__
    #pragma omp parallel for private(i1, i2) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){
        for (i2=0; i2<max_threads2; i2++){ hist3[i1] += hist3n[i1 + i2*256]; } }

  }
  /* just copy the former histogram with no change because there was no saturation processing */
  else {
    for (i1=0; i1<256; i1++){ hist3[i1] = hist2[i1]; }
  }


/* The test part only works if there is a math library, because here I use trigonometric functions */

/* ------------------------------------------------------------------------------ */
/* TEST: Show Histograms by Drawing them on Image */
/* ------------------------------------------------------------------------------ */
/* */
    /* Testing is only possible with a normal RGB array */
    if ((format_flag == 0) && (test_flag == 1)) {
        /* The starting value of the max is 1, to avoid division by zero */
        long hist1_max = 1;
        long hist2_max = 1;
        long hist3_max = 1;
        long histS_max = 1;

        /* Be aware! Running the routine multiple times on the same image */
        /* produces unexpected results if the histograms are also drawn, */
        /* because then the test image is also included in the calculation */
        for (i1=0; i1<256; i1++){ if (hist1_max < hist1[i1]) { hist1_max = hist1[i1]; } }
        for (i1=0; i1<256; i1++){ if (hist2_max < hist2[i1]) { hist2_max = hist2[i1]; } }
        for (i1=0; i1<256; i1++){ if (hist3_max < hist3[i1]) { hist3_max = hist3[i1]; } }
        for (i1=0; i1<256; i1++){ if (histS_max < hist_satur[i1]) { histS_max = hist_satur[i1]; } }

        for (y=0; y<=bh-1; y++){
            addr2 = y * bw * 3 + y * addr_offset;
            for (x=0; x<=bw-1; x++){
                addr = addr2 + x * 3;
                xm = x;
                ym = y;
                if (format_flag == 1) { ym = bh - 1 - ym; }
                if (format_flag == 2) { ym = bh - 1 - ym; }

                /* Drawing a frame around the histograms with 1-pixel width */
                if (((xm == 256) && (ym <= 601)) || ((ym == 601) && (xm <= 256))) {
                        col = color_black;
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                }

                /* Drawing the Histograms */
                if ((xm >= 0) && (xm <= 255)) {
                    double rad, outline, pi;

                    /* 1. histogram: ORIGINAL IMAGE STATE WITH CONTRAST ADJUSTMENT */
                    if (ym <= 99) {
                        i2 = hist1[xm] * 99 / hist1_max;
                        xma1 = (long)(bp * 255);
                        xmb1 = (long)(wp * 255);
                        xma2 = (long)(bp_end * 255);
                        xmb2 = (long)(wp_end * 255);
                        /* show the average limit with a horizontal line */
                        if (y == (99 - (long)(hist_cut_limit) * 99 / hist_max)){
                            col = color_yellow;
                        }
                        else{
                            if (99 - ym >= i2) {
                                if (((xm > xma2) && (xm < xma1)) || ((xm < xmb2) && (xm > xmb1))) {
                                    col = color_brown;
                                }
                                else {
                                    col = color_black;
                                }
                            }
                            else {
                                col = color_green;
                            }
                        }
                        if (xm == xma1) { col = color_white; }
                        if (xm == xmb1) { col = color_white; }
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                    }

                    /* 2. histogram: AFTER CONTRAST WITH GAMMA ADJUSTMENT */
                    if ((ym >= 100) && (ym <= 199)) {
                        i2 = hist2[xm] * 99 / hist2_max;
                        xma1 = (long)(gamma_weight_low * 255);
                        xmb1 = (long)(gamma_weight_high * 255);
                        xma2 = (long)(gamma_interval_low * 255);
                        xmb2 = (long)(gamma_interval_high * 255);
                        if (199 - ym >= i2) {
                            if ((xm < xma2) && (xm > xma1))
                            {
                                col = color_brown;
                            }
                            else{
                                /* I mark it with red instead of brown, indicating that this would be a backward gamma adjustment,
                                but this is not calculated */
                                if ((xm > xma2) && (xm < xma1))
                                {
                                    col = color_red;
                                }
                                else{
                                    col = color_black;
                                }
                            }
                        }
                        else {
                            col = color_green;
                        }
                        /* Lower gamma center of gravity */
                        if (xm == xma1) {
                            col = color_white;
                        }
                        /* Upper gamma center of gravity */
/*                        if (xm == (long)(gamma_weight_high * 255)) {
                            col = color_white;
                        }*/
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                    }

                    /* 3. histogram: WITH SATURATION ADJUSTMENT */
                    if ((ym >= 200) && (ym <= 299)) {
                        i2 = hist_satur[xm] * 99 / histS_max;
                        xma1 = (long)((hist_satur_low)* 255);
                        xma2 = (long)(satur_max * 255);
                        if (299 - ym >= i2) {
                            if (((xm < xma2) && (xm > xma1)) || ((xm > xma2) && (xm < xma1)))
                            {
                                col = color_brown;
                            }
                            else{
                                col = color_black;
                            }
                        }
                        else {
                            col = color_blue;
                        }
                        if (xm == xma1) { col = color_white; }
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                    }

                    /* 4. histogram: AFTER CONTRAST, GAMMA, AND SATURATION ADJUSTMENT (FINAL) */
                    if ((ym >= 300) && (ym <= 399)) {
                        i2 = hist3[xm] * 99 / hist3_max;
                        if (399 - ym >= i2) {
                            col = color_black;
                        }
                        else {
                            col = color_green;
                        }
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                    }

                    /* 5.: Drawing the color direction of black and white points in a circle */
                    /* Constant for the radius of the circle */
                    rad = 100;
                    /* Constant for the thickness of the circle line */
                    outline = 10;
                    /* Value of pi */
                    pi = 3.1415926535897932;
                    if ((ym >= 400) && (ym <= 400+rad*2)) {
                        double rr, rx, ry;

                        rx = xm - rad;
                        ry = ym - 400 - rad;
                        rr = sqrt(rx*rx + ry*ry);

                        /* DRAWING THE CIRCLE ARC */
                        if ((rr <= rad) && (rr >= rad - outline)) {
                                H = 0;
                                if ((rx >= 0) && (ry <  0)) { H = asin(rx / rr) / pi / 2 + 0.00; }
                                if ((rx >= 0) && (ry >= 0)) { H = asin(ry / rr) / pi / 2 + 0.25; }
                                if ((rx <  0) && (ry >= 0)) { H = asin(-rx / rr) / pi / 2 + 0.50; }
                                if ((rx <  0) && (ry <  0)) { H = asin(-ry / rr) / pi / 2 + 0.75; }
                                S = 1;
                                L = 0.5;
                                L = rad - rr;
                                if (L > outline / 2) { L = outline - L; }
                                L = L / outline;
                                HSL_TO_RGB(H, S, L, &col_r2, &col_g2, &col_b2);
                                col_r = (long)(col_r2 * 255);
                                col_g = (long)(col_g2 * 255);
                                col_b = (long)(col_b2 * 255);
                                image_buffer[addr + 0] = col_r & 0xff;
                                image_buffer[addr + 1] = col_g & 0xff;
                                image_buffer[addr + 2] = col_b & 0xff;
                        }
                        else {
                            double px, py, lx, ly, d;

                            col = color_black;
                            px = 0; py = 0;
                            lx = 0; ly = 0;

                            /* DRAWING THE COLOR DIRECTION LINE OF THE BLACK POINT */
                            /* The line does not appear in the circle for a gray point */
                            RGB_TO_HSL(bp_r, bp_g, bp_b, &H, &S, &L);

                            if ((H >= 0.00) && (H < 0.25)) { lx =  sin((H - 0.00) * 2 * pi) * rr; ly =  cos((H - 0.00) * 2 * pi) * rr; }
                            if ((H >= 0.25) && (H < 0.50)) { lx =  cos((H - 0.25) * 2 * pi) * rr; ly = -sin((H - 0.25) * 2 * pi) * rr; }
                            if ((H >= 0.50) && (H < 0.75)) { lx = -sin((H - 0.50) * 2 * pi) * rr; ly = -cos((H - 0.50) * 2 * pi) * rr; }
                            if ((H >= 0.75) && (H < 1.00)) { lx = -cos((H - 0.75) * 2 * pi) * rr; ly =  sin((H - 0.75) * 2 * pi) * rr; }
                            px = lx * S;
                            py = ly * S;

                            d = ly * rx + lx * ry;
                            d = d / sqrt(lx * lx + ly * ly);
                            if (d < 0) { d = 0 - d; }

                            if ((d < 1) && (rr < rad) && (sgn(rx) == sgn(px)) && (sgn(ry) == -sgn(py))) {
                                if (rr < rad * S) { col = color_white; }
                                else { col = color_gray; }
                            }

                            /* DRAWING THE COLOR DIRECTION LINE OF THE WHITE POINT */
                            RGB_TO_HSL(wp_r, wp_g, wp_b, &H, &S, &L);

                            if ((H >= 0.00) && (H < 0.25)) { lx =  sin((H - 0.00) * 2 * pi) * rr; ly =  cos((H - 0.00) * 2 * pi) * rr; }
                            if ((H >= 0.25) && (H < 0.50)) { lx =  cos((H - 0.25) * 2 * pi) * rr; ly = -sin((H - 0.25) * 2 * pi) * rr; }
                            if ((H >= 0.50) && (H < 0.75)) { lx = -sin((H - 0.50) * 2 * pi) * rr; ly = -cos((H - 0.50) * 2 * pi) * rr; }
                            if ((H >= 0.75) && (H < 1.00)) { lx = -cos((H - 0.75) * 2 * pi) * rr; ly =  sin((H - 0.75) * 2 * pi) * rr; }
                            px = lx * S;
                            py = ly * S;

                            d = ly * rx + lx * ry;
                            d = d / sqrt(lx * lx + ly * ly);
                            if (d < 0) { d = 0 - d; }

                            if ((d < 1) && (rr < rad) && (sgn(rx) == sgn(px)) && (sgn(ry) == -sgn(py))) {
                                if (rr < rad * S) { col = color_white; }
                                else { col = color_gray; }
                            }

                            image_buffer[addr + 2] = (col >> 0)  & 0xff;
                            image_buffer[addr + 1] = (col >> 8)  & 0xff;
                            image_buffer[addr + 0] = (col >> 16) & 0xff;

                        }
                    }
                }
            }
        }
    }


exit:

    if (hist_saturn){ free(hist_saturn); }
    if (hist1n)     { free(hist1n);      }
    if (hist2n)     { free(hist2n);      }
    if (hist3n)     { free(hist3n);      }
    if (col_r3n)    { free(col_r3n);     }
    if (col_g3n)    { free(col_g3n);     }
    if (col_b3n)    { free(col_b3n);     }

}



/* ------------------------------ */
/* ----- MAIN AARGB ENTRIES ----- */
/* ------------------------------ */

void AARGB_NORMAL(
    unsigned char *image_buffer,
    int image_width,
    int image_height)
{
    AARGB_MAIN(
        image_buffer,
        image_width,
        image_height,
        0,
        0,
        image_width-1,
        image_height-1,
        0,
        0,
        0);
}

void AARGB_NORMAL_SEL(
    unsigned char *image_buffer,
    int image_width,
    int image_height,
    int x1,
    int y1,
    int x2,
    int y2,
    int apply_on_selection)
{
    AARGB_MAIN(
        image_buffer,
        image_width,
        image_height,
        x1,
        y1,
        x2,
        y2,
        0,
        apply_on_selection,
        0);
}

void AARGB_DIB(
    unsigned char *image_buffer,
    int image_width,
    int image_height)
{
    AARGB_MAIN(
        image_buffer,
        image_width,
        image_height,
        0,
        0,
        image_width-1,
        image_height-1,
        1,
        0,
        0);
}

void AARGB_DIB_SEL(
    unsigned char *image_buffer,
    int image_width,
    int image_height,
    int x1,
    int y1,
    int x2,
    int y2,
    int apply_on_selection)
{
    AARGB_MAIN(
        image_buffer,
        image_width,
        image_height,
        x1,
        y1,
        x2,
        y2,
        1,
        apply_on_selection,
        0);
}

void AARGB_BMP(
    unsigned char *image_buffer)
{
	long f_bm;
	long f_bitcount;
	long f_compressed;
	long f_offs;
	long f_width;
	long f_height;

	f_bm = 0;
	f_bm += image_buffer[0] << 0;
	f_bm += image_buffer[1] << 8;

	f_bitcount = 0;
	f_bitcount += image_buffer[28] << 0;
	f_bitcount += image_buffer[29] << 8;

	f_compressed = 0;
	f_compressed += image_buffer[30] << 0;
	f_compressed += image_buffer[31] << 8;
	f_compressed += image_buffer[32] << 16;
	f_compressed += image_buffer[33] << 24;

	if ((f_bm == 0x00004d42) && (f_bitcount == 24) && (f_compressed == 0)) {

		f_offs = 0;
		f_offs += image_buffer[10] << 0;
		f_offs += image_buffer[11] << 8;
		f_offs += image_buffer[12] << 16;
		f_offs += image_buffer[13] << 24;

		f_width = 0;
		f_width += image_buffer[18] << 0;
		f_width += image_buffer[19] << 8;
		f_width += image_buffer[20] << 16;
		f_width += image_buffer[21] << 24;

		f_height = 0;
		f_height += image_buffer[22] << 0;
		f_height += image_buffer[23] << 8;
		f_height += image_buffer[24] << 16;
		f_height += image_buffer[25] << 24;

		AARGB_MAIN(
			image_buffer + f_offs,
			f_width,
			f_height,
			0,
			0,
			f_width-1,
			f_height-1,
            		2,
			0,
            		0);
	}
}

void AARGB_BMP_SEL(
    unsigned char *image_buffer,
    int x1,
    int y1,
    int x2,
    int y2,
    int apply_on_selection)
{
	long f_bm;
	long f_bitcount;
	long f_compressed;
	long f_offs;
	long f_width;
	long f_height;

	f_bm = 0;
	f_bm += image_buffer[0] << 0;
	f_bm += image_buffer[1] << 8;

	f_bitcount = 0;
	f_bitcount += image_buffer[28] << 0;
	f_bitcount += image_buffer[29] << 8;

	f_compressed = 0;
	f_compressed += image_buffer[30] << 0;
	f_compressed += image_buffer[31] << 8;
	f_compressed += image_buffer[32] << 16;
	f_compressed += image_buffer[33] << 24;

	if ((f_bm == 0x00004d42) && (f_bitcount == 24) && (f_compressed == 0)) {

		f_offs = 0;
		f_offs += image_buffer[10] << 0;
		f_offs += image_buffer[11] << 8;
		f_offs += image_buffer[12] << 16;
		f_offs += image_buffer[13] << 24;

		f_width = 0;
		f_width += image_buffer[18] << 0;
		f_width += image_buffer[19] << 8;
		f_width += image_buffer[20] << 16;
		f_width += image_buffer[21] << 24;

		f_height = 0;
		f_height += image_buffer[22] << 0;
		f_height += image_buffer[23] << 8;
		f_height += image_buffer[24] << 16;
		f_height += image_buffer[25] << 24;

		AARGB_MAIN(
			image_buffer + f_offs,
			f_width,
			f_height,
			x1,
			y1,
			x2,
			y2,
            		2,
			apply_on_selection,
            		0);
	}
}




/*
int main (){
	return 0;
}
*/
