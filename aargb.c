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
/* ----------- András Horváth (C) 2006-2013 ---------- */
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
/* Kezdõ értékek és konstansok megadása a kontraszt és gamma mûveletekhez */

/* Kontraszt konstans megadása: megadja, hogy mekkora lehet az automata */
/* kontraszt állítás maximum értéke (0...1-ig terjedhet az értéke), */
/* alapértelmezett = 0.1 */
/* Gamma konstans megadása: megadja, hogy mekkora lehet az automatikus */
/* gamma állítás maximum értéke (1...10-ig ajánlott), */
/* alapértelmezett = 1.5 */

/* Gamma állításhoz a maximum elõfordulás értékének megadása (occur_max), */
/* amely megadja, hogy a gamma érték számításánál ha ennél kisebb az */
/* elõfordulása a nagy mértékben elõforduló színeknek, akkor kihagyjuk õket a */
/* számításból, vagyis a nagy mértékben elõforduló színeket azért nem vesszük */
/* be a számításba, mert inkább a részletek látszódjanak jól, */
/* vagyis a részletekhez legyen inkább kiszámolva a megfelelõ fényerõ */


/* ------------------------------ */
/* Konstans értékek beállítása */
/* ------------------------------ */
    /* maximális kontraszt mértéke (0..1) */
    cont_max = 0.066666;

    /* maximális gamma korrekció mértéke (1..10) */
    gamma_max = 1.5;
    gamma_interval_low = 0.333;
    gamma_interval_high = 1;

    /* maximális színtelítettség limit (0..1) */
    satur_max = 0.333;

    bw = image_width;
    bh = image_height;

    /* Kijelölt területhez a koordináták határértékeinek vizsgálata */
    /* és megfelelõ beállítása */
    /* A kijelölés célja, hogy egy meghatározott képrészletet szeretnénk jól */
    /* láthatóvá tenni (nem pedig az egészet arányaiban) */

    /* a szélesség konvertálása abszolút koordinátává, jelenleg nem él */
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
    /* A DIB formátum 4-byte-os igazításához az eltolás értékének kiszámítása */
    /* Normál tömbnél erre nincs szükség, ekkor a format_flag értéke = 0 */
    /* egyébként a DIB formátumú tömb függõlegesen fordított sorokat tartalmaz, */
    /* és a sorvégek 4 byte-tal vannak igazítva */
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
    /* Kijelölés koordinátáinak felcserélése, ha nem jó sorrendben adták meg, */
    /* vagyis a bal felsõ sarok az x1 és y1, a jobb alsó pedig az x2 és y2 */
    if (x1 > x2) { i1 = x1; x1 = x2; x2 = i1; }
    if (y1 > y2) { i1 = y1; y1 = y2; y2 = i1; }


/* ------------------------------------------------------------------------------ */
/* Create Histogram and average RGB colors for the Image */
/* ------------------------------------------------------------------------------ */
/* Hisztogram generálása a kép színeibõl, ami a kép feényerõ eloszlását adja meg */
/* plusz a színegyensúly beállításához az egy fényerejû színek átlag RGB */
/* értékeinek letárolása, megelõlegezve egy késõbbi rutin munkáját */

/* A kép minden egyes pontjának átlag fényereje (szürkéje) bekerül egy */
/* 256 elemû tömbbe, ahol a fényerejük értéke az azonos indexû elem értékét */
/* 1-gyel növeli */

/* Így ez a hisztogram tömb pontos leírást ad a feketétõl a fehérig terjedõ */
/* színskálájáról a képnek */

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
                /* szürke hisztogramm */
                hist1n[cc + N*256]++;
                /* Átlag RGB értékek letárolása a fekete és fehér pont */
                /* átlag RGB-jének későbbi megállapításához */
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
/* Az automata kontraszt beállításához a Fekete és Fehér pont megállapítása */

/* A hisztogram bal és jobb oldaláról elkezdem beolvasni a szürke értékek */
/* nagyságát és addig olvasom be, amíg az nem nagyobb egy elõre meghatározott */
/* értéknél, ekkor megkapom a fekete és fehér pont helyzetét */

/* Ez az elõre meghatározott érték a kontraszt konstans és a hisztogram */
/* összes tömbeleme átlagának a szorzata */
/* Az átlag szorzat egyensúlyt teremt a határérték elérésénél, mert egyébként */
/* ha ez az érték mondjuk a hisztogram maximum értéke lenne, akkor */
/* drasztikus kontraszt túlállítás jellemezné a funkciót */

/* Vagyis ez az érték azt adja meg, hogy a kép átlag fényerejének hány */
/* százaléka az az érték, amely a továbbiakban megadja, hogy az ekkora */
/* százalék alatt található feketék és fehérek lesznek kihúzva a határig (le lesznek vágva)*/

    hist_min = bw * bh;
    hist_max = 0;
    hist_avg = 0;

    /* Átlag fényerõ kiszámítása: ezt úgy kapom meg, hogy összeadom a histogramm
    összes oszlopát és osztom 256-al (oszlopok száma), vagyis matematikai átlaga.
    Ezt még leosztom egy konstanssal (10%-ára alapból) úgy, hogy nullánál kisebbel szorzok,
    az így kapott értéket nevezem itt limit-nek.

    (a maximum és minimum érték számítása csak tesztelési céllal él)
    */
    for (i1=0; i1<256; i1++){
        temp1 = hist1[i1];
        if (hist_min > (long)temp1){ hist_min = (long)temp1; }
        if (hist_max < (long)temp1){ hist_max = (long)temp1; }
        hist_avg = hist_avg + temp1;
    }
    /* histogram teljes összege */
/*    hist_sum = hist_avg; */

    /* histogram matematikai átlag értéke */
    hist_avg = hist_avg / 256;

/*
    hist_min_test = hist_min;
    hist_max_test = hist_max;
    hist_avg_test = hist_avg;
*/

    /* ez lesz itt a limit */
    temp1 = hist_avg * cont_max;
    hist_cut_limit = temp1;



    /* fehér és feketepont keresése 0 (nulla) maximumig */

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
    Az automatikus kontraszt állító algoritmus néha olyan bemenetet is kaphat,
    ahol a részletek a képen túl kis mennyiségben vannak jelen, és nagy terjedelmű
    egy színű részek a jellemzőek - ekkor az volt a jelenség, hogy az algoritmus
    túl nagyot vág le a hisztogram széléből, és az ezen a képen fontosabb kis mennyiségű
    részlet veszik el.
    
    Ezért beillesztettem az algoritmusba egy önszabályozó mechanizmust, amelynél ha
    felmerül a fenti eset, akkor jó esély van rá, hogy a részletek maradnak meg inkább
    a kontraszt állítás után.
    
    Ez úgy történik, hogy megnézem, mennyi oszlop esik a limit érték alá, és ennek függvényében
    tovább csökkentem ezt a limit-et, mivel minél több rész esik az átlag alá (oszlopok száma a hisztogramban),
    ez azt jelenti hogy annál több a képen az olyan érdekes részlet, amely levágásra kerülne,
    ezért a sok és meghatározó nagy üres részek kevésbé kellene hogy számítsanak.

    A limit értéket az alábbi módon húzom tovább lefelé: újra fogom kalkulálni az átlagot
    és ezt a limit értéket is úgy, hogy a kalkuláció során az ez fölé eső oszlop értékek
    kevésbé számítsanak, ha minél több az előzőkben a limit érték alatti részek oszlopainak száma,
    és ha ezen oszlopok magasságának összegei is minél kisebbek,
    mert így tudjuk meg, hogy egyre kisebb területen van az értékes részlet.
    
    Ez így egy önfékező folyamatot eredményez a nagy üres területek átlaga okozta túl nagy limit értékhez,
    és ezzel a túl nagy kontraszt levágáshoz, ahol is pont a részlet veszik el.

    Ezek után újra kalkulálom ezt a limit értéket az előző módon,
    és ezzel lesz tovább kalkulálva a histogram szélek levágása.
    
    Hogy mennyire számítson a kis kép részlet, azt is analóg módon akarom eldönteni,
    tehát úgy, hogy az eredeti limit érték módosulása ne szakaszos módon történjen,
    viszont mivel úgy gondolom, hogy olyan görbére van szükségem, amely a 0-1 intervallumon
    a feléig nagyon kicsit emelkedik, majd innét drasztikusabban,
    ezért az 5. hatványt szorozva 3-mal találtam a megfelelőnek.
    
    f(x) = x^5*3
    Wolphram Aplha link a szemléltetéshez:
    http://www.wolframalpha.com/input/?i=x^5+*+3+from+0+to+1
    
    Ennek mentén, ha a kis és nagy részlet egyensúly felborul, akkor
    a 0.5 től felfelé kezd a kontraszt számítás egyre kevésbé drasztikusba átmenni,
    és ezzel érem el, hogy a kisebb értékek lehetőleg minél kevésbé,
    míg a nagyobb értékek egyre jobban folyásolják be ezt.
    -------------------------------------------------------------------------------------- */


    /* azért csak a fekete és fehér pont közötti szakaszt vizsgálom
       (ez a hisztogramm megmaradó szélessége a bal és jobb oldali levágás után)
       mert a tervezett végleges eredményen akarom vizsgálni a részlet mennyiségének arányát
       a teljes terjedelemhez képest (amennyi maradna belőle).
    */

    hist_cut_columns = 0;
    hist_cut_weight = 0;
    for (i1=bp; i1<wp; i1++){
        if (hist1[i1] < temp1){
            hist_cut_columns++;
            hist_cut_weight += hist1[i1];
        }
    }


    /* temp1 mutatja az eredeti levágandó limit értéket,
       hist_cut_columns mutatja a limit alatti oszlopok számát (ezt nevezem limit alatti részletnek),
       hist_cut_weight mutatja ezen oszlopok értékének összegét (súlyát),

       temp2-t pedig úgy állítom be, hogy minél nagyobb súly oszlik el kevesebb oszlop számon
       (vagyis minél kevesebb oszlop van a limit alatt és ezeknek a súlya minél nagyobb),
       úgy ennek is annál nagyobb lesz az értéke - vagyis ha temp2-nek nagyobb az értéke,
       az azt jelenti hogy annál több értékes részlet van a limit alatt.
       
       azért nevezem a limit alatti oszlopokat értékesebb részletnek, mert a teljes hisztogramm átlagot
       elhúzzák felfelé a nagy egyszínű részek, amelyek kevés oszlopok nagy súllyal,
       vagyis ezek nyilván nagyobb terjedelmű üres részek - tehát ezek maguk a "nem" részletek,
       míg ezen átlag alattiakat veszem a részletnek.
       
       mivel a limit alatti oszlopokat nézzük, ezért leosztva magával a limit értékkel, egy
       0..1 intervallumos arány értéket kapok.
    */
    if ((hist_cut_columns == 0) || (temp1 == 0)){ temp2 = 0; }
    else { temp2 = (double)(hist_cut_weight) / hist_cut_columns / temp1; }


    /* itt temp3 értékét úgy határozom meg, hogy az eredeti limit alatti oszlopok számát
       osztom a histogram középső (levágás utáni) megmaradt szélességével (wp-bp)
       (de csak a megmaradandó szakaszon, ezért az eredmény 0..1 közötti lesz),

       vagyis minél nagyobb rész esik le bal és jobb oldalt, annál kisebb értékkel osztunk,
       és ha a limit alatti oszlopok száma egyre több, akkor ezt minél kisebb értékkel osztva
       annál nagyobb számot kapunk, ezért annál drasztikusabbnak vehetjük az eredetileg alkalmazandó kontrasztot,
       és ezért ezt az értéket az 5. hatványra emelve és szorozva 3-mal - olyan értéket eredményez,
       mely 50% fölött egyre nagyobb értéket ad vissza, és itt minél nagyobb az érték,
       annál jobban csökkentem az eredeti tervezett kontraszt (levágás) mértékét.
    */
    temp3 = (double)(hist_cut_columns) / (wp-bp);

    /* itt a lineáris értéket hatványra emeléssel görbítem, hogy a kisebb értékek kevésbé,
    míg a nagyobb értékek egyre jobban befolyásolják az eredményt */
    temp3 = temp3 * temp3 * temp3 * temp3 * temp3 * 3;
    if (temp3 > 1){ temp3 = 1; }

    /* itt temp1 (eredeti limit érték) alatt keletkezett temp2 (új limit érték)
       limit vonalat visszahúzom felfelé az eredeti temp1 felé a temp3-as görbe alapján.
       
       vagyis a drasztikus kontrasztot lecsökkentettem, majd vissza engedem a görbe alapján
       (amelynél kb. 50% után görbül drasztikusan).
    */
    temp3 = temp1 - ((temp1 - temp2) * temp3);
    /* itt nem engedem hogy az új csökkentett limit érték az eredeti 10% alá menjen,
       ez csupán egy alsó korlát a kontraszt csökkentés mértékéhez */
    if (temp3 < temp1 * 0.1){ temp1 = temp1 * 0.1; }
    else{ temp1 = temp3; }

    hist_cut_limit = temp1;

/* ------------------------------------------------------------ */


    bp = 255;
    flag1 = 0;

    /* histogram jobb oldaláról lépkedek és vágom majd le addig az értékig,
    amely még kisebb mint az átlag * cont_max (átlag 10 %-a) */
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

    /* Határértékek beállítása és korrekciója */
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
/* A színegyensúly beállításához az átlag RGB értékek kiszámítása */
/* a fekete és fehér pont értéke alapján */

/* Itt keletkezik egy átlag RGB érték a fekete és fehér pontokhoz egyaránt */
/* Ez az érték azt adja meg, hogy az automatikus kontraszt állításakor */
/* minden egyes szín milyen irányba torzul lefelé (fekete pont RGB átlaga) */
/* és felfelé (fehér pont RGB átlaga) */

/* A fehér pont feletti összes szín átlagának RGB-je lesz a viszonyítási pont */
/* a fehér írányába való torzításhoz, */
/* vagyis ez lesz kihúzva a tökéletes fehérbe */

/* Ez valóságban a kép színegyensúlyát állítja be megfelelõen úgy, hogy a */
/* levágandó mértékû fehérek színátlaga lesz a tökéletes fehér, */
/* ezért ha ezek átlaga nem tökéletes fehér, akkor az ettõl eltérõ nagyságot */
/* minden színnél arányosan el kell tolni a tökéletes fehér irányába, */
/* ugyanez a fekete estében */

    bp_r = 0;
    bp_g = 0;
    bp_b = 0;
    wp_r = 0;
    wp_g = 0;
    wp_b = 0;

    i3 = 0;
    /* fekete pont alatti összes szín RGB átlagának kiszámítása */
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
    /* fehér pont feletti összes szín RGB átlagának kiszámítása */
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

    /* skálázás 255-ről a [0..1] intervallumra */
    bp_r = bp_r / 255;
    bp_g = bp_g / 255;
    bp_b = bp_b / 255;
    wp_r = wp_r / 255;
    wp_g = wp_g / 255;
    wp_b = wp_b / 255;

    /* A kapott átlag RGB érték fényerejének visszaállítása a fehér pont szintjére. */
    /* Mivel ugye nem csak a fehér pont fényerejével azonos színeknek kalkuláltuk ki */
    /* az átlag színét, hanem az attól világosabb összes színnek, ezért a kapott */
    /* átlag szín fényereje nagyobb vagy egyenlő lesz, mint a kiindulási fehér pont */
    /* ezért a korrekcióhoz visszaállítjuk az RGB érték fényerejét */
    /* de az R, G és B komponensek arányának a megtartásával */
    /* és ugyanez a fekete pont esetében */
    RGB_TO_HSL (bp_r, bp_g, bp_b, &H, &S, &L);
    L = bp;
    HSL_TO_RGB (H, S, L, &bp_r, &bp_g, &bp_b);
    RGB_TO_HSL (wp_r, wp_g, wp_b, &H, &S, &L);
    L = wp;
    HSL_TO_RGB (H, S, L, &wp_r, &wp_g, &wp_b);

    /* A fekete és fehér pont célpontjának kiszámítása. */
    /* Ez mutatja meg, hogy a fekete és fehér pont átlag RGB-jét */
    /* hova kell húzni úgy, hogy az RGB kockában a két csúcsot */
    /* összekötő 'szürke' egyenessel párhuzamosan tolva a távolsága */
    /* a 'szürke' egyenestől és a színiránya megmaradjon, */
    /* de a lehető legsötétebb- vagy legvilágosabb legyen */
    /* */
    /* Másképpen fogalmazva eltoljuk a szürke egyenes mentén addig, */
    /* amíg az RGB kocka falába nem ütközünk (mindkét iránynál) */
    /* */
    /* Ez  annyiban változtatás az előzőkhöz képest, hogy a fekete pontot */
    /* most már nem a tökéletes feketébe húzzuk, hanem az annak megfelelő */
    /* olyan legsötétebb pontba, ahol a maximum a színtelítettség */
    /* és aminek színe megegyezik a fekete pontéval, */
    /* ezzel a rossz színegyensúlyt és nem megfelelő kontrasztot küszöbölöm ki. */

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
/* ---- RGB TÉR ---- */
/* ----------------- */
/* A teljes RGB teret egy szabályos 3D-s kocka foglalja magába, */
/* amelynek 1-1 éle jelenti a R, a G és a B tengelyt */
/* és egyik csúcsában található a tökéletes fehér szín, */
/* a másik (ezzel szemköti) csúcsában pedig a tökéletes fekete */
/* és ezt a két csúcsot összekötő egyenes tartalmazza a */
/* feketétől fehérig terjedő teljes szürke skálát. */
/* A színeltolás mértéke pedig nem más, mint az adott szín */
/* távolsága merőleges írányban szürke skála egyenesétől */
/* (amit egy 0 és 1 közötti érték jellemezhet, ahol a 0 */
/* azt jelenti, hogy a szín szürke, vagyis az egyenesen található) */
/* */
/* ------------------------- */
/* ---- RGB SZÍN IRÁNYA ---- */
/* ------------------------- */
/* A fekete és fehér pont átlag eltolási RGB értékét összehasonlítom, */
/* hogy megállapítsam, vajon megegyező irányban vannak-e eltolva, */
/* mivel ha nem jó a kép színegyensúlya, akkor feltételezem, hogy */
/* a kép összes színe a színegyensúly felborulását okozó tényező miatt */
/* megegyező irányban tolódik el. Ha nem megegyező irányba mutat */
/* az eltolásuk értéke, akkor feltételezem, hogy ez nem azért van, */
/* mert a színegyensúly felborult. Ekkor a kontraszt műveletnél nem */
/* alkalmazok színegyensúly kiegyenlítést (vagyis a kép színeinek */
/* a tökéletes fehér és tökéletes fekete irányába való RGB korrekcióját). */
/* */
/* Egy adott szín irányán az RGB kockában található pontjából kiinduló */
/* merőleges szakasz körülforgási szögét értem a szürke egyenesre nézve. */
/* Ennek értéke -180 és 180 fok közé kell hogy essen. */
/* Így a fekete és fehér pont átlag RGB színeinek iránya megad két szöget. */
/* Ennek különbsége adja meg, hogy milyen mértékkel kell színegyensúly */
/* kompenzációt végezni. Minél jobban egyírányba mutatnak, annál erősebb */
/* színkompenzáció szükséges. */



/* ------------------------------------------------------------------------------------ */
/* Get RGB color directions of the White and Black points and change average RGB colors */
/* ------------------------------------------------------------------------------------ */
/* Megállapítjuk a fekete és fehér pont átlag RGB értékeinek irányát */
/* Ez két szöget ad vissza, és ennek a különbségét vizsgáljuk, */
/* minél kevésbé eltérő, az átlag RGB értékeket annál jobban lecseréljük */
/* a tökéletes fekete és fehér értékre, így a kontraszt állításnál */
/* jobban keletkezik színegyensúly korrekció */
/* */
/* Mivel az irányuk egy 360˚-os szöget zár be, és az eltérő színek 60˚-onként vannak, */
/* ezért a 60˚-nál nagyobb eltérést teljesen különbözőnek vesszük. */
/* Vagyis a 60˚-nál kisebb eltérésnél toljuk csak el arányosan a fekete és */
/* fehér pont célpontját a tökéletes fekete és a tökéletes fehér írányába */
/* (vagyis egyező irány esetén teljes színkorrekció lép fel) */

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
    /* if the angle is greater then 60 degree, that means the colors are totally different,
       so i check the amount of difference on this 1/6 intervall only from 0 to 60 degrees. */
    temp2 = temp2 * 6;
    if (temp2 > 1){ temp2 = 1; }

    /* raise the value (0..1) of angle difference to 3th power to make color balance a bit more aggressive */
    temp2 = temp2 * temp2 * temp2;

    /* Ezzel megvan az iránykülönbség értéke egy [0..1] intervallumon, */
    /* ahol a 0 a teljes egyezést mutatja */
    /* most az egész 'kör' hatod részét vizsgálom csak és */
    /* abból alakítok ki egy értéket a [0..1] intervallumon, */
    /* hogy majd ezzel szorozni tudjam a fekete pont fényerejét, */
    /* vagyis ha egyeznek az írányok, akkor tökéletes feketébe megy el */
    /* ugyanez a fehér pont esetében */
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
/* Eredeti hisztogramból a fekete és fehér pont alapján megváltoztatott */
/* hisztogram létrehozása */

/* Nem az egész kép újraanalizálása, hanem csak az eredeti hisztogramé, */
/* mert így csak 256 értéket kell feldolgozni a kép összes pontjai számának */
/* helyett. Ez a hisztogram az automatikus kontraszt állítás utáni állapotát */
/* mutatja a képnek */

/* Ez a hisztogram lesz felhasználva a gamma súlypont megállapításához */

    #ifdef __OPENMP__
    #pragma omp parallel for private(temp2, temp3, cc, N) num_threads(max_threads2)
    #endif
    for (i1=0; i1<256; i1++){

	temp2 = (double)(i1) / 255;

        /* bp-től és wp-től viszonyított nullára húzással a kontraszt számolás az alábbi */
        /*temp2 = bp + ((temp2 - bp) * (1 - bp) / (wp - bp)); */
        /*temp2 = 1 - (1 - temp2) * 1 / (1 - bp); */

        /* teljes intervallumon számolt húzással a kontraszt számolás az alábbi */
        /*temp2 = temp2 * wp_end / wp; */
        /* ---> bp = bp * wp_end / wp; */
        /* ---> bp_end = bp_end * wp_end / wp; */
        /*temp2 = 1 - (1 - temp2) * (1 - bp_end * wp_end / wp) / (1 - bp * wp_end / wp); */

        /* bp_end-től és wp_end-től viszonyított bp-ből és wp-ből húzással a kontraszt számolás az alábbi */
        if ((temp2 > bp_end) && (wp > bp_end)) {
		temp2 = bp_end + (temp2 - bp_end) * (wp_end - bp_end) / (wp - bp_end); }
        /* itt a bp fekete pontot is fel kell szorozni a következő számításhoz, */
        /* mert a bp_end -től viszonyítva nyújtjuk a skálát jobbra a fehér irányába */
        /* és ezért elmászik a bp */
	temp3 = 0;
        if (wp > bp_end) {
		temp3 = bp_end + (bp - bp_end) * (wp_end - bp_end) / (wp - bp_end); }
        if ((temp2 < wp_end) && (wp_end != temp3)) {
		temp2 = wp_end - (wp_end - temp2) * (wp_end - bp_end) / (wp_end - temp3); }
        /*az 'if' utasításoknál mindenhol vizsgálom hogy ne lehessen nullával való osztás */

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
/* Gamma súlypont megállapítása a második hisztogram alapján, ami már az */
/* állított kontraszt utáni helyzetet mutatja */
/* A végeredmény azt adja meg, hogy mennyire kell világosítani, vagy éppen */
/* sötétíteni a képet, hogy az össz fényereje a képnek egyensúlyban legyen */

/* A gamma súlypont az az érték, amely a hisztogramban azt mutatja, */
/* hogy tõle balra és jobbra egyaránt egyforma számú képpont található */
/* (vagyis fele a kép összes pontjainak) */

/* Ezt úgy kapjuk meg, hogy elkezdjük olvasni a hisztogram értékeit az */
/* egyik oldalról befelé, és közben össze adjuk a kapott értékeket */
/* Ha ez az érték elérte vagy túllépte a kép összes pontjainak a számának */
/* felét, akkor megállunk és a tömb aktuális indexe adja meg */
/* a gamma súlypont megfelelõ értékét */

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

    /* Hisztogramm súlyának megállapítása */
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

    /* Hisztogramm súlypont megállapítása */
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

    /* Súlypont eltolás szükségességének megállapítása */
    gamma_low = 1;
    gamma_high = 1;

    /* gammát csak egyírányban toljuk el, */
    /* vagyis csak világosítunk ha szükséges, de soha sem sötétítünk */
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
/* Kép színeinek újrakalkulálása (kontraszt, színegyensúly és gamma korrekció) */
/* A kép összes pontja újrakalkulálódik és visszaíródik a pufferba */

/* A gamma állításnál a kép színének RGB-jét a gammához mérten nem külön */
/* színcsatornánként, hanem a fényerejükhöz mérten egyben vannnak állítva */
/* A gamma úgy állítódik, hogy a súlypont el van tolva 128-ba */
/* (tehát ha kisebb az értéke akkor növekszik, ha nagyobb, akkor meg csökken) */
/* és ez magával húzza nyújtásos módon arányosan a közép színeket */
/* a 128 irányába */


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
            /* bp_end-től és wp_end-től viszonyított bp-ből és wp-ből húzással a kontraszt számolás az alábbi */
            if ((col_r2 > bp_r_end) && (wp_r > bp_r_end)) {
                col_r2 = bp_r_end + (col_r2 - bp_r_end) * (wp_r_end - bp_r_end) / (wp_r - bp_r_end); }
            if ((col_g2 > bp_g_end) && (wp_g > bp_g_end)) {
                col_g2 = bp_g_end + (col_g2 - bp_g_end) * (wp_g_end - bp_g_end) / (wp_g - bp_g_end); }
            if ((col_b2 > bp_b_end) && (wp_b > bp_b_end)) {
                col_b2 = bp_b_end + (col_b2 - bp_b_end) * (wp_b_end - bp_b_end) / (wp_b - bp_b_end); }

            /* itt a bp fekete pontot is fel kell szorozni a következő számításhoz, */
            /* mert a bp_end -től viszonyítva nyújtjuk a skálát jobbra a fehér irányába */
            /* és ezért elmászik a bp */
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

            /* határérték ellenőrzés és visszaírás */
            if (col_r2 > 1){ col_r2 = 1; }
            if (col_g2 > 1){ col_g2 = 1; }
            if (col_b2 > 1){ col_b2 = 1; }
            if (col_r2 < 0){ col_r2 = 0; }
            if (col_g2 < 0){ col_g2 = 0; }
            if (col_b2 < 0){ col_b2 = 0; }


            /* GAMMA CORRECTION */
            /* megjegyzés: ez a gamma felhúzásos módszer színesebb végeredményt ad */
            /* mint amelyiknél külön - külön toljuk a színeket, nem arányosan */
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

            /* határérték ellenőrzés és visszaírás */
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
                    /* Az 'S' értékét a szürke egyenes közepétől távolodva kisebbnek veszem itt, */
                    /* mert az optikailag egyre kevésbé tűnik színesnek, */
                    /* és ennél a résznél optikailag vizsgálok */
                    if (L > 0.5) { L = 1 - L; }
                    S = S * L * 2;
                    /* Ha S = 0, vagyis szürke a szín, akkor nem adom hozzá */
                    /* a színtelítettség hisztogramjához értelemszerűen */
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
/* Recalculate the RGB values by setting the SATURAION */
/* ------------------------------------------------------------------------------ */

    /* Hisztogram átlag értékének megállapítása */
    hist_satur_avg = 0;

    #ifdef __OPENMP__
    #pragma omp parallel for reduction(+:hist_satur_avg) num_threads(max_threads2)
    #endif
    for (i1=255; i1>=0; i1--){
        hist_satur_avg += hist_satur[i1];
    }
    hist_satur_avg = hist_satur_avg / 255;
    /* levágási limit megállapítása, ez az átlag 10%-a bevált a kontrasztnál is,
    ezzel az értékkel megfelelően állítja be önmagát a hisztogram,
    és a megfelelő nagyságú szélek esnek le */
    temp1 = hist_satur_avg * 0.1;

    /* Hisztogram szélének keresése a színek felhúzásához */
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

    /* Határérték ellenőrzés */
    hist_satur_ok = 1;
    if (hist_satur_low > satur_max){ hist_satur_low = satur_max; }
    if (hist_satur_low > 0){ hist_satur_ok = log(satur_max) / log(hist_satur_low); }

  /* run saturation recalculation only if necessary */
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

            /* apply changes ONLY on selected area of the image */
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

            /* szín telítettségének hatványos emelése a gammához hasonlóan,
               vagyis hatványra emelem az eredeti színtelítettség értékét,
               és beszorzom a 0.5-től való távolságának kétszeresével,
               hogy a 0.5 pontban lévők maximálisan nővekedjenek, míg az ettől
               távolabbra lévők egyre kevesebb mértékben, a 0 és 1 helyen lévők pedig semennyire,
               
               másképpen fogalmazva, exponenciálisan állítom, de úgy, hogy a közepéhez
               közelebb lévők jobban állítódjanak, míg az ettől egyre távolabb esők kevésbé -
               erre azért van így szükség, mert a sima hatványra emelésnél túl drasztikus
               a görbe és az alsóbb részei a hisztogramnak is túl nagyot ugranak,
               így viszont megfelelő lesz a változás, mindegy mekkora az.
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

            /*Mellékesen a legvégső Hisztogram létrehozása is */
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


/* a test rész csak akkor működik, ha van math library, mert itt használok szögfüggvényeket */

/* ------------------------------------------------------------------------------ */
/* TEST: Show Histograms by Drawing them on Image */
/* ------------------------------------------------------------------------------ */
/* */
    /* csak normál RGB tömbnél élhet a tesztelés */
    if ((format_flag == 0) && (test_flag == 1)) {
        /* a max érték kezdõértéke 1, hogy 0-val való osztás ne fordulhasson elõ */
        long hist1_max = 1;
        long hist2_max = 1;
        long hist3_max = 1;
        long histS_max = 1;

        /* Figyelem! Többszöri lefuttatása a rutinnak ugyanazon a képen */
        /* nem várt eredményt produkál ha a hisztogrammok is ki vannak rajzolva, */
        /* mert akkor már a teszt képet is beleveszi a számításba */
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

                /* Keret rajzolása a hisztogrammok köré 1 pixel szélességben */
                if (((xm == 256) && (ym <= 601)) || ((ym == 601) && (xm <= 256))) {
                        col = color_black;
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                }

                /* Hisztogrammok kirajzolása */
                if ((xm >= 0) && (xm <= 255)) {
			double rad, outline, pi;

                    /* 1. hisztogramm: EREDETI KÉP ÁLLAPOTA KONTRASZT ÁLLÍTÁSSAL*/
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

                    /* 2. hisztogramm: KONTRASZT UTÁN GAMMA ÁLLÍTÁSSAL */
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
			        /* vörössel jelölöm barna helyett, hogy ez vissza irányú gamma állítás lenne,
			        ez viszont nem kerül számításra */
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
                        /* Gamma alsó súlypont */
                        if (xm == xma1) {
                            col = color_white;
                        }
                        /* Gamma felső súlypont */
/*                        if (xm == (long)(gamma_weight_high * 255)) {
                            col = color_white;
                        }*/
                        image_buffer[addr + 2] = (col >> 0)  & 0xff;
                        image_buffer[addr + 1] = (col >> 8)  & 0xff;
                        image_buffer[addr + 0] = (col >> 16) & 0xff;
                    }

                    /* 3. hisztogramm SZÍNTELÍTETTSÉG ÁLLÍTÁSSAL */
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

                    /* 4. hisztogramm: KONTRASZT, GAMMA és SZÍNTELÍTETTSÉG UTÁN (VÉGSŐ) */
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

                    /* 5. Fekete és fehér pont színirányának kirajzolása egy körbe */
                    /* konstans a kör sugarához */
                    rad = 100;
                    /* konstans a körvonal vastagságához */
                    outline = 10;
                    /* pi értéke */
                    pi = 3.1415926535897932;
                    if ((ym >= 400) && (ym <= 400+rad*2)) {
                        double rr, rx, ry;

                        rx = xm - rad;
                        ry = ym - 400 - rad;
                        rr = sqrt(rx*rx + ry*ry);

                        /* KÖRÍV MEGRAJZOLÁSA */
                        if ((rr <= rad) && (rr >= rad-outline)) {
                                H = 0;
                                if ((rx >= 0) && (ry <  0)) { H = asin( rx / rr) / pi / 2 + 0.00; }
                                if ((rx >= 0) && (ry >= 0)) { H = asin( ry / rr) / pi / 2 + 0.25; }
                                if ((rx <  0) && (ry >= 0)) { H = asin(-rx / rr) / pi / 2 + 0.50; }
                                if ((rx <  0) && (ry <  0)) { H = asin(-ry / rr) / pi / 2 + 0.75; }
                                S = 1;
                                L = 0.5;
                                L = rad - rr;
                                if (L > outline / 2) { L = outline - L; }
                                L = L / outline;
                                HSL_TO_RGB (H, S, L, &col_r2, &col_g2, &col_b2);
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

                            /* FEKETE PONT SZÍNIRÁNYA EGYENESÉNEK MEGHÚZÁSA */
                            /* Szürke pontnál nem jelenik meg egyenes a körben */
                            RGB_TO_HSL (bp_r, bp_g, bp_b, &H, &S, &L);

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

                            /* FEHÉR PONT SZÍNIRÁNYA EGYENESÉNEK MEGHÚZÁSA */
                            RGB_TO_HSL (wp_r, wp_g, wp_b, &H, &S, &L);

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

    if (hist_saturn){ free (hist_saturn); }
    if (hist1n)     { free (hist1n);      }
    if (hist2n)     { free (hist2n);      }
    if (hist3n)     { free (hist3n);      }
    if (col_r3n)    { free (col_r3n);     }
    if (col_g3n)    { free (col_g3n);     }
    if (col_b3n)    { free (col_b3n);     }

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
