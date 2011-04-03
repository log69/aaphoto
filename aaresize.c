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


int RESIZE(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    unsigned long new_width,
    unsigned long new_height,
    int alpha_flag)
{
	unsigned int startx, endx;
	unsigned int starty, endy;
	unsigned int i, j, x, y;

	unsigned char *bitmap_buffer2;
	/* if there is an Alpha channel beside the RGB too */
	/* then i need to allocate more memory for Alpha channel */
	/* not 3x but 5x more than the number of pixels in image */
	if (alpha_flag) { alpha_flag = 2; }
	/*allocate memory */
	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
	bitmap_buffer2 = malloc(new_width * new_height * (3 + alpha_flag) * sizeof (*bitmap_buffer2));
	if (bitmap_buffer2 == 0) return 1;

	/* print info */ STRING_PRINTV("resizing image\n");

	#ifdef __OPENMP__
	#pragma omp parallel for private(i, j, x, y, startx, starty, endx, endy) num_threads(max_threads)
	#endif
	for (j=0; j<new_height; j++){

		for (i=0; i<new_width; i++){

			long addr = 0;
			long rgb_avg_r = 0;
			long rgb_avg_g = 0;
			long rgb_avg_b = 0;
			long rgb_avg_a = 0;
			long rgb_counter = 0;

			startx = *image_width  * i     / new_width;
			endx   = *image_width  * (i+1) / new_width;
			starty = *image_height * j     / new_height;
			endy   = *image_height * (j+1) / new_height;

			for (y=starty; y<endy; y++){
				for (x=startx; x<endx; x++){

					addr = (x + y * *image_width);
					rgb_avg_r += (*image_buffer)[3 * addr + 0];
					rgb_avg_g += (*image_buffer)[3 * addr + 1];
					rgb_avg_b += (*image_buffer)[3 * addr + 2];
					if (alpha_flag) {
						rgb_avg_a += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr + 0];
					}
					rgb_counter++;

				}
			}
			if (rgb_counter > 0){
				rgb_avg_r = rgb_avg_r / rgb_counter;
				rgb_avg_g = rgb_avg_g / rgb_counter;
				rgb_avg_b = rgb_avg_b / rgb_counter;
				if (alpha_flag) { rgb_avg_a = rgb_avg_a / rgb_counter; }
			}

			addr = (i + j * new_width);
			bitmap_buffer2[3 * addr + 0] = rgb_avg_r;
			bitmap_buffer2[3 * addr + 1] = rgb_avg_g;
			bitmap_buffer2[3 * addr + 2] = rgb_avg_b;
			if (alpha_flag) {
				bitmap_buffer2[new_width * new_height * 3 + addr + 0] = rgb_avg_a;
			}

		}
	}

	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(*image_buffer);
	*image_buffer = bitmap_buffer2;
	*image_width = new_width;
	*image_height = new_height;

	return 0;
}



int ROTATE90(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    int alpha_flag)
{
	unsigned char *bitmap_buffer2;
	unsigned int x, y;

	unsigned long new_width;
    	unsigned long new_height;
	new_width = *image_width;
	new_height = *image_height;

	/*allocate memory */
	if (alpha_flag) { alpha_flag = 2; }
	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
	bitmap_buffer2 = malloc(new_width * new_height * (3 + alpha_flag) * sizeof (*bitmap_buffer2));
	if (bitmap_buffer2 == 0) return 1;

	/* print info */ STRING_PRINTV("rotating image with 90\n");

	#ifdef __OPENMP__
	#pragma omp parallel for private(x, y) num_threads(max_threads)
	#endif
	for (y=0; y<new_height; y++){

		for (x=0; x<new_width; x++){

			long addr = 0;
			long rgb_avg_r = 0;
			long rgb_avg_g = 0;
			long rgb_avg_b = 0;
			long rgb_avg_a = 0;

			addr = (x + y * new_width);
			rgb_avg_r += (*image_buffer)[3 * addr + 0];
			rgb_avg_g += (*image_buffer)[3 * addr + 1];
			rgb_avg_b += (*image_buffer)[3 * addr + 2];
			if (alpha_flag){ rgb_avg_a += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr + 0]; }

			addr = ((new_height - 1 - y) + x * new_height);
			bitmap_buffer2[3 * addr + 0] = rgb_avg_r;
			bitmap_buffer2[3 * addr + 1] = rgb_avg_g;
			bitmap_buffer2[3 * addr + 2] = rgb_avg_b;
			if (alpha_flag){ bitmap_buffer2[new_width * new_height * 3 + addr + 2] = rgb_avg_a; }
		}
	}

	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(*image_buffer);
	*image_buffer = bitmap_buffer2;
	*image_width = new_height;
	*image_height = new_width;

	return 0;
}



int ROTATE180(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    int alpha_flag)
{
	int x, y, x_max, y_max, y_odd;

	unsigned long new_width;
	unsigned long new_height;
	new_width = *image_width;
	new_height = *image_height;

	/* print info */ STRING_PRINTV("rotating image with 180\n");

	/* is the number of lines an odd number? */
	y_odd = new_height & 1;
	y_max = new_height/2 + y_odd;

	#ifdef __OPENMP__
	#pragma omp parallel for private(x, y, x_max) num_threads(max_threads)
	#endif
	for (y=0; y<y_max; y++){

		x_max = new_width;
		/* if the height is an odd number, then the flipping should go only to the middle of the width */
		/* in the middle line */
		if ((y_odd) && (y == y_max-1)){ x_max = new_width/2; }

		for (x=0; x<x_max; x++){

			long addr1 = 0;
			long addr2 = 0;
			long rgb_avg_r1 = 0;
			long rgb_avg_g1 = 0;
			long rgb_avg_b1 = 0;
			long rgb_avg_a1 = 0;
			long rgb_avg_r2 = 0;
			long rgb_avg_g2 = 0;
			long rgb_avg_b2 = 0;
			long rgb_avg_a2 = 0;

			addr1 = (x + y * new_width);
			addr2 = ((new_width - 1 - x) + (new_height - 1 - y) * new_width);

			rgb_avg_r1 += (*image_buffer)[3 * addr1 + 0];
			rgb_avg_g1 += (*image_buffer)[3 * addr1 + 1];
			rgb_avg_b1 += (*image_buffer)[3 * addr1 + 2];
			if (alpha_flag){ rgb_avg_a1 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0]; }
			rgb_avg_r2 += (*image_buffer)[3 * addr2 + 0];
			rgb_avg_g2 += (*image_buffer)[3 * addr2 + 1];
			rgb_avg_b2 += (*image_buffer)[3 * addr2 + 2];
			if (alpha_flag){ rgb_avg_a2 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0]; }

			(*image_buffer)[3 * addr1 + 0] = rgb_avg_r2;
			(*image_buffer)[3 * addr1 + 1] = rgb_avg_g2;
			(*image_buffer)[3 * addr1 + 2] = rgb_avg_b2;
			if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0] = rgb_avg_a2; }
			(*image_buffer)[3 * addr2 + 0] = rgb_avg_r1;
			(*image_buffer)[3 * addr2 + 1] = rgb_avg_g1;
			(*image_buffer)[3 * addr2 + 2] = rgb_avg_b1;
			if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0] = rgb_avg_a1; }
		}
	}

	return 0;
}



int ROTATE270(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    int alpha_flag)
{
	unsigned char *bitmap_buffer2;
	unsigned int x, y;

	unsigned long new_width;
    	unsigned long new_height;
	new_width = *image_width;
	new_height = *image_height;

	/*allocate memory */
	if (alpha_flag) { alpha_flag = 2; }
	/* print info */ STRING_PRINTV("allocating memory for uncompressed bitmap image\n");
	bitmap_buffer2 = malloc(new_width * new_height * (3 + alpha_flag) * sizeof (*bitmap_buffer2));
	if (bitmap_buffer2 == 0) return 1;

	/* print info */ STRING_PRINTV("rotating image with 270\n");

	#ifdef __OPENMP__
	#pragma omp parallel for private(x, y) num_threads(max_threads)
	#endif
	for (y=0; y<new_height; y++){

		for (x=0; x<new_width; x++){

			long addr = 0;
			long rgb_avg_r = 0;
			long rgb_avg_g = 0;
			long rgb_avg_b = 0;
			long rgb_avg_a = 0;

			addr = (x + y * new_width);
			rgb_avg_r += (*image_buffer)[3 * addr + 0];
			rgb_avg_g += (*image_buffer)[3 * addr + 1];
			rgb_avg_b += (*image_buffer)[3 * addr + 2];
			if (alpha_flag){ rgb_avg_a += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr + 0]; }

			addr = (y + (new_width - 1 - x) * new_height);
			bitmap_buffer2[3 * addr + 0] = rgb_avg_r;
			bitmap_buffer2[3 * addr + 1] = rgb_avg_g;
			bitmap_buffer2[3 * addr + 2] = rgb_avg_b;
			if (alpha_flag){ bitmap_buffer2[(*image_width) * (*image_height) * 3 + addr + 0] = rgb_avg_a; }
		}
	}

	/* print info */ STRING_PRINTV("freeing uncompressed bitmap buffer\n");
	free(*image_buffer);
	*image_buffer = bitmap_buffer2;
	*image_width = new_height;
	*image_height = new_width;

	return 0;
}



int FLIPX(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    int alpha_flag)
{
	unsigned int x, y;

	unsigned long new_width;
    	unsigned long new_height;
	new_width = *image_width;
	new_height = *image_height;

	/* print info */ STRING_PRINTV("flipping image horizontally\n");

	#ifdef __OPENMP__
	#pragma omp parallel for private(x, y) num_threads(max_threads)
	#endif
	for (y=0; y<new_height; y++){

		for (x=0; x<(new_width/2); x++){

			long addr1 = 0;
			long addr2 = 0;
			long rgb_avg_r1 = 0;
			long rgb_avg_g1 = 0;
			long rgb_avg_b1 = 0;
			long rgb_avg_a1 = 0;
			long rgb_avg_r2 = 0;
			long rgb_avg_g2 = 0;
			long rgb_avg_b2 = 0;
			long rgb_avg_a2 = 0;

			addr1 = (x + y * new_width);
			addr2 = ((new_width - 1 - x) + y * new_width);

			rgb_avg_r1 += (*image_buffer)[3 * addr1 + 0];
			rgb_avg_g1 += (*image_buffer)[3 * addr1 + 1];
			rgb_avg_b1 += (*image_buffer)[3 * addr1 + 2];
                        if (alpha_flag){ rgb_avg_a1 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0]; }
			rgb_avg_r2 += (*image_buffer)[3 * addr2 + 0];
			rgb_avg_g2 += (*image_buffer)[3 * addr2 + 1];
			rgb_avg_b2 += (*image_buffer)[3 * addr2 + 2];
                        if (alpha_flag){ rgb_avg_a2 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0]; }

			(*image_buffer)[3 * addr1 + 0] = rgb_avg_r2;
			(*image_buffer)[3 * addr1 + 1] = rgb_avg_g2;
			(*image_buffer)[3 * addr1 + 2] = rgb_avg_b2;
			if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0] = rgb_avg_a2; }
			(*image_buffer)[3 * addr2 + 0] = rgb_avg_r1;
			(*image_buffer)[3 * addr2 + 1] = rgb_avg_g1;
			(*image_buffer)[3 * addr2 + 2] = rgb_avg_b1;
                        if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0] = rgb_avg_a1; }
		}
	}

	return 0;
}



int FLIPY(
    unsigned char **image_buffer,
    unsigned long *image_width,
    unsigned long *image_height,
    int alpha_flag)
{
	unsigned int x, y;

	unsigned long new_width;
    	unsigned long new_height;
	new_width = *image_width;
	new_height = *image_height;

	/* print info */ STRING_PRINTV("flipping image vertically\n");

	#ifdef __OPENMP__
	#pragma omp parallel for private(x, y) num_threads(max_threads)
	#endif
	for (y=0; y<(new_height/2); y++){

		for (x=0; x<new_width; x++){

			long addr1 = 0;
			long addr2 = 0;
			long rgb_avg_r1 = 0;
			long rgb_avg_g1 = 0;
			long rgb_avg_b1 = 0;
			long rgb_avg_a1 = 0;
			long rgb_avg_r2 = 0;
			long rgb_avg_g2 = 0;
			long rgb_avg_b2 = 0;
			long rgb_avg_a2 = 0;

			addr1 = (x + y * new_width);
			addr2 = (x + (new_height - 1 - y) * new_width);

			rgb_avg_r1 += (*image_buffer)[3 * addr1 + 0];
			rgb_avg_g1 += (*image_buffer)[3 * addr1 + 1];
			rgb_avg_b1 += (*image_buffer)[3 * addr1 + 2];
                        if (alpha_flag){ rgb_avg_a1 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0]; }
			rgb_avg_r2 += (*image_buffer)[3 * addr2 + 0];
			rgb_avg_g2 += (*image_buffer)[3 * addr2 + 1];
			rgb_avg_b2 += (*image_buffer)[3 * addr2 + 2];
                        if (alpha_flag){ rgb_avg_a2 += (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0]; }

			(*image_buffer)[3 * addr1 + 0] = rgb_avg_r2;
			(*image_buffer)[3 * addr1 + 1] = rgb_avg_g2;
			(*image_buffer)[3 * addr1 + 2] = rgb_avg_b2;
			if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr1 + 0] = rgb_avg_a2; }
			(*image_buffer)[3 * addr2 + 0] = rgb_avg_r1;
			(*image_buffer)[3 * addr2 + 1] = rgb_avg_g1;
			(*image_buffer)[3 * addr2 + 2] = rgb_avg_b1;
                        if (alpha_flag){ (*image_buffer)[(*image_width) * (*image_height) * 3 + addr2 + 0] = rgb_avg_a1; }
		}
	}

	return 0;
}

