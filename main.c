#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif
#define min3(x, y, z) min(min(x, y), z)

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#define max3(x, y, z) max(max(x, y), z)

#define CLAMP2BYTE(v) (((unsigned) (v)) < 255 ? (v) : (v < 0) ? 0 : 255)

unsigned int detect(uint8_t *pixel, int width, int height, int channels)
{
    int stride = width * channels;
    int last_col = width * channels - channels;
    int last_row = height * stride - stride;

    unsigned int sum = 0;
    for (int y = 0; y < height; y++) {
        int cur_row = stride * y;
        int next_row = min(cur_row + stride, last_row);
        uint8_t *next_scanline = pixel + next_row;
        uint8_t *cur_scanline = pixel + cur_row;
        for (int x = 0; x < width; x++) {
            int cur_col = x * channels;
            int next_col = min(cur_col + channels, last_col);
            uint8_t *c00 = cur_scanline + cur_col;
            uint8_t *c10 = cur_scanline + next_col;
            uint8_t *c01 = next_scanline + cur_col;
            uint8_t *c11 = next_scanline + next_col;
            int r_avg = ((c00[0] + c10[0] + c01[0] + c11[0])) >> 2;
            int g_avg = ((c00[1] + c10[1] + c01[1] + c11[1])) >> 2;
            int b_avg = ((c00[2] + c10[2] + c01[2] + c11[2])) >> 2;

            /* TODO: detect appropriate RGB values */
            if (r_avg >= 60 && g_avg >= 40 && b_avg >= 20 && r_avg >= b_avg &&
                (r_avg - g_avg) >= 10)
                if (max3(r_avg, g_avg, b_avg) - min3(r_avg, g_avg, b_avg) >= 10)
                    sum++;
        }
    }
    return sum;
}

void compute_offset(int *out, int len, int left, int right, int step)
{
    assert(out);
    assert((len >= 0) && (left >= 0) && (right >= 0));

    for (int x = -left; x < len + right; x++) {
        int pos = x;
        int len2 = 2 * len;
        if (pos < 0) {
            do {
                pos += len2;
            } while (pos < 0);
        } else if (pos >= len2) {
            do {
                pos -= len2;
            } while (pos >= len2);
        }
        if (pos >= len)
            pos = len2 - 1 - pos;
        out[x + left] = pos * step;
    }
}

void denoise(uint8_t *out,
             uint8_t *in,
             int width,
             int height,
             int channels,
             int radius)
{
    assert(in && out);
    assert(width > 0 && height > 0 && radius > 0);
    assert(channels == 1 || channels == 3);

    const int smoothing_level = 10; /* FIXME: should be adjustable */
    int window_size = (2 * radius + 1) * (2 * radius + 1);

    int *col_pow = malloc(width * channels * sizeof(int));
    int *col_val = malloc(width * channels * sizeof(int));
    int *row_pos = malloc((width + 2 * radius) * channels * sizeof(int));
    int *col_pos = malloc((height + 2 * radius) * channels * sizeof(int));

    int stride = width * channels;
    int smooth_table[256] = {0};
    float ii = 0.f;
    for (int i = 0; i <= 255; i++, ii -= 1.) {
        smooth_table[i] = (expf(ii * (1.0f / (smoothing_level * 255.0f))) +
                           (smoothing_level * (i + 1)) + 1) /
                          2;
        smooth_table[i] = max(smooth_table[i], 1);
    }

    compute_offset(row_pos, width, radius, radius, channels);
    compute_offset(col_pos, height, radius, radius, stride);

    int *row_off = row_pos + radius;
    int *col_off = col_pos + radius;
    for (int y = 0; y < height; y++) {
        uint8_t *scan_in_line = in + y * stride;
        uint8_t *scan_out_line = out + y * stride;
        if (y == 0) {
            for (int x = 0; x < stride; x += channels) {
                int col_sum[3] = {0};
                int col_sum_pow[3] = {0};
                for (int z = -radius; z <= radius; z++) {
                    uint8_t *sample = in + col_off[z] + x;
                    for (int c = 0; c < channels; ++c) {
                        col_sum[c] += sample[c];
                        col_sum_pow[c] += sample[c] * sample[c];
                    }
                }
                for (int c = 0; c < channels; ++c) {
                    col_val[x + c] = col_sum[c];
                    col_pow[x + c] = col_sum_pow[c];
                }
            }
        } else {
            uint8_t *last_col = in + col_off[y - radius - 1];
            uint8_t *next_col = in + col_off[y + radius];
            for (int x = 0; x < stride; x += channels) {
                for (int c = 0; c < channels; ++c) {
                    col_val[x + c] -= last_col[x + c] - next_col[x + c];
                    col_pow[x + c] -= last_col[x + c] * last_col[x + c] -
                                      next_col[x + c] * next_col[x + c];
                }
            }
        }

        int prev_sum[3] = {0}, prev_sum_pow[3] = {0};
        for (int z = -radius; z <= radius; z++) {
            int index = row_off[z];
            for (int c = 0; c < channels; ++c) {
                prev_sum[c] += col_val[index + c];
                prev_sum_pow[c] += col_pow[index + c];
            }
        }

        for (int c = 0; c < channels; ++c) {
            int mean = prev_sum[c] / window_size;
            int diff = mean - scan_in_line[c];
            int edge = CLAMP2BYTE(diff);
            int masked_edge =
                (edge * scan_in_line[c] + (256 - edge) * mean) >> 8;
            int var = (prev_sum_pow[c] - mean * prev_sum[c]) / window_size;
            int out = masked_edge -
                      diff * var / (var + smooth_table[scan_in_line[c]]);
            scan_out_line[c] = CLAMP2BYTE(out);
        }

        scan_in_line += channels, scan_out_line += channels;
        for (int x = 1; x < width; x++) {
            int last_row = row_off[x - radius - 1];
            int next_row = row_off[x + radius];
            for (int c = 0; c < channels; ++c) {
                prev_sum[c] -= col_val[last_row + c] - col_val[next_row + c];
                prev_sum_pow[c] = prev_sum_pow[c] - col_pow[last_row + c] +
                                  col_pow[next_row + c];
                int mean = prev_sum[c] / window_size;
                int diff = mean - scan_in_line[c];
                int edge = CLAMP2BYTE(diff);
                int masked_edge =
                    (edge * scan_in_line[c] + (256 - edge) * mean) >> 8;
                int var = (prev_sum_pow[c] - mean * prev_sum[c]) / window_size;
                int out = masked_edge -
                          diff * var / (var + smooth_table[scan_in_line[c]]);
                scan_out_line[c] = CLAMP2BYTE(out);
            }

            scan_in_line += channels, scan_out_line += channels;
        }
    }

    free(col_pow);
    free(col_val);
    free(row_pos);
    free(col_pos);
}

static void die(char *msg)
{
    fprintf(stderr, "Fatal: %s\n", msg);
    exit(-1);
}

int main(int argc, char **argv)
{
    if (argc < 2)
        return -1;

    int width = 0, height = 0, channels = 0;
    uint8_t *in = stbi_load(argv[1], &width, &height, &channels, 0);
    if (!in)
        die("Fail to load input file");

    int dimension = width * height;
    uint8_t *out = malloc(dimension * channels);
    if (!out)
        die("Out of memory");

    /* Separation between skin and non-skin pixels */
    float rate = detect(in, width, height, channels) / (float) dimension * 100;

    /* Perform edge detection, resulting in an edge map for further denoise */
    denoise(out, in, width, height, channels, min(width, height) / rate + 1);

    if (!stbi_write_jpg("out.jpg", width, height, channels, out, 100))
        die("Fail to generate");

    free(out);
    free(in);
    return 0;
}
