#include "utils.h"

void mpi_print(const char *c, const uint64_t *a, int len)
{
    int i;

    printf("%s", c);
    for (i = len - 1; i > 0; i--)
        printf("%016lX", a[i]);
    printf("%016lX\n", a[0]);
}

void mpi_conv_64to52(uint64_t *r, const uint64_t *a, int rlen, int alen)
{
    int i, j, shr_pos, shl_pos;
    uint64_t word, temp;

    i = j = 0;
    shr_pos = 64;
    shl_pos = 0;
    temp = 0;
    while ((i < rlen) && (j < alen))
    {
        word = ((temp >> shr_pos) | (a[j] << shl_pos));
        r[i] = (word & HT_BMASK);
        shr_pos -= 12, shl_pos += 12;
        if ((shr_pos > 0) && (shl_pos < 64))
            temp = a[j++];
        if (shr_pos <= 0)
            shr_pos += 64;
        if (shl_pos >= 64)
            shl_pos -= 64;
        // Any shift past 63 is undefined!
        if (shr_pos == 64)
            temp = 0;
        i++;
    }
    if (i < rlen)
        r[i++] = ((temp >> shr_pos) & HT_BMASK);
    for (; i < rlen; i++)
        r[i] = 0;
}

void mpi_conv_52to64(uint64_t *r, const uint64_t *a, int rlen, int alen)
{
    int i, j, bits_in_word, bits_to_shift;
    uint64_t word;

    i = j = 0;
    bits_in_word = bits_to_shift = 0;
    word = 0;
    while ((i < rlen) && (j < alen))
    {
        word |= (a[j] << bits_in_word);
        bits_to_shift = (64 - bits_in_word);
        bits_in_word += 52;
        if (bits_in_word >= 64)
        {
            r[i++] = word;
            word = ((bits_to_shift > 0) ? (a[j] >> bits_to_shift) : 0);
            bits_in_word = ((bits_to_shift > 0) ? (52 - bits_to_shift) : 0);
        }
        j++;
    }
    if (i < rlen)
        r[i++] = word;
    for (; i < rlen; i++)
        r[i] = 0;
}

__m512i set_vector(const uint64_t a7, const uint64_t a6, const uint64_t a5, const uint64_t a4,
                   const uint64_t a3, const uint64_t a2, const uint64_t a1, const uint64_t a0)
{
    __m512i r;

    ((uint64_t *)&r)[0] = a0;
    ((uint64_t *)&r)[1] = a1;
    ((uint64_t *)&r)[2] = a2;
    ((uint64_t *)&r)[3] = a3;
    ((uint64_t *)&r)[4] = a4;
    ((uint64_t *)&r)[5] = a5;
    ((uint64_t *)&r)[6] = a6;
    ((uint64_t *)&r)[7] = a7;

    return r;
}

void get_channel_8x1w(uint64_t *r, const htfe_t a, const int ch)
{
    int i;

    for (i = 0; i < HT_NWORDS; i++)
    {
        r[i] = ((uint64_t *)&a[i])[ch];
    }
}
