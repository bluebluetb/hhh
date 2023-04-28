#include "random.h"

unsigned int M_MT[624];
unsigned int M_ind = 0;

unsigned int MT[624];
unsigned int ind = 0;

void random_init(unsigned int seed) {
    int i;
    MT[0] = seed;
    for (i = 1; i < 624; i++) {
        MT[i] = (unsigned int)(1812433253u * (((unsigned long)MT[i-1]) ^ (((unsigned long)MT[i-1]) << 30)) + i);
    }
    ind=1;
}

void generateNext() {
    int i;
    for (i = 0; i < 624; i++) {
        unsigned int y = (MT[i] & (1 << 31)) + (MT[(i+1) % 624] & ~(1 << 31));
        MT[i] = MT[(i + 397) % 624] ^ (y >> 1);
        if (y % 2 != 0) {
            MT[i] ^= (2567483615u);
        }
    }
}

unsigned int random_int() {
    if (ind == 0) {
        generateNext();
    }

    unsigned int y = MT[ind];
    y = y ^ (y >> 11);
    y = y ^ ((y << 7) & (2636928640u));
    y = y ^ ((y << 15) & (4022730752u));
    y = y ^ (y >> 18);

    ind = (ind + 1) % 624;
    return y;
}

unsigned int random_count_bits(unsigned int low, unsigned int high) {
    unsigned int result = 0;
    while (1u << result <= high - low) {
        result++;
    }
    return result;
}

unsigned int random_range(unsigned int low, unsigned int high, int bits) {
    unsigned int guess;
    do {
        guess = random_int() % (1 << bits);
    } while (guess > high - low);
    return guess + low;
}

const double maxInt = (double)0xffffffffu;
double random_range_d(double low, double high) {
    return random_int() * ((high - low) / maxInt) + low;
}

void random_mark_state() {
    int i;
    M_ind = ind;
    for (i = 0; i < 624; i++) M_MT[i] = MT[i];
}

void random_restore_state() {
    int i;
    ind = M_ind;
    for (i = 0; i < 624; i++) MT[i] = M_MT[i];
}
