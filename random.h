#ifndef RANDOM_H2
#define RANDOM_H2

#define STATIC_ASSERT(COND,MSG) typedef char static_assertion_##MSG[(COND)?1:-1]
STATIC_ASSERT(sizeof(unsigned int)==4,int_is_4_bytes);

extern unsigned int M_MT[624];
extern unsigned int M_ind;

extern unsigned int MT[624];
extern unsigned int ind;

/**
 * Initializes the random number generator.
 * Implementations of rand() vary between operating systems. 
 * Using this random generator will guarantee that the results can be reproduced on other machines.
 */
void random_init(unsigned int seed);

/**
 * Will produce a random unsigned integer.
 * This integer is in the range 0..2^32-1.
 */
unsigned int random_int();

/**
 * Generates a random integer value uniformly selected between low and high.
 * The bits argument is computed using random_count_bits(low,high).
 */
unsigned int random_range(unsigned int low, unsigned int high, int bits);

/**
 * Computes the bit argument for the random_range function.
 */
unsigned int random_count_bits(unsigned int low, unsigned int high);

/**
 * Stores the current state of the random number generator. 
 * Can later be restored using random_restore_state().
 * Not that this replaces the previously marked state.
 */
void random_mark_state();

/**
 * Loads the state of the random generator marked by random_mark_state().
 * After a call to this method, the random functions will produce the 
 * same sequence of random numbers as after the last call to random_mark_state().
 */
void random_restore_state();

/**
 * Computes a random floating point value uniformly selected between low and high. 
 */
double random_range_d(double low, double high);



#endif
