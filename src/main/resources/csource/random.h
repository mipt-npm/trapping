/* $Id$ */
/*
 * random.h
 *
 * written by Sebastian Vц╤cking <seb.voeck@uni-muenster.de>
 *
 * The random module is a frontend for different random generators.
 */

/*
 * The different generators
 */
typedef enum
{
  RANDOM_STDLIB, /* standard c random generator (the default one)*/
  RANDOM_CW,
  RANDOM_JAMES
} RandomMethode;

/*
 * Selects a generator
 */
void random_set_method(RandomMethode method);

/*
 * Returns a number between 0 and 1 generated with the selected generator.
 */
double random_get();

/*
 * Seeds the random generator if it supports it.
 */
void random_seed(int seed);
