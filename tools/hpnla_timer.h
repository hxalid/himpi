#ifndef timer_H_
#define timer_H_

#include "config.h"
#if defined HAVE_LIBRT
# include <time.h>
#else
# include <sys/time.h>
#endif

/*!
 * \defgroup timer Timing functions
 * Generic interface to time some parts of the code
 * \{
 */

inline double get_microsecond(struct timespec *res);
inline double get_nanosecond(struct timespec *res);
void init_timer();

int get_res(struct timespec *res);

int get_time(struct timespec *tp);
double get_timediff(struct timespec *start, struct timespec * end);

/*!
 * \}
 */

#endif /*timer_H_*/
