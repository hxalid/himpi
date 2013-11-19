#ifndef HPNLA_DEBUG_H_
#define HPNLA_DEBUG_H_

#include <stdio.h>
#include <libgen.h>

/*!
 * \defgroup debug Debug macros
 * Macros to print debug information
 * \{
 */

#define DEBUG_MPI 0
#define DEBUG 0
#define CHECK_25D 1
#define CHECK_CANNON 1
#define NOT_USED(x) ( (void)(x) )

#define info_print(level, row, col, fmt, args...) \
  do { if (level <= DEBUG) fprintf(stderr, "%s:%d:%s(): processor (%zu,%zu) " fmt,\
      basename(__FILE__),                  \
      __LINE__, __func__, row, col,   \
##args);                   \
  } while (0)

#if DEBUG
#define debug_print(level, row, col, fmt, args...) \
  do {
  } while (0)
#else

#define debug_print(level, row, col, fmt, args...) \
  do { if (level <= DEBUG) fprintf(stderr, "%s:%d:%s(): processor (%zu,%zu) " fmt,\
      basename(__FILE__),                  \
      __LINE__, __func__, row, col,   \
##args);                   \
  } while (0)
#endif

#define debug_print_rank(level, rank, fmt, args...) \
  do { if (level <= DEBUG) fprintf(stderr, "%s:%d:%s(): processor (%zu) " fmt,\
      basename(__FILE__),                  \
      __LINE__, __func__, rank,   \
##args);                   \
  } while (0)

/*!
 * \}
 */

#endif /*HPNLA_DEBUG_H_*/
