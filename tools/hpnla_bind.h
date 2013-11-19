#ifndef HPNLA_BIND_H_
#define HPNLA_BIND_H_

#ifdef __cplusplus
extern "C" {
#endif
/*!
 * we provide a function bind which bind a process on a processor
 * this module could be extend to the memory too.
 */
int hpnla_bind_process(char* str_mask);

#ifdef __cplusplus
}
#endif

#endif /*HPNLA_BIND_H_*/
