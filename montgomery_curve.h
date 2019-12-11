#ifndef MONTGOMERY_CURVE_H
#define MONTGOMERY_CURVE_H

#include "ak_curves.h"

#include <sys/types.h>

#include <ak_mpzn.h>

extern int ak_wpoint_to_mpoint(ak_wpoint wp, ak_wcurve ec);

extern int ak_mpoint_to_wpoint(ak_wpoint mp, ak_wcurve ec);

extern int ak_mpoint_set_as_unit(ak_wpoint mp, ak_wcurve ec);

extern void ak_wpoint_double_montgomery(ak_wpoint wp, ak_wcurve ec);

extern void ak_mpoint_double(ak_wpoint mp, ak_wcurve ec);

extern void ak_mpoint_dadd(ak_wpoint p, ak_wpoint q, ak_wpoint pmq, ak_wcurve ec);

extern void ak_wpoint_pow_montgomery(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec);

extern void ak_wpoint_pow_montgomery_x(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec);

#endif // MONTGOMERY_CURVE_H
