#ifndef MONTGOMERY_CURVE_H
#define MONTGOMERY_CURVE_H

#include <ak_mpzn.h>

struct wtomutils
{
    ak_uint64 A[ak_mpzn512_size];
    ak_uint64 B[ak_mpzn512_size];
    ak_uint64 malphas[ak_mpzn512_size];
    ak_uint64 Binv[ak_mpzn512_size];
    ak_uint64 A3B[ak_mpzn512_size];
    ak_uint64 pxm[ak_mpzn512_size];
    ak_uint64 pym[ak_mpzn512_size];
};

#endif // MONTGOMERY_CURVE_H
