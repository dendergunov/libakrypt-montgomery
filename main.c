#define LIBAKRYPT_HAVE_SYSTYPES_H

#include <sys/types.h>

#include <libakrypt.h>
#include <ak_curves.h>
#include <ak_parameters.h>
#include <ak_mpzn.h>

#include "montgomery_curve.h"

int main()
{

    ak_mcurve mc = &wtomutils256;
    ak_wcurve wc = &id_tc26_gost_3410_2012_256_paramSetA;

    struct wpoint point_, point2_;
    ak_wpoint point = &point_, point2 = &point2_;

    ak_wpoint_set(point, wc); //right
    ak_wpoint_set(point2, wc);


    ak_wpoint_wmw(point2, wc, mc);

    ak_wpoint_double_montgomery(point2, wc, mc);
    ak_wpoint_double(point, wc);

    ak_wpoint_reduce(point2, wc);
    ak_wpoint_reduce(point, wc);


    //ak_uint64 A24[ak_mpzn256_size] = {0x53b3e54bcadb1ba2, 0x3cc58dc21b53907f, 0x0fc95ab8b9c2278c, 0x2af35aad0b3de2f6};
    //ak_

    return 0;
}
