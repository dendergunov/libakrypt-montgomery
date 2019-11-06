#define LIBAKRYPT_HAVE_SYSTYPES_H

#include "ak_curves.h"
#include "ak_parameters.h"

#include <sys/types.h>

#include <libakrypt.h>


#include <ak_mpzn.h>

#include "montgomery_curve.h"

int main()
{

    ak_wcurve wc = &id_tc26_gost_3410_2012_256_paramSetA;

    struct wpoint point_, point2_, minf_;
    ak_wpoint point = &point_, point2 = &point2_, minf = &minf_;

    ak_wpoint_set(point, wc); //right
    ak_wpoint_set(point2, wc);
    ak_wpoint_set(minf, wc);

    ak_wpoint_to_mpoint(point, wc);
    ak_wpoint_to_mpoint(point2, wc);
    ak_wpoint_to_mpoint(minf, wc);

    ak_mpoint_double(minf, wc);

    ak_mpoint_dadd(point, point, minf, wc);

    ak_mpoint_to_wpoint(point, wc);
    ak_wpoint_reduce(point, wc);


    return 0;
}
