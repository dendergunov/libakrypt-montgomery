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


    //check that functions work right (it doesn't)
    ak_wpoint_w_to_m(point, wc, mc);
    ak_wpoint_w_to_m(point2, wc, mc);

    ak_wpoint_m_to_w(point, wc, mc);
    ak_wpoint_m_to_w(point2, wc, mc);
    //-------------------------------


    //
    ak_wpoint_w_to_m(point2, wc, mc);
    ak_wpoint_double_montgomery(point2, wc, mc);
    ak_wpoint_reduce(point2, wc);

    ak_wpoint_double(point, wc);
    //ak_wpoint_w_to_m(point, wc, mc);
    ak_wpoint_reduce(point, wc);

    //ak_wpoint_m_to_w(point, wc, mc);
    ak_wpoint_m_to_w(point2, wc, mc);
    //ok in montgomery they are right


    return 0;
}
