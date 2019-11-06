#define LIBAKRYPT_HAVE_SYSTYPES_H

#include "montgomery_curve.h"
#include "ak_curves.h"
#include "ak_parameters.h"

#include <libakrypt.h>
#include <ak_mpzn.h>

#include <sys/types.h>
#include <stdio.h>

#include <time.h>


int main()
{

    ak_wcurve wc = &id_tc26_gost_3410_2012_256_paramSetA;

    struct wpoint point_, point2_, minf_;
    ak_wpoint point = &point_, point2 = &point2_, minf = &minf_;

    ak_wpoint_set(point, wc); //right

    ak_mpzn256 q;
    ak_mpzn_set(q, wc->q, wc->size);
    ak_mpzn256 k = ak_mpzn256_one;

    ak_mpzn_add(q, q, k, wc->size);

    clock_t start, end;
    double weir_time, mont_time;

    start = clock();
    ak_wpoint_pow(point2, point, q, wc->size, wc);
    end = clock();
    weir_time = ((double) end - start) / CLOCKS_PER_SEC;


    ak_wpoint_reduce(point2, wc);

    start = clock();
    ak_wpoint_pow_montgomery(minf, point, q, wc->size, wc);
    end = clock();
    mont_time = ((double) end - start) / CLOCKS_PER_SEC;

    ak_wpoint_reduce(minf, wc);

    printf("%f %f\n", weir_time, mont_time);
    printf("x-coord: %d\n", ak_mpzn_cmp(point2->x, minf->x, wc->size));

    for(int i = 1; i < 100; ++i){
        ak_mpzn_set_ui(k, wc->size, i);

        ak_wpoint_pow_montgomery(point2, point, k, wc->size, wc);

        ak_wpoint_pow(minf, point, k, wc->size, wc);

        ak_wpoint_reduce(point2, wc);
        ak_wpoint_reduce(minf, wc);

        printf("%d: x-coord: %d\n", i, ak_mpzn_cmp(point2->x, minf->x, wc->size));
    }

        return 0;
}
