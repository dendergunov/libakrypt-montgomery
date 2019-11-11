#define LIBAKRYPT_HAVE_SYSTYPES_H

#include "montgomery_curve.h"
#include "ak_curves.h"
#include "ak_parameters.h"

#include <libakrypt.h>
#include <ak_mpzn.h>

#include <sys/types.h>
#include <stdio.h>

#include <time.h>

int ak_wpoint_cmp_x(ak_wpoint p1, ak_wpoint p2, ak_wcurve wc);

int main()
{

    ak_wcurve wc;

    struct wpoint point_, point2_, minf_;
    ak_wpoint point = &point_, point2 = &point2_, minf = &minf_;


    clock_t start, end;
    double weir_time, mont_time;

    for(int i =0; i <2; ++i){
        if(i){
            wc = &id_tc26_gost_3410_2012_512_paramSetC;
            printf("\n\nid_tc26_gost_3410_2012_512_paramSetC curve\n");
        }else{
            wc = &id_tc26_gost_3410_2012_256_paramSetA;
            printf("id_tc26_gost_3410_2012_256_paramSetA curve\n");
        }

        ///////////////////////////////////////////////////////////////
        printf("\nTest 0: Weir(Mont(point)) == point\n");
        ak_wpoint_set(point, wc);

        ak_wpoint_set(point2, wc);
        ak_wpoint_to_mpoint(point2, wc);
        ak_mpoint_to_wpoint(point2, wc);

        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("Ok!\n"):printf("Fail!\n");

        ///////////////////////////////////////////////////////////////
        printf("\nTest 1: double point (weierstrass & montgomery)\n");
        ak_wpoint_set(point, wc);
        ak_wpoint_set(point2, wc);

        start = clock();
        ak_wpoint_double(point, wc);
        end = clock();
        weir_time = ((double) end - start) / CLOCKS_PER_SEC;

        start = clock();
        ak_wpoint_double_montgomery(point2, wc);
        end = clock();
        mont_time = ((double) end - start) / CLOCKS_PER_SEC;

        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("x(2P)_w = x(2P)_m\n"):printf("x(2P)_w != x(2P)_m\n");
        printf("Weierstrass time: %f secs, Montgomery time %f secs\n", weir_time, mont_time);

        ///////////////////////////////////////////////////////////////
        printf("\nTest 2: add point (weierstrass & montgomery)\n");
        ak_wpoint_set(point, wc);
        ak_wpoint_set(point2, wc);
        ak_wpoint_set_as_unit(minf, wc);
        ak_wpoint_to_mpoint(minf, wc);

        start = clock();
        ak_wpoint_add(point, point, wc);
        end = clock();
        weir_time = ((double) end - start) / CLOCKS_PER_SEC;

        start = clock();
        ak_wpoint_to_mpoint(point2, wc);
        ak_mpoint_dadd(point2, point2, minf, wc);
        ak_mpoint_to_wpoint(point2, wc);
        end = clock();
        mont_time = ((double) end - start) / CLOCKS_PER_SEC;

        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("x(P+P)_w = x(P+P)_m\n"):printf("x(P+P)_w != x(P+P)_m\n");
        printf("Weierstrass time: %f secs, Montgomery time %f secs\n", weir_time, mont_time);

        ///////////////////////////////////////////////////////////////
        printf("\nTest 3: point^q (weierstrass & montgomery)\n");
        ak_wpoint_set(point, wc);
        ak_wpoint_set(point2, wc);

        start = clock();
        ak_wpoint_pow(point, point, wc->q, wc->size, wc);
        end = clock();
        weir_time = ((double) end - start) / CLOCKS_PER_SEC;

        start = clock();
        ak_wpoint_pow_montgomery(point2, point2, wc->q, wc->size, wc);
        end = clock();
        mont_time = ((double) end - start) / CLOCKS_PER_SEC;

        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("x(q*P)_w = x(q*P)_m\n"):printf("x(q*P)_w != x(q*P)_m\n");
        printf("Weierstrass time: %f secs, Montgomery time %f secs\n", weir_time, mont_time);

        ///////////////////////////////////////////////////////////////
        printf("\nTest 4: point^(q) montgomery & unit weierstrass)\n");
        ak_wpoint_set_as_unit(point, wc);
        ak_wpoint_set(point2, wc);

        ak_wpoint_pow_montgomery(point2, point2, wc->q, wc->size, wc);
        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("x(unit)_w = x((q)*P)_m\n"):printf("x(unit)_w != x((q)*P)_m\n");

        ///////////////////////////////////////////////////////////////
        printf("\nTest 5: point^(q+1) montgomery & point weierstrass)\n");
        ak_wpoint_set(point, wc);
        ak_wpoint_set(point2, wc);

        ak_mpzn512 k, one = ak_mpzn512_one;
        ak_mpzn_set(k, wc->q, wc->size);
        ak_mpzn_add(k, k, one, wc->size);

        ak_wpoint_pow_montgomery(point2, point2, k, wc->size, wc);
        (ak_wpoint_cmp_x(point, point2, wc) == ak_true)?printf("x(P)_w = x((q+1)*P)_m\n"):printf("x(P)_w != x((q+1)*P)_m\n");
    }

    return 0;
}

int ak_wpoint_cmp_x(ak_wpoint p1, ak_wpoint p2, ak_wcurve wc)
{
    ak_mpzn512 a,b;

    bool_t p1isinf, p2isinf;
    p1isinf = ak_mpzn_cmp_ui(p1->z, wc->size, 0);
    p2isinf = ak_mpzn_cmp_ui(p2->z, wc->size, 0);

    if(p1isinf^p2isinf)
        return ak_false;

    ak_mpzn_mul_montgomery(a, p1->x, p2->z, wc->p, wc->n, wc->size);
    ak_mpzn_mul_montgomery(b, p2->x, p1->z, wc->p, wc->n, wc->size);

    int cmp = ak_mpzn_cmp(a, b, wc->size);

    if(cmp == 0)
        return ak_true;
    else
        return ak_false;
}
