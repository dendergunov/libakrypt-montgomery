#define LIBAKRYPT_HAVE_SYSTYPES_H


#include "montgomery_curve.h"
#include "ak_curves.h"

#include <ak_mpzn.h>

int ak_wpoint_to_mpoint(ak_wpoint wp, ak_wcurve ec)
{
    if( wp == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                      "using null pointer to elliptic curve point");

    if( ec->mc == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");

    //check for infinity point
    if( ak_mpzn_cmp_ui(wp->z, ec->size, 0) == ak_true){
        ak_mpzn_set_ui(wp->x, ec->size, 1);
        return ak_error_ok;
    }

    //x-coordinate
    ak_mpzn_mul_montgomery(wp->x, wp->x, ec->mc->b, ec->p, ec->n, ec->size);

    ak_mpzn512 m_alpha_b;
    ak_mpzn_mul_montgomery(m_alpha_b, wp->z, ec->mc->m_alpha_b_r, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(wp->x, wp->x, m_alpha_b, ec->p, ec->size);

    //y-coordinate
    ak_mpzn_mul_montgomery(wp->y, wp->y, ec->mc->b, ec->p, ec->n, ec->size);

    return ak_error_ok;
}

int ak_mpoint_to_wpoint(ak_wpoint mp, ak_wcurve ec)
{
    if( mp == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                      "using null pointer to elliptic curve point");

    if( ec->mc == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");

    //check for infinity point
    if (ak_mpzn_cmp_ui(mp->z, ec->size, 0) == ak_true){
        ak_mpzn_set_ui(mp->x, ec->size, 0);
        return ak_error_ok;
    }

    //x-coordinate
    ak_mpzn_mul_montgomery(mp->x, mp->x, ec->mc->b_inv, ec->p, ec->n, ec->size);

    ak_mpzn512 alpha;
    ak_mpzn_mul_montgomery(alpha, mp->z, ec->mc->alpha_r, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(mp->x, mp->x, alpha, ec->p, ec->size);

    //y-coordinate
    ak_mpzn_mul_montgomery(mp->y, mp->y, ec->mc->b_inv, ec->p, ec->n, ec->size);

    return ak_error_ok;
}

int ak_mpoint_set_as_unit(ak_wpoint mp, ak_wcurve ec)
{
    if( mp == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                      "using null pointer to elliptic curve point");

    if( ec->mc == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");

    ak_mpzn_set_ui(mp->x, ec->size, 1);
    ak_mpzn_set_ui(mp->z, ec->size, 0);

    return ak_error_ok;
}

void ak_wpoint_double_montgomery(ak_wpoint wp, ak_wcurve ec)
{
    if( wp == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                      "using null pointer to elliptic curve point");

    if( ec->mc == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");

    if( ak_mpzn_cmp_ui(wp->z, ec->size, 0) == ak_true) return;

    //Turn x-coordinate to Montgomery curve
    ak_mpzn_mul_montgomery(wp->x, wp->x, ec->mc->b, ec->p, ec->n, ec->size);

    ak_mpzn512 m_alpha_b;
    ak_mpzn_mul_montgomery(m_alpha_b, ec->mc->m_alpha_b_r, wp->z, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(wp->x, wp->x, m_alpha_b, ec->p, ec->size);

    //Now it's Montgomery point

    ak_mpzn512 A, AA, B, BB, C;

    ak_mpzn_add_montgomery(A, wp->x, wp->z, ec->p, ec->size);
    ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

    ak_mpzn_sub(B, ec->p, wp->z, ec->size);
    ak_mpzn_add_montgomery(B, B, wp->x, ec->p, ec->size);
    ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

    ak_mpzn_sub(C, ec->p, BB, ec->size);
    ak_mpzn_add_montgomery(C, C, AA, ec->p, ec->size);

    //now it's x-coordinate of 2P in Montgomery
    ak_mpzn_mul_montgomery(wp->x, AA, BB, ec->p, ec->n, ec->size);

    ak_mpzn_mul_montgomery(wp->z, C, ec->mc->a24, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(wp->z, wp->z, BB, ec->p, ec->size);
    ak_mpzn_mul_montgomery(wp->z, wp->z, C, ec->p, ec->n, ec->size);

    //turn back into weierstrass form
    ak_mpzn_mul_montgomery(wp->x, wp->x, ec->mc->b_inv, ec->p, ec->n, ec->size);

    //take allocated m_alpha_b to compute alpha*wp->z
    ak_mpzn_mul_montgomery(m_alpha_b, wp->z, ec->mc->alpha_r, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(wp->x, wp->x, m_alpha_b, ec->p, ec->size);
    //now wpoint contains (2P.x, P.y, 2P.z) (it's suitable only to check x-coordinates
}

/*
void ak_wpoint_double_montgomery(ak_wpoint p, ak_wcurve ec, ak_mcurve mc)
{
    //Now I accept point in Montgomery, that's not right
    //TODO: add check for infinity point

    //ak_mpzn_mul_montgomery(p->x, p->x, mc->B, ec->p, ec->n, ec->size);
    //ak_mpzn_add_montgomery(p->x, p->x, mc->malphaB, ec->p, ec->size);

    //точка в монтгомери
    ak_mpznmax A, AA, B, BB, C;

    ak_mpzn_add_montgomery(A, p->x, p->z, ec->p, ec->size);
    ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);//A^2*R^{-1}

    ak_mpzn_set(B, ec->p, ec->size);
    ak_mpzn_sub(B, B, p->z, ec->size);
    ak_mpzn_add_montgomery(B, B, p->x, ec->p, ec->size);
    ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);//B^2*R^{-1}

    ak_mpzn_set(C, ec->p, ec->size);
    ak_mpzn_sub(C, C, BB, ec->size);
    ak_mpzn_add_montgomery(C, C, AA, ec->p, ec->size);//(AA-BB)*R^{-1}

    ak_mpzn_mul_montgomery(p->x, AA, BB, ec->p, ec->n, ec->size);//AA*BB*R^{-3}

    ak_mpzn_mul_montgomery(p->z, C, mc->A24, ec->p, ec->n, ec->size);//C*a24*R^{-1}
    ak_mpzn_add_montgomery(p->z, p->z, BB, ec->p, ec->size);//(BB+C*a24)*R^{-1}
    ak_mpzn_mul_montgomery(p->z, p->z, C, ec->p, ec->n, ec->size);//C*(BB+C*a24)*R^{-3}

    //ak_mpzn_set(ma, ec->p, ec->size);
    //ak_mpzn_sub(ma, ma, mc->alphamR3, ec->size);

    //ak_mpzn_mul_montgomery(p->x, p->x, mc->Binv, ec->p, ec->n, ec->size);
    //ak_mpzn_add_montgomery(p->x, p->x, ma, ec->p, ec->size);
}

void ak_mpoint_dadd(ak_wpoint p, ak_wpoint q, ak_wpoint pmq, ak_wcurve ec, ak_mcurve mc);
void ak_wpoint_pow_montgomery(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec);
*/
