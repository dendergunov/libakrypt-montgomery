#define LIBAKRYPT_HAVE_SYSTYPES_H


#include "montgomery_curve.h"

#include <ak_mpzn.h>
#include <ak_curves.h>


void ak_wpoint_w_to_m(ak_wpoint p, ak_wcurve ec, ak_mcurve mc)
{
    ak_mpzn_mul_montgomery(p->x, p->x, mc->B, ec->p, ec->n, ec->size);

    ak_mpznmax malphaB;
    ak_mpzn_mul_montgomery(malphaB, mc->malphaBR, p->z, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(p->x, p->x, malphaB, ec->p, ec->size);

    ak_mpzn_mul_montgomery(p->y, p->y, mc->B, ec->p, ec->n, ec->size);
}

void ak_wpoint_m_to_w(ak_wpoint p, ak_wcurve ec, ak_mcurve mc)
{
    ak_mpzn_mul_montgomery(p->x, p->x, mc->Binv, ec->p, ec->n, ec->size);

    ak_mpznmax alpha;
    ak_mpzn_mul_montgomery(alpha, p->z, mc->alphaR, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(p->x, p->x, mc->alpha, ec->p, ec->size);

    ak_mpzn_mul_montgomery(p->y, p->y, mc->Binv, ec->p, ec->n, ec->size);
}



void ak_wpoint_double_montgomery(ak_wpoint p, ak_wcurve ec, ak_mcurve mc)
{
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
