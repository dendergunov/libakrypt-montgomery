#define LIBAKRYPT_HAVE_SYSTYPES_H


#include "montgomery_curve.h"
#include "ak_curves.h"

#include <ak_mpzn.h>
#include <string.h>

int ak_wpoint_to_mpoint(ak_wpoint wp, ak_wcurve ec)
{
    if( wp == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                      "using null pointer to elliptic curve point");

    if( ec->mc == NULL) ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return ak_error_wrong_option;
    }
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

    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return ak_error_wrong_option;
    }

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

    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return ak_error_wrong_option;
    }

    ak_mpzn_set_ui(mp->x, ec->size, 1);
    ak_mpzn_set_ui(mp->z, ec->size, 0);

    return ak_error_ok;
}

void ak_wpoint_double_montgomery(ak_wpoint wp, ak_wcurve ec)
{
    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return;
    }

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

void ak_mpoint_double(ak_wpoint mp, ak_wcurve ec)
{
    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return;
    }

    if( ak_mpzn_cmp_ui(mp->z, ec->size, 0) == ak_true) return;

    ak_mpzn512 A, AA, B, BB, C;

    ak_mpzn_add_montgomery(A, mp->x, mp->z, ec->p, ec->size);
    ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

    ak_mpzn_sub(B, ec->p, mp->z, ec->size);
    ak_mpzn_add_montgomery(B, B, mp->x, ec->p, ec->size);
    ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

    ak_mpzn_sub(C, ec->p, BB, ec->size);
    ak_mpzn_add_montgomery(C, C, AA, ec->p, ec->size);

    //now it's x-coordinate of 2P in Montgomery
    ak_mpzn_mul_montgomery(mp->x, AA, BB, ec->p, ec->n, ec->size);

    ak_mpzn_mul_montgomery(mp->z, C, ec->mc->a24, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(mp->z, mp->z, BB, ec->p, ec->size);
    ak_mpzn_mul_montgomery(mp->z, mp->z, C, ec->p, ec->n, ec->size);
}

void ak_mpoint_dadd(ak_wpoint p, ak_wpoint q, ak_wpoint pmq, ak_wcurve ec)
{
    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return;
    }

    if(ak_mpzn_cmp_ui(pmq->z, ec->size, 0) == ak_true){
        ak_mpoint_to_wpoint(p, ec);
        ak_wpoint_double_montgomery(p, ec);
        ak_wpoint_to_mpoint(p, ec);
        return;
    }

    ak_mpzn512 A, B, C, D, DA, CB;

    ak_mpzn_add_montgomery(A, p->x, p->z, ec->p, ec->size); //X1 + Z1

    ak_mpzn_sub(B, ec->p, p->z, ec->size);
    ak_mpzn_add_montgomery(B, B, p->x, ec->p, ec->size); //X1 - Z1

    ak_mpzn_add_montgomery(C, q->x, q->z, ec->p, ec->size); //X2 + Z2
    ak_mpzn_sub(D, ec->p, q->z, ec->size);
    ak_mpzn_add_montgomery(D, D, q->x, ec->p, ec->size); //X2 - Z2

    ak_mpzn_mul_montgomery(DA, D, A, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(CB, C, B, ec->p, ec->n, ec->size);

    ak_mpzn512 X3, Z3;
    ak_mpzn_set(X3, pmq->x, ec->size);
    ak_mpzn_set(Z3, pmq->z, ec->size);

    ak_mpzn_add_montgomery(p->x, DA, CB, ec->p, ec->size);
    ak_mpzn_mul_montgomery(p->x, p->x, p->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(p->x, p->x, Z3, ec->p, ec->n, ec->size);

    ak_mpzn_sub(p->z, ec->p, CB, ec->size);
    ak_mpzn_add_montgomery(p->z, p->z, DA, ec->p, ec->size);
    ak_mpzn_mul_montgomery(p->z, p->z, p->z, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(p->z, p->z, X3, ec->p, ec->n, ec->size);

}

void ak_wpoint_pow_montgomery(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec)
{
    //check for errors and unnecessary calculations
    if( ec->mc == NULL){
        ak_error_message( ak_error_null_pointer, __func__,
                                          "elliptic curve does not have Montgomery alternative");
        return;
    }

    if( ak_mpzn_cmp_ui(p->z, ec->size, 0) == ak_true){
        ak_wpoint_set_as_unit(q, ec);
        return;
    }

    if( ak_mpzn_cmp_ui(k, size, 0) == ak_true){
        ak_wpoint_set_as_unit(q, ec);
        return;
    }

    //store p in montgomery form
    struct wpoint point_;
    ak_wpoint point = &point_;

    ak_wpoint_set_wpoint(point, p, ec);
    ak_wpoint_to_mpoint(point, ec);

    struct wpoint left_, right_;
    ak_wpoint left = &left_, right = &right_;

    //find first bit equals one
    long long first_ll_with_one;
    bool_t first_one_found = ak_false;

    for(first_ll_with_one = size-1; first_ll_with_one >= 0 && !first_one_found; --first_ll_with_one){
        if( k[first_ll_with_one] != 0){
            first_one_found = ak_true;
            break;
        }
    }

    ak_mpzn512 A, AA, B, BB, C, D, DA, CB, E;

    ak_uint64 uk = k[first_ll_with_one];
    long long i, j;
    first_one_found = ak_false;
    for(j = 0; j < 64 && !first_one_found; ++j){
        if( uk&0x8000000000000000LL ){
            first_one_found = ak_true;

            memcpy(left->x, point->x, ec->size*sizeof( ak_uint64 ));
            memcpy(left->z, point->z, ec->size*sizeof( ak_uint64 )); //left now contains P

            ak_mpzn_add_montgomery(A, left->x, left->z, ec->p, ec->size);
            ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

            ak_mpzn_sub(B, ec->p, left->z, ec->size);
            ak_mpzn_add_montgomery(B, B, left->x, ec->p, ec->size);
            ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

            ak_mpzn_sub(C, ec->p, BB, ec->size);
            ak_mpzn_add_montgomery(C, C, AA, ec->p, ec->size);

            //now it's x-coordinate of 2P in Montgomery
            ak_mpzn_mul_montgomery(right->x, AA, BB, ec->p, ec->n, ec->size);

            ak_mpzn_mul_montgomery(right->z, C, ec->mc->a24, ec->p, ec->n, ec->size);
            ak_mpzn_add_montgomery(right->z, right->z, BB, ec->p, ec->size);
            ak_mpzn_mul_montgomery(right->z, right->z, C, ec->p, ec->n, ec->size); //right now contains 2P
        }

        uk<<=1;
    }

    //let's compute the rest
    //cycle to compute both double and differential addition
    for( ; j < 64; ++j){ //computing in same long long

        if( uk&0x8000000000000000LL ){ //1: left -> add, right -> double

            ak_mpzn_add_montgomery(A, right->x, right->z, ec->p, ec->size);
            ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

            ak_mpzn_sub(B, ec->p, right->z, ec->size);
            ak_mpzn_add_montgomery(B, B, right->x, ec->p, ec->size);
            ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

            ak_mpzn_sub(E, ec->p, BB, ec->size);
            ak_mpzn_add_montgomery(E, AA, E, ec->p, ec->size);

            ak_mpzn_add_montgomery(C, left->x, left->z, ec->p, ec->size);

            ak_mpzn_sub(D, ec->p, left->z, ec->size);
            ak_mpzn_add_montgomery(D, left->x, D, ec->p, ec->size);

            ak_mpzn_mul_montgomery(DA, D, A, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(CB, C, B, ec->p, ec->n, ec->size);

            ak_mpzn_mul_montgomery(right->x, AA, BB, ec->p, ec->n, ec->size);

            ak_mpzn_mul_montgomery(right->z, E, ec->mc->a24, ec->p, ec->n, ec->size);
            ak_mpzn_add_montgomery(right->z, right->z, BB, ec->p, ec->size);
            ak_mpzn_mul_montgomery(right->z, right->z, E, ec->p, ec->n, ec->size);


            ak_mpzn_add_montgomery(left->x, DA, CB, ec->p, ec->size);
            ak_mpzn_mul_montgomery(left->x, left->x, left->x, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(left->x, left->x, point->z, ec->p, ec->n, ec->size);

            ak_mpzn_sub(left->z, ec->p, CB, ec->size);
            ak_mpzn_add_montgomery(left->z, left->z, DA, ec->p, ec->size);
            ak_mpzn_mul_montgomery(left->z, left->z, left->z, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(left->z, left->z, point->x, ec->p, ec->n, ec->size);


        } else { //0: left -> double, right-> add
            ak_mpzn_add_montgomery(A, left->x, left->z, ec->p, ec->size);
            ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

            ak_mpzn_sub(B, ec->p, left->z, ec->size);
            ak_mpzn_add_montgomery(B, B, left->x, ec->p, ec->size);
            ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

            ak_mpzn_sub(E, ec->p, BB, ec->size);
            ak_mpzn_add_montgomery(E, AA, E, ec->p, ec->size);

            ak_mpzn_add_montgomery(C, right->x, right->z, ec->p, ec->size);

            ak_mpzn_sub(D, ec->p, right->z, ec->size);
            ak_mpzn_add_montgomery(D, right->x, D, ec->p, ec->size);

            ak_mpzn_mul_montgomery(DA, D, A, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(CB, C, B, ec->p, ec->n, ec->size);

            ak_mpzn_mul_montgomery(left->x, AA, BB, ec->p, ec->n, ec->size);

            ak_mpzn_mul_montgomery(left->z, E, ec->mc->a24, ec->p, ec->n, ec->size);
            ak_mpzn_add_montgomery(left->z, left->z, BB, ec->p, ec->size);
            ak_mpzn_mul_montgomery(left->z, left->z, E, ec->p, ec->n, ec->size); //left now contains 2*left

            ak_mpzn_add_montgomery(right->x, DA, CB, ec->p, ec->size);
            ak_mpzn_mul_montgomery(right->x, right->x, right->x, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(right->x, point->z, right->x, ec->p, ec->n, ec->size);

            ak_mpzn_sub(right->z, ec->p, CB, ec->size);
            ak_mpzn_add_montgomery(right->z, DA, right->z, ec->p, ec->size);
            ak_mpzn_mul_montgomery(right->z, right->z, right->z, ec->p, ec->n, ec->size);
            ak_mpzn_mul_montgomery(right->z, point->x, right->z, ec->p, ec->n, ec->size); //right now conatins 2*right

        }

        uk<<=1;

    }

    --first_ll_with_one;

    for(i = first_ll_with_one; i >= 0; --i){
        uk = k[i];
        for(j = 0; j < 64; ++j){
            if( uk&0x8000000000000000LL ){ //1: left -> add, right -> double

                ak_mpzn_add_montgomery(A, right->x, right->z, ec->p, ec->size);
                ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

                ak_mpzn_sub(B, ec->p, right->z, ec->size);
                ak_mpzn_add_montgomery(B, B, right->x, ec->p, ec->size);
                ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

                ak_mpzn_sub(E, ec->p, BB, ec->size);
                ak_mpzn_add_montgomery(E, AA, E, ec->p, ec->size);

                ak_mpzn_add_montgomery(C, left->x, left->z, ec->p, ec->size);

                ak_mpzn_sub(D, ec->p, left->z, ec->size);
                ak_mpzn_add_montgomery(D, left->x, D, ec->p, ec->size);

                ak_mpzn_mul_montgomery(DA, D, A, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(CB, C, B, ec->p, ec->n, ec->size);

                ak_mpzn_mul_montgomery(right->x, AA, BB, ec->p, ec->n, ec->size);

                ak_mpzn_mul_montgomery(right->z, E, ec->mc->a24, ec->p, ec->n, ec->size);
                ak_mpzn_add_montgomery(right->z, right->z, BB, ec->p, ec->size);
                ak_mpzn_mul_montgomery(right->z, right->z, E, ec->p, ec->n, ec->size);


                ak_mpzn_add_montgomery(left->x, DA, CB, ec->p, ec->size);
                ak_mpzn_mul_montgomery(left->x, left->x, left->x, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(left->x, left->x, point->z, ec->p, ec->n, ec->size);

                ak_mpzn_sub(left->z, ec->p, CB, ec->size);
                ak_mpzn_add_montgomery(left->z, left->z, DA, ec->p, ec->size);
                ak_mpzn_mul_montgomery(left->z, left->z, left->z, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(left->z, left->z, point->x, ec->p, ec->n, ec->size);


            } else { //0: left -> double, right-> add
                ak_mpzn_add_montgomery(A, left->x, left->z, ec->p, ec->size);
                ak_mpzn_mul_montgomery(AA, A, A, ec->p, ec->n, ec->size);

                ak_mpzn_sub(B, ec->p, left->z, ec->size);
                ak_mpzn_add_montgomery(B, B, left->x, ec->p, ec->size);
                ak_mpzn_mul_montgomery(BB, B, B, ec->p, ec->n, ec->size);

                ak_mpzn_sub(E, ec->p, BB, ec->size);
                ak_mpzn_add_montgomery(E, AA, E, ec->p, ec->size);

                ak_mpzn_add_montgomery(C, right->x, right->z, ec->p, ec->size);

                ak_mpzn_sub(D, ec->p, right->z, ec->size);
                ak_mpzn_add_montgomery(D, right->x, D, ec->p, ec->size);

                ak_mpzn_mul_montgomery(DA, D, A, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(CB, C, B, ec->p, ec->n, ec->size);

                ak_mpzn_mul_montgomery(left->x, AA, BB, ec->p, ec->n, ec->size);

                ak_mpzn_mul_montgomery(left->z, E, ec->mc->a24, ec->p, ec->n, ec->size);
                ak_mpzn_add_montgomery(left->z, left->z, BB, ec->p, ec->size);
                ak_mpzn_mul_montgomery(left->z, left->z, E, ec->p, ec->n, ec->size); //left now contains 2*left

                ak_mpzn_add_montgomery(right->x, DA, CB, ec->p, ec->size);
                ak_mpzn_mul_montgomery(right->x, right->x, right->x, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(right->x, point->z, right->x, ec->p, ec->n, ec->size);

                ak_mpzn_sub(right->z, ec->p, CB, ec->size);
                ak_mpzn_add_montgomery(right->z, DA, right->z, ec->p, ec->size);
                ak_mpzn_mul_montgomery(right->z, right->z, right->z, ec->p, ec->n, ec->size);
                ak_mpzn_mul_montgomery(right->z, point->x, right->z, ec->p, ec->n, ec->size); //right now conatins 2*right

            }

            uk<<=1;
        }
    }

    memcpy(q->x, left->x, ec->size*sizeof( ak_uint64));
    memcpy(q->z, left->z, ec->size*sizeof( ak_uint64));

    //switch back to weierstrass
    ak_mpoint_to_wpoint(q, ec);

}


