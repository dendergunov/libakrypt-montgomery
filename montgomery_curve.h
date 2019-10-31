#ifndef MONTGOMERY_CURVE_H
#define MONTGOMERY_CURVE_H

#include <ak_mpzn.h>
#include <ak_curves.h>

struct wtomutils;
typedef struct wtomutils* ak_mcurve;

struct wtomutils
{
    ak_uint64 A[ak_mpzn512_size];
    ak_uint64 B[ak_mpzn512_size];
    ak_uint64 Binv[ak_mpzn512_size];
    ak_uint64 malphaB[ak_mpzn512_size]; //w -> m x -> x*B + malphas
    ak_uint64 alpha[ak_mpzn512_size];
    ak_uint64 A24[ak_mpzn512_size];
};

/*--------------------------------------------*/
//a = "C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335"
//b = "295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513"
//p = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97"
//q = "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67"
//px = "91E38443A5E82C0D880923425712B2BB658B9196932E02C78B2582FE742DAA28"
//py = "32879423AB1A0375895786C4BB46E9565FDE0B5344766740AF268ADB32322E5C"
//A = "abcd6ab42cf78bd83f256ae2e7089e30f31637086d4e41fd4ecf952f2b6c6e86"
//B = "abcd6ab42cf78bd83f256ae2e7089e30f31637086d4e41fd4ecf952f2b6c6e88"
//alpha = "100fe73f595ff158e974b44d478d9588744fe5c192ac47ea63075dce7a14aaa"
//alphas = "394478e6b9a7d948150c78f64d02df65a65cbd02cf1a15ff1a4531ba63cecf82"
//malphas = "c6bb8719465826b7eaf38709b2fd209a59a342fd30e5ea00e5bace459c312e15"
//Binv = "7e7e82520f9f015faa1d0f18c14ab9fb35188275da3fd94206b74f34a48e0ecd"
//pxm = "cbb8f5ebd80486b923ebfb17e5464173144cac7b0447717b0ea8de20545a6a23"
//pym = "370e3a4d38005921ef122701d68f401c8b685c09767ba6448af94c29df1aa555"
//A24 = "2af35aad0b3de2f60fc95ab8b9c2278c3cc58dc21b53907f53b3e54bcadb1ba2"
/*-------------------------------------------*/
static struct wtomutils wtomutils256 = {
{0xf24e8eafa85a46c4, 0xe08aa34f719d0b82, 0x312e9cded3c543f9, 0x12162c4060a00c30},//A (montgomery representation)
{0xf24e8eafa85a4b96, 0xe08aa34f719d0b82, 0x312e9cded3c543f9, 0x12162c4060a00c30},//B (montgomery representation)
{0x2fcde5e09a6488c5, 0xf8126e0b03e2a022, 0x962a9dd1a3e72,    0xdee817c7a63a4f91},//Binv (montgomery representation)
{0xe5bace459c312e15, 0x59a342fd30e5ea00, 0xeaf38709b2fd209a, 0xc6bb8719465826b7},//malphaB
{0xa63075dce7a14aaa, 0x8744fe5c192ac47e, 0x8e974b44d478d958, 0x0100fe73f595ff15},//alpha
{0xbc93a3abea1691b1, 0x7822a8d3dc6742e0, 0xc4ba737b4f150fe,  0x84858b101828030c}//A24(montgomery representation)
};

void ak_wpoint_wmw(ak_wpoint p, ak_wcurve ec, ak_mcurve mc);
void ak_wpoint_double_montgomery(ak_wpoint p, ak_wcurve ec, ak_mcurve mc);
#endif // MONTGOMERY_CURVE_H
