#ifndef MONTGOMERY_CURVE_H
#define MONTGOMERY_CURVE_H

#include "ak_curves.h"

#include <sys/types.h>

#include <ak_mpzn.h>

/*! \brief Преобразование точки в форме Вейерштрасса в точку в форме Монтгомери */
extern int ak_wpoint_to_mpoint(ak_wpoint wp, ak_wcurve ec);
/*! \brief Преобразование точки в форме Монтгомери в точку в форме Вейерштрасса */
extern int ak_mpoint_to_wpoint(ak_wpoint mp, ak_wcurve ec);
/*! \brief Инициализация и присвоение контексту значения бесконечно удаленной точки эллиптической кривой в форме Монтгомери. */
extern int ak_mpoint_set_as_unit(ak_wpoint mp, ak_wcurve ec);
/*! \brief Удвоение с помощью кривой в форме Монтгомери точки эллиптической кривой, заданной в короткой форме Вейерштрасса. */
extern void ak_wpoint_double_montgomery(ak_wpoint wp, ak_wcurve ec);
/*! \brief Удвоение точки эллиптической кривой, заданной в форме Монтгомери */
extern void ak_mpoint_double(ak_wpoint mp, ak_wcurve ec);
/*! \brief Прибавление к одной точке другой в форме Монтгомери по их разности */
extern void ak_mpoint_dadd(ak_wpoint p, ak_wpoint q, ak_wpoint pmq, ak_wcurve ec);
/*! \brief Вычисление кратной точки эллиптической кривой в короткой форме Вейерштрасса, с помощью кривой Монтгомери. */
extern void ak_wpoint_pow_montgomery(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec);
/*! \brief Вычисление х-координаты кратной точки эллиптической кривой в короткой форме Вейерштрасса, с помощью кривой Монтгомери. */
extern void ak_wpoint_pow_montgomery_x(ak_wpoint q, ak_wpoint p, ak_uint64* k, size_t size, ak_wcurve ec);

#endif // MONTGOMERY_CURVE_H
