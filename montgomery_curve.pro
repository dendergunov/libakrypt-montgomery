TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ak_sign.c \
        main.c \
        montgomery_curve.c

unix:!macx: LIBS += -L$$PWD/../../../../usr/local/lib/ -lakrypt-static

INCLUDEPATH += $$PWD/../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../usr/local/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../usr/local/lib/libakrypt-static.a

HEADERS += \
    ak_curves.h \
    ak_parameters.h \
    ak_sign.h \
    montgomery_curve.h
