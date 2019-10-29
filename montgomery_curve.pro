TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.c

unix:!macx: LIBS += -L$$PWD/../../../../usr/local/lib/ -lakrypt-static

INCLUDEPATH += $$PWD/../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../usr/local/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../usr/local/lib/libakrypt-static.a
