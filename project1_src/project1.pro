TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    function.cpp \
    singlestate.cpp \
    basis.cpp \
    hermitepolynomials.cpp \
    unittests.cpp \
    gaussianhermitequadrature.cpp

HEADERS += \
    singlestate.h \
    basis.h \
    hermitepolynomials.h \
    unittests.h \
    gaussianhermitequadrature.h
