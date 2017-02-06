TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    functions.cpp \
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
    gaussianhermitequadrature.h \
    functions.h

LIBS += -llapack -lblas -larmadillo

#INCLUDEPATH += -I/usr/local/include
#INCLUDEPATH += -L/usr/local/lib
#compileCommand-I/usr/local/include -L/usr/local/lib -llapack -lblas -larmadillo
