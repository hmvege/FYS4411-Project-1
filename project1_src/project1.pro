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
    gaussianhermitequadrature.cpp \
    hartreefock.cpp \
    quantumdot.cpp \
    Coulomb_Functions.cpp

HEADERS += \
    singlestate.h \
    basis.h \
    hermitepolynomials.h \
    unittests.h \
    gaussianhermitequadrature.h \
    functions.h \
    hartreefock.h \
    quantumdot.h \
    Coulomb_Functions.h

LIBS += -llapack -lblas -larmadillo


# Following used to make code usable on mac
LIBS += -L/usr/local/lib -larmadillo
INCLUDEPATH += /usr/local/include

#INCLUDEPATH += -I/usr/local/include
#INCLUDEPATH += -L/usr/local/lib
#compileCommand-I/usr/local/include -L/usr/local/lib -llapack -lblas -larmadillo
