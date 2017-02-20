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

# Following to make openmp usable on linux
QMAKE_LFLAGS += -fopenmp

# Following to make openmp usable on mac
#QMAKE_CFLAGS_RELEASE += -fopenmp
#QMAKE_LDFLAGS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib

# Following used to make armadillo usable on mac
#LIBS += -L/usr/local/lib -larmadillo
#INCLUDEPATH += /usr/local/include

#INCLUDEPATH += -I/usr/local/include
#INCLUDEPATH += -L/usr/local/lib
#compileCommand-I/usr/local/include -L/usr/local/lib -llapack -lblas -larmadillo
