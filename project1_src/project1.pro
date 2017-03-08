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


# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK


# Following to make openmp usable on linux
#QMAKE_LFLAGS += -fopenmp

# Following to make openmp usable on mac
#QMAKE_LDFLAGS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib

# Following used to make armadillo usable on mac
LIBS += -L/usr/local/lib -larmadillo
INCLUDEPATH += /usr/local/include

#INCLUDEPATH += -I/usr/local/include
#INCLUDEPATH += -L/usr/local/lib
#compileCommand-I/usr/local/include -L/usr/local/lib -llapack -lblas -larmadillo
