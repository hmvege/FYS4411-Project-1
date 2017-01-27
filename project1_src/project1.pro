TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    function.cpp \
    state.cpp \
    singlestate.cpp

HEADERS += \
    state.h \
    singlestate.h
