#-------------------------------------------------
#
# Project created by QtCreator 2015-12-17T17:24:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Ploter
TEMPLATE = app


SOURCES += main.cpp\
    axis.cpp \
    ploter.cpp \
    plot.cpp \
    utils.cpp \
    numeth1.cpp \
    numeth2.cpp \
    numeth3.cpp

HEADERS  += \
    axis.h \
    ploter.h \
    plot.h \
    utils.h \
    numeth1.h \
    numeth2.h \
    numeth3.h

CONFIG += c++11
