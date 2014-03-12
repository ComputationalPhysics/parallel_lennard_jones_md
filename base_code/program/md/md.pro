TEMPLATE = app
CONFIG += console
CONFIG -= qt

release {

}

SOURCES += \
    system.cpp \
    statisticssampler.cpp \
    md.cpp \
    random.cpp \
    thermostat.cpp \
    cutil.cpp \
    unitconverter.cpp \
    settings.cpp \
    mdio.cpp \
    mdtimer.cpp

HEADERS += \
    system.h \
    statisticssampler.h \
    random.h \
    thermostat.h \
    cutil.h \
    cinifile.h \
    unitconverter.h \
    settings.h \
    mdio.h \
    mdtimer.h \
    atom_types.h \
    potential_lennard_jones.h

mac {
    CONFIG -= app_bundle
}

QMAKE_CXX = mpic++
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

# MPI Settings
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
