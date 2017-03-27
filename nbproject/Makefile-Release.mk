#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/ChowLiu.o \
	${OBJECTDIR}/Common.o \
	${OBJECTDIR}/Factor.o \
	${OBJECTDIR}/LearnSPGM.o \
	${OBJECTDIR}/Log.o \
	${OBJECTDIR}/SPGM.o \
	${OBJECTDIR}/SPGM_mixture.o \
	${OBJECTDIR}/SPGMnodes.o \
	${OBJECTDIR}/kmeans-master/kmeans.o \
	${OBJECTDIR}/main1.o \
	${OBJECTDIR}/runs.o \
	${OBJECTDIR}/tests.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-O3
CXXFLAGS=-O3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spgm_0.2

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spgm_0.2: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spgm_0.2 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/ChowLiu.o: ChowLiu.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ChowLiu.o ChowLiu.cpp

${OBJECTDIR}/Common.o: Common.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Common.o Common.cpp

${OBJECTDIR}/Factor.o: Factor.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Factor.o Factor.cpp

${OBJECTDIR}/LearnSPGM.o: LearnSPGM.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LearnSPGM.o LearnSPGM.cpp

${OBJECTDIR}/Log.o: Log.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Log.o Log.cpp

${OBJECTDIR}/SPGM.o: SPGM.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SPGM.o SPGM.cpp

${OBJECTDIR}/SPGM_mixture.o: SPGM_mixture.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SPGM_mixture.o SPGM_mixture.cpp

${OBJECTDIR}/SPGMnodes.o: SPGMnodes.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SPGMnodes.o SPGMnodes.cpp

${OBJECTDIR}/kmeans-master/kmeans.o: kmeans-master/kmeans.cpp 
	${MKDIR} -p ${OBJECTDIR}/kmeans-master
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/kmeans-master/kmeans.o kmeans-master/kmeans.cpp

${OBJECTDIR}/main1.o: main1.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main1.o main1.cpp

${OBJECTDIR}/runs.o: runs.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/runs.o runs.cpp

${OBJECTDIR}/tests.o: tests.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -DNDEBUG -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tests.o tests.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/spgm_0.2

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
