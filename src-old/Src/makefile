##################################
# My makefile. Edit conveniently!
##################################

EXE = 

.SUFFIXES:
.SUFFIXES: $(EXE) .o .h .cpp

#EFENCEPATH = "D:\electric-fence-2.1.13"

OBJS = IDException.o IDOptions.o ioutil.o linalg.o dae.o dyn.o tra.o Matrix.o TPolyn.o matrixseq.o daeIndexDet.o
#time.o

all: daeIndexDet$(EXE)

clean:
	-rm -f $(OBJS) daeIndexDet$(EXE)

.cpp.o:
#	g++ -g -I ${HOME}/adolc_base/include -c $*.cpp
	g++ -I"D:\Temp\temp\adolc_base\include" -c $*.cpp
	     
daeIndexDet$(EXE): $(OBJS)
#	g++  -o $@ $(OBJS) -L ${HOME}/adolc_base/lib -Xlinker -rpath -Xlinker ${HOME}/adolc_base/lib -ladolc -lm
#	g++  -o $@ $(OBJS) ${HOME}/adolc_base/lib/libadolc.a -lm
	g++  -o $@ $(OBJS) -L"D:\Temp\adolc_base\lib" -llibadolc
#-L$(EFENCEPATH) -lefence -lpthread
