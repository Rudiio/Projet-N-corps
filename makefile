CC=g++ 
CFLAGS = -Wall -ggdb3 -O3
SRC = ./src/*.c
OBJ = BHTree.o IIntegrator.o IModel.o IntegratorEuler.o ModelNBody.o SDLWnd.o NBodyWnd.o Types.o Vector.o gettimes.o IntegratorLeapFrog.o 
FLAGS = -L /lib64 
LIBS = -l pthread -fopenmp 
LIBS_FLAGS = -lSDL -lGL -lGLU -lglut -lX11 -lSDL_gfx
EXEC = main
EXEC_MESURE = main_mesure
EXEC_TESTS = tests 
# -lusb-1.0
	
.PHONY : $(EXEC) $(EXEC_MESURE) $(EXEC_TESTS) $(OBJ) mv_object clean

all : $(EXEC) $(EXEC_MESURE) $(EXEC_TESTS)

BHTree.o : ./src/BHTree.cpp ./src/BHTree.h
	$(CC) $(CFLAGS) -c $< 

IIntegrator.o : ./src/IIntegrator.cpp ./src/IIntegrator.h
	$(CC) $(CFLAGS) -c $< 

IModel.o : ./src/IModel.cpp ./src/IModel.h
	$(CC) $(CFLAGS) -c $< 

IntegratorEuler.o :  ./src/IntegratorEuler.cpp ./src/IntegratorEuler.h
	$(CC) $(CFLAGS) -c $< 

IntegratorLeapFrog.o :  ./src/IntegratorLeapFrog.cpp ./src/IntegratorLeapFrog.h
	$(CC) $(CFLAGS) -c $< 

ModelNBody.o : ./src/ModelNBody.cpp ./src/ModelNBody.h
	$(CC) $(CFLAGS) $(LIBS) -c $< 

SDLWnd.o : ./src/SDLWnd.cpp ./src/SDLWnd.h
	$(CC) $(CFLAGS) -c $< -lSDL -DOpenGL_GL_PREFERENCE=LEGACY

NBodyWnd.o : ./src/NBodyWnd.cpp ./src/NBodyWnd.h
	$(CC) $(CFLAGS) -c $< -lSDL -DOpenGL_GL_PREFERENCE=LEGACY

Types.o :  ./src/Types.cpp ./src/Types.h
	$(CC) $(CFLAGS) -c $< 

Vector.o : ./src/Vector.cpp ./src/Vector.h
	$(CC) $(CFLAGS) -c $< 

# graphique.pyc: ./src/graphique.py
# 	python3 ./src/graphique.py $<

mv_object :
	mkdir -p OBJ && mv *.o OBJ

$(EXEC) : $(OBJ) mv_object 
	$(CC) $(CFLAGS) $(LIBS) ./src/main.cpp ./OBJ/*.o -o $@ $(LIBS_FLAGS)

$(EXEC_MESURE) :  
	$(CC) $(CFLAGS) $(LIBS) ./src/main_mesure.cpp ./OBJ/*.o -o $@ $(LIBS_FLAGS) 

$(EXEC_TESTS):
	$(CC) $(CFLAGS) ./src/lance_tests.cpp  -o $@ 

clean:
	rm -f prog *.o *.mod ./OBJ/*.o $(EXEC) $(EXEC_MESURE) $(EXEC_TESTS)

# default:
# 	$(GC) -O3 -c BHTree.cpp
# 	$(GC) -O3 -c IIntegrator.cpp
# 	$(GC) -O3 -c IModel.cpp
# 	$(GC) -O3 -c IntegratorEuler.cpp
# 	$(GC) -O3 $(LIBS) -c ModelNBody.cpp 
# 	$(GC) -O3 -c SDLWnd.cpp -lSDL -DOpenGL_GL_PREFERENCE=LEGACY
# 	$(GC) -O3 -c NBodyWnd.cpp -lSDL -DOpenGL_GL_PREFERENCE=LEGACY
# 	$(GC) -O3 -c Types.cpp
# 	$(GC) -O3 -c Vector.cpp
# 	$(GC) -O3 $(LIBS) main.cpp BHTree.o IIntegrator.o IModel.o IntegratorEuler.o ModelNBody.o SDLWnd.o NBodyWnd.o Types.o Vector.o -o main -lSDL -lGL -lGLU -lglut -lX11 -lSDL_gfx