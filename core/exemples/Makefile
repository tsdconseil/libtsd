# Exemple de Makefile minimal pour utiliser libtsd

LDFLAGS  =  -ltsd -lfreetype -lpng -lz -lfmt -L/usr/local/lib -fopenmp 

# Changer éventuellement le dossier "-I../core/include"
CXXFLAGS = -g -march=native -fPIC -std=c++2a `pkg-config --cflags eigen3` -I/usr/local/include `pkg-config freetype2 --cflags` 


all: 
	g++ exemple.cc $(CXXFLAGS) $(LDFLAGS) -o exemple.exe
	g++ example.cc $(CXXFLAGS) $(LDFLAGS) -o example.exe




