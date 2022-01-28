# Exemple de Makefile minimal pour utiliser libtsd

LDFLAGS  =  -ltsd -lfreetype -lpng -lz -lfmt 

# Changer éventuellement le dossier "-I../core/include"
CXXFLAGS = -g -march=native -fPIC -std=c++2a `pkg-config --cflags eigen3`


all: 
	g++ exlibtsd.cc $(CXXFLAGS) $(LDFLAGS)




