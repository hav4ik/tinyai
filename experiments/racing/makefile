TARGET = racing

CC = g++
CFLAGS = -O2 -std=c++11

INCAI = -I../../src
INCXML = -ITINY

INCSFML = -ISFML/include
LIBSFML = -LSFML/lib

INC = $(INCAI) $(INCXML) $(INCSFML)
LIB = $(LIBSFML) -lsfml-graphics -lsfml-system -lsfml-window


SOURCES = \
	main.cpp \
	level.cpp \
	TINY/tinyxmlerror.cpp \
    TINY/tinyxml.cpp \
    TINY/tinystr.cpp \
    TINY/tinyxmlparser.cpp \
    TINY/tinyxml2.cpp 

HEADERS = \
	level.h \
	physics.h \
	TINY/tinyxml.h \
    TINY/tinystr.h \
    TINY/tinyxml2.h 


platformer:
	$(CC) $(CFLAGS) $(INC) $(LIB) $(SOURCES) -o $(TARGET)

clean:
	rm racing