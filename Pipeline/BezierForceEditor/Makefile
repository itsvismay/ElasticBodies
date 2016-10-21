CC=clang++
CFLAGS+= -std=c++11
PROD=bezEditor
SRCS= $(shell ls *.cpp)
OBJS= $(patsubst %.cpp, %.o, $(SRCS))
all: $(PROD)

$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(PROD): $(OBJS)
	$(CC) -o $(PROD) $^ $(LIBS)

clean:
	rm -f $(PROD)
	rm -f $(OBJS)
