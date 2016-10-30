CC=clang++
CFLAGS+= -std=c++11
PROD=bezEditor
SRCS= $(shell ls *.cpp)
OBJS= $(patsubst %.cpp, %.o, $(SRCS))
STRINGFIX = -ldl -lm -lXext -lX11 -lXxf86vm -lXrandr -lXi -lXcursor -lXinerama -L/usr/lib/nvidia-340-updates
GLFW = -lglfw3 -lGLU -lGL -lpthread
LIBS = $(STRINGFIX) $(GLFW)
all: $(PROD)

$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(PROD): $(OBJS)
	$(CC) -o $(PROD) $^ $(LIBS)

clean:
	rm -f $(PROD)
	rm -f $(OBJS)
