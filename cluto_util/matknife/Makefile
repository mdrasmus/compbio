APP=matknife
CC = g++
OBJS = \
        matknife.o 

COMPILE = -c -g -Wall
LINK = 

all: $(APP)

# basic compile rule
$(OBJS): %.o: %.cpp
	$(CC) $(COMPILE) -o $@ $<

$(APP): $(OBJS)
	$(CC) $(OBJS) $(LINK) -o $(APP)

clean:
	rm -rf $(APP) $(OBJS)
