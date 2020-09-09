OBJDIR = Obj
SRCDIR = Src
INCDIR = Inc

HEADERS = $(wildcard $(INCDIR)/*.h)
OBJECTS = $(subst $(INCDIR), $(OBJDIR), $(HEADERS:.h=.o))

CC = g++
CFLAGS = -O2 -std=c++11 `root-config --libs --cflags` \
		 -I$(INCDIR)

LDFLAGS = `root-config --glibs` 

all: $(OBJECTS) analyze
testbeam: $(OBJECTS) analyze


analyze: $(OBJECTS) analyze.cpp 
	@echo "Building executable 'analyze'..."
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/Detector.o: $(SRCDIR)/Detector.cpp $(INCDIR)/Detector.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/DetectorSetup.o: $(SRCDIR)/DetectorSetup.cpp $(INCDIR)/DetectorSetup.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@


$(OBJDIR)/ExpSetup.o: $(SRCDIR)/ExpSetup.cpp $(INCDIR)/ExpSetup.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/FileReader.o: $(SRCDIR)/FileReader.cpp $(INCDIR)/FileReader.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@


$(OBJDIR)/Channel.o: $(SRCDIR)/Channel.cpp $(INCDIR)/Channel.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/global.o: $(SRCDIR)/global.cpp $(INCDIR)/global.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm $(OBJDIR)/*.o
