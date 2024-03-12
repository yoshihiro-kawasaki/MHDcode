PROGNAME := briowu
# PROGNAME := OrszagTang
OUTDIR := build
SRCDIR := src
EOSDIR := eos
RECONDIR := reconstruction
PROBDIR := problem
TARGET := $(OUTDIR)/$(PROGNAME)
SRCS := $(wildcard $(SRCDIR)/*.cpp) \
        $(wildcard $(SRCDIR)/$(EOSDIR)/*.cpp) \
        $(wildcard $(SRCDIR)/$(RECONDIR)/*.cpp) \
		$(SRCDIR)/$(PROBDIR)/$(PROGNAME).cpp
OBJS := $(addprefix $(OUTDIR)/,$(patsubst %.cpp,%.o,$(SRCS)))

CPP = g++
# CFLAGS = -Wall -g -c -O3 -std=c++17
CFLAGS = -g -c -O3 -std=c++17

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) -o $@ $^

$(OUTDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CPP) $(CFLAGS) -o $@ -c $<

# $(OUTDIR)/$(EOSDIR)/%.o: $(SRCDIR)/$(EOSDIR)/%.cpp
# 	@mkdir -p $(dir $@)
# 	$(CPP) $(CFLAGS) -o $@ -c $<

# $(OUTDIR)/$(RECONDIR)/%.o: $(SRCDIR)/$(RECONDIR)/%.cpp
# 	@mkdir -p $(dir $@)
# 	$(CPP) $(CFLAGS) -o $@ -c $<

# $(OUTDIR)/$(RECONDIR)/%.o: $(SRCDIR)/$(RECONDIR)/%.cpp
# 	@mkdir -p $(dir $@)
# 	$(CPP) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf $(OUTDIR) $(TARGET)

# PROGNAME := mhd
# OUTDIR := build
# SRCDIR := src
# EOSDIR := eos 
# RECONDIR := reconstruction
# TARGET := $(OUTDIR)/$(PROGNAME)
# SRCS := $(wildcard $(SRCDIR)/*.cpp) \
# 		$(wildcard $(SRCDIR)/$(EOSDIR)/*.cpp) \
# 		$(wildcard $(SRCDIR)/$(RECONDIR)/*.cpp)
# OBJS := $(addprefix $(OUTDIR)/, $(patsubst %.cpp,%.o,$(SRCS)))
# # OBJS := $(addprefix $(OUTDIR)/, $(notdir $(SRCS:.cpp=.o)))

# CPP = g++
# CFLAGS = -Wall -g -c -O3 -std=c++17

# .PHONY: all clean
# all: $(TARGET)

# $(TARGET): $(OBJS)
# 	$(CPP) $(CFLAGS) -o $@ $^

# $(OUTDIR)/%.o:%.cpp
# 	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
# 	$(CPP) $(CFLAGS) -o $@ -c $<

# # $(OUTDIR)/$(EOSDIR)/%.o: $(SRCDIR)/$(EOSDIR)/%.cpp
# # 	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
# # 	$(CPP) $(CFLAGS) -o $@ -c $<

# $(OUTDIR)/$(EOSDIR)/%.o:%.cpp
# 	@echo `dirname $@`
# 	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
# 	$(CPP) $(CFLAGS) -o $@ -c $<



# clean:
# 	rm -rf $(OUTDIR)