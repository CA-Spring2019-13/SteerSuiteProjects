# GNU Make project makefile autogenerated by Premake
ifndef config
  config=debug
endif

ifndef verbose
  SILENT = @
endif

ifndef CC
  CC = gcc
endif

ifndef CXX
  CXX = g++
endif

ifndef AR
  AR = ar
endif

ifeq ($(config),debug)
  OBJDIR     = obj/Debug/steersimlib
  TARGETDIR  = ../lib
  TARGET     = $(TARGETDIR)/libsteersimlib.so
  DEFINES   += -DENABLE_GUI -DENABLE_GLFW -DDEBUG
  INCLUDES  += -I../../steerlib/include -I../../steersimlib/include -I../../external -I../../util/include
  CPPFLAGS  += -MMD -MP $(DEFINES) $(INCLUDES)
  CFLAGS    += $(CPPFLAGS) $(ARCH) -Wall -g -fPIC -std=c++0x -ggdb `pkg-config --cflags gl` `pkg-config --cflags glu` -fPIC
  CXXFLAGS  += $(CFLAGS) 
  LDFLAGS   += -shared -Wl,-rpath,/home/krishna/Documents/CA/SteerSuiteProjects/A3/build/lib `pkg-config --libs gl` `pkg-config --libs glu` -fPIC -L../lib
  LIBS      += -lsteerlib -lutil -lglfw -lXrandr -lX11 -ldl -lpthread -ltinyxml
  RESFLAGS  += $(DEFINES) $(INCLUDES) 
  LDDEPS    += ../lib/libsteerlib.so ../lib/libutil.so ../lib/libglfw.so ../lib/libtinyxml.so
  LINKCMD    = $(CXX) -o $(TARGET) $(OBJECTS) $(LDFLAGS) $(RESOURCES) $(ARCH) $(LIBS)
  define PREBUILDCMDS
  endef
  define PRELINKCMDS
  endef
  define POSTBUILDCMDS
  endef
endif

ifeq ($(config),release)
  OBJDIR     = obj/Release/steersimlib
  TARGETDIR  = ../lib
  TARGET     = $(TARGETDIR)/libsteersimlib.so
  DEFINES   += -DENABLE_GUI -DENABLE_GLFW -DNDEBUG
  INCLUDES  += -I../../steerlib/include -I../../steersimlib/include -I../../external -I../../util/include
  CPPFLAGS  += -MMD -MP $(DEFINES) $(INCLUDES)
  CFLAGS    += $(CPPFLAGS) $(ARCH) -Wall -g -O2 -fPIC -std=c++0x -ggdb `pkg-config --cflags gl` `pkg-config --cflags glu` -fPIC
  CXXFLAGS  += $(CFLAGS) 
  LDFLAGS   += -shared -Wl,-rpath,/home/krishna/Documents/CA/SteerSuiteProjects/A3/build/lib `pkg-config --libs gl` `pkg-config --libs glu` -fPIC -L../lib
  LIBS      += -lsteerlib -lutil -lglfw -lXrandr -lX11 -ldl -lpthread -ltinyxml
  RESFLAGS  += $(DEFINES) $(INCLUDES) 
  LDDEPS    += ../lib/libsteerlib.so ../lib/libutil.so ../lib/libglfw.so ../lib/libtinyxml.so
  LINKCMD    = $(CXX) -o $(TARGET) $(OBJECTS) $(LDFLAGS) $(RESOURCES) $(ARCH) $(LIBS)
  define PREBUILDCMDS
  endef
  define PRELINKCMDS
  endef
  define POSTBUILDCMDS
  endef
endif

OBJECTS := \
	$(OBJDIR)/ModuleManagerWidget.o \
	$(OBJDIR)/ClockWidget.o \
	$(OBJDIR)/GLWidget.o \
	$(OBJDIR)/ConsoleWidget.o \
	$(OBJDIR)/RecFilePlayerWidget.o \
	$(OBJDIR)/QtEngineDriver.o \
	$(OBJDIR)/QtEngineController.o \
	$(OBJDIR)/GLFWEngineDriver.o \
	$(OBJDIR)/SteerSim.o \
	$(OBJDIR)/CommandLineEngineDriver.o \
	$(OBJDIR)/TestCasePlayerWidget.o \
	$(OBJDIR)/GlobalEventFilter.o \

RESOURCES := \

SHELLTYPE := msdos
ifeq (,$(ComSpec)$(COMSPEC))
  SHELLTYPE := posix
endif
ifeq (/bin,$(findstring /bin,$(SHELL)))
  SHELLTYPE := posix
endif

.PHONY: clean prebuild prelink

all: $(TARGETDIR) $(OBJDIR) prebuild prelink $(TARGET)
	@:

$(TARGET): $(GCH) $(OBJECTS) $(LDDEPS) $(RESOURCES)
	@echo Linking steersimlib
	$(SILENT) $(LINKCMD)
	$(POSTBUILDCMDS)

$(TARGETDIR):
	@echo Creating $(TARGETDIR)
ifeq (posix,$(SHELLTYPE))
	$(SILENT) mkdir -p $(TARGETDIR)
else
	$(SILENT) mkdir $(subst /,\\,$(TARGETDIR))
endif

$(OBJDIR):
	@echo Creating $(OBJDIR)
ifeq (posix,$(SHELLTYPE))
	$(SILENT) mkdir -p $(OBJDIR)
else
	$(SILENT) mkdir $(subst /,\\,$(OBJDIR))
endif

clean:
	@echo Cleaning steersimlib
ifeq (posix,$(SHELLTYPE))
	$(SILENT) rm -f  $(TARGET)
	$(SILENT) rm -rf $(OBJDIR)
else
	$(SILENT) if exist $(subst /,\\,$(TARGET)) del $(subst /,\\,$(TARGET))
	$(SILENT) if exist $(subst /,\\,$(OBJDIR)) rmdir /s /q $(subst /,\\,$(OBJDIR))
endif

prebuild:
	$(PREBUILDCMDS)

prelink:
	$(PRELINKCMDS)

ifneq (,$(PCH))
$(GCH): $(PCH)
	@echo $(notdir $<)
	-$(SILENT) cp $< $(OBJDIR)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
endif

$(OBJDIR)/ModuleManagerWidget.o: ../../steersimlib/src/ModuleManagerWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/ClockWidget.o: ../../steersimlib/src/ClockWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/GLWidget.o: ../../steersimlib/src/GLWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/ConsoleWidget.o: ../../steersimlib/src/ConsoleWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/RecFilePlayerWidget.o: ../../steersimlib/src/RecFilePlayerWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/QtEngineDriver.o: ../../steersimlib/src/QtEngineDriver.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/QtEngineController.o: ../../steersimlib/src/QtEngineController.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/GLFWEngineDriver.o: ../../steersimlib/src/GLFWEngineDriver.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/SteerSim.o: ../../steersimlib/src/SteerSim.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/CommandLineEngineDriver.o: ../../steersimlib/src/CommandLineEngineDriver.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/TestCasePlayerWidget.o: ../../steersimlib/src/TestCasePlayerWidget.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"
$(OBJDIR)/GlobalEventFilter.o: ../../steersimlib/src/GlobalEventFilter.cpp
	@echo $(notdir $<)
	$(SILENT) $(CXX) $(CXXFLAGS) -o "$@" -c "$<"

-include $(OBJECTS:%.o=%.d)
