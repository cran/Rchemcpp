
CXXFLAGS += -w -fPIC -I../src/

objects := constant.o node.o atom.o bond.o cerror.o chemcpp.o datacontainer.o descriptor.o elements.o jlpioutils.o kcfmolecule.o molecule.o moleculeset.o moleculeutils.o ring.o stringutils.o

all: ../../libchemcpp.dll ../../libchemcpp.so

../../libchemcpp.dll: $(objects) 
	$(CXX) -shared -o ../../libchemcpp.dll $(CXXFLAGS) $(objects)

UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
../../libchemcpp.so: $(objects) 
	$(CXX) -dynamiclib -o ../../libchemcpp.so $(CXXFLAGS) $(objects) -Wl,-install_name -Wl,@loader_path/libchemcpp.so
else 
../../libchemcpp.so: $(objects) 
	$(CXX) -shared -o ../../libchemcpp.so $(CXXFLAGS) $(objects)
endif
