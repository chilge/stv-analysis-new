all: analyzer univmake

univmake: univmake.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer_chilgenb.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

.PHONY: clean

clean:
	$(RM) univmake analyzer
