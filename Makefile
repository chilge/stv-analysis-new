all: analyzer analyzer_3d univmake

univmake: univmake.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer_3d: analyzer_3d.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

.PHONY: clean

clean:
	$(RM) univmake analyzer analyzer_3d
