build:
	rm -f application.exe
	g++ -std=c++11 -Wall -g application.cpp dist.cpp osm.cpp tinyxml2.cpp -o application.exe

run:
	./application.exe

buildtest:
	rm -f testing.exe
	g++ -std=c++11 -Wall -g testing.cpp -o testing.exe

runtest:
	./testing.exe
	

valgrind:
	valgrind --tool=memcheck --leak-check=yes ./application.exe
