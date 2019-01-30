#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int main() {
	int i;
	istringstream ss;
	ifstream inctl("test.txt",ios::in);
	char opt[1024], null[1024];

	//inctl>>opt;
	while (!inctl.eof()) {
		inctl.getline(opt,1024);
		printf("%s\n",opt);
	}
}

