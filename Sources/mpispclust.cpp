#include<iostream>
#include<cstring>
using namespace std;

std::string exec(const char* cmd, bool realtimeOutput=false) {
	/*
	This function executes an external shell command and retrieves its output in a string.
	Input:
	-cmd : The command to execute.
	Output:
	-result : A string holding the command output.
	*/
	char buffer[128];
	std::string result = "";
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	FILE* pipe = popen(cmd, "r");
#else
	FILE* pipe = _popen(cmd, "r");
#endif
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
		while (!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL){
				if(realtimeOutput){
					cout<<buffer;
				}
				result += buffer;
			}
		}
	}
	catch (...) {
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
		pclose(pipe);
#else			
		_pclose(pipe);
#endif
		throw;
	}
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	pclose(pipe);
#else			
	_pclose(pipe);
#endif
	return result;
}

bool containsSpaces(char *inStr){
	for(int i=0; i<strlen(inStr); i++)
		if(isspace(inStr[i]))
			return true;
	return false;
}

int main(int argc, char* argv[]){
	char cmd[300];
	string n;
#if !(defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)) //Not under Windows OS
	strcpy(cmd, "mpirun -n ");
	n=exec("nproc --all");
#else
	strcpy(cmd, "mpiexec -n ");
	n=exec("echo %NUMBER_OF_PROCESSORS%");
#endif
	strcat(cmd, n.c_str());
	for(int i=0; i<strlen(cmd); i++)
		if(cmd[i]=='\r'||cmd[i]=='\n')
			cmd[i]='\0';
	strcat(cmd, " spclust");
	for(int i=1; i<argc; i++){
		strcat(cmd, " ");
		if(containsSpaces(argv[i])){
			strcat(cmd, "\"");
			strcat(cmd, argv[i]);
			strcat(cmd, "\"");
		} else
			strcat(cmd, argv[i]);
	}
	//cout<<cmd<<endl;
	//n=exec(cmd);
	//cout<<n.c_str();
	exec(cmd,true);
	return 0;
}