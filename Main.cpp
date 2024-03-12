#include "Declarations.h"
#include "Headers/GlobalParameters.h"


using namespace std;

int main()
{
	srand(time(0));
	vector<ofstream> outputs(2);
	outputs[0].open("Outputs/output_u.txt");
	outputs[1].open("Outputs/output_v.txt");
	for(int o=0;o<outputs.size();o++)
		outputs[o].precision(10);

	// simulation of the gray-scott model 
	GrayScottSimulation(&outputs);
	
	for(int o=0;o<outputs.size();o++)
	{	
		outputs[o].flush();
		outputs[o].close();
	}
	return 0;
}