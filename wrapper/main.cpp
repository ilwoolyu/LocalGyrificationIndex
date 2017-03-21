#include "Gyrification.h"
#include "wrapperCLP.h"

int main(int argc, char *argv[])
{
	PARSE_ARGS;

	if (argc < 2 || (output.empty() && gmap.empty()) || input.empty())
	{
		std::cout << "Usage: " << argv[0] << " " << "--help" << std::endl;
		return EXIT_FAILURE;
	}
	
	Gyrification gi;
	
	const char *cfile = NULL;
	if (!corr.empty()) cfile = corr.c_str();

	if (gpoint.empty())
		gi.open(input.c_str(), scurve.c_str(), gcurve.c_str(), outer.c_str(), cfile);
	else
	{
		gi.open(input.c_str(), gpoint.c_str(), outer.c_str(), cfile);
		if (!scurve.empty()) gi.loadScurve(scurve.c_str());
		else if (!spoint.empty()) gi.loadSpoint(spoint.c_str());
	}

	gi.setPopulationArea(populationArea);
	gi.setOutputFileName(output.c_str());
	gi.setKernelInterval(intvArea);
	gi.setMaxKernelSize(maxArea);
	gi.setSpeed(speed);
	
	if (map.empty()) gi.run(rad);
	else gi.run(map.c_str());

	if (!output.empty()) gi.saveGI(output.c_str());
	if (!gmap.empty()) gi.saveGMap(gmap.c_str());
	if (karea1) gi.saveKMap1(output.c_str());
	if (karea2) gi.saveKMap2(output.c_str());

	return EXIT_SUCCESS;
}
