#include <cstring>
#include <vector>
#include "CLI11.hpp"

std::string input;
std::string gcurve;
std::string scurve;
std::string output;
std::string gmap;
std::string map;
std::string outer;
std::string corr;
bool karea1 = false;
bool karea2 = false;
double rad = 30;
double maxArea = 316;
double intvArea = 0;
double speed = 0.2;
double refCortexArea = 0;
double refHullArea = 77100;
int nThreads = 0;

void PARSE_ARGS(int argc, char **argv)
{
    
    std::string desc("Local Gyrification Index "
					 LGI_VERSION "\n"
					 "Author: Ilwoo Lyu\n"
					 "Please refer to the following papers for details:\n"
					 "[1] Lyu et al., A Cortical Shape-Adaptive Approach to Local Gyrification Index, Medical Image Analysis, 48, 244-258, 2018.\n"
					 "[2] Lyu et al. Novel Local Shape-Adaptive Gyrification Index with Application to Brain Development, LNCS10433, 31-39, MICCAI 2017.\n"
					 );

    CLI::App app(desc);

	app.add_option("-i,--inputSurface", input, "Specify input surface")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("--outer", outer, "Specify outer surface - cerebral hull")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-g,--gyrus", gcurve, "Specify input gyral region: *.gcurve|*.gcurve.bary")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-s,--sulcus", scurve, "Specify input sulcal region: *.scurve|*.scurve.bary")->required()->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("-c,--corr", corr, "Specify correspondence between outer and input")->check(CLI::ExistingFile)->group("Inputs");
	app.add_option("--inputMap", map, "Specify input velocity map")->check(CLI::ExistingFile)->group("Inputs");

	app.add_option("-o,--output", output, "Specify output file name without extension")->required()->group("Outputs");
	app.add_option("--geodesicMap", gmap, "Write output geodesic distance map")->group("Outputs");
	app.add_flag("--kArea1", karea1, "Write kernel area of the surface per vertex")->group("Outputs");
	app.add_flag("--kArea2", karea2, "Write kernel area of the hull per vertex")->group("Outputs");
	app.add_option("-t,--intvArea", intvArea, "Write outputs per area interval", true)->check(CLI::NonNegativeNumber)->group("Outputs");

	app.add_option("-r,--radius", rad, "Specify a search radius for the closest sulcal/gyral regions", true)->check(CLI::PositiveNumber)->group("Kernel parameters");
	app.add_option("-m,--maxArea", maxArea, "Specify maximum propagation area (rho)", true)->check(CLI::PositiveNumber)->group("Kernel parameters");
	app.add_option("--speed", speed, "Specify speed at sulcal/gyral points (eta)", true)->check(CLI::PositiveNumber)->group("Kernel parameters");
	app.add_option("--refCortexArea", refCortexArea, "Specify expected reference cortex area", true)->check(CLI::PositiveNumber)->group("Kernel parameters");
	app.add_option("--refHullArea", refHullArea, "Specify expected reference cerebral hull area", true)->check(CLI::PositiveNumber)->group("Kernel parameters");

	app.add_option("--nThreads", nThreads, "Specify the nubmer of OpenMP threads")->check(CLI::NonNegativeNumber)->group("Multi-threading");

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		exit(app.exit(e));
	}
}

