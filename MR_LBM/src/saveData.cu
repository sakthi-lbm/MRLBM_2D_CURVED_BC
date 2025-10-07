#include "saveData.cuh"
#include "globalFunctions.h"
#include <iostream>

__host__ void saveMacr(
	dfloat *h_fMom, dfloat *rho, dfloat *ux, dfloat *uy, unsigned int nSteps)
{
	int x, y;

	// linearize
	size_t indexMacr;
	// double uSum = 0;
	// double t_star = 0;

	// printf("\n--------------------------- Save macr %d ---------------------------\n", step);

	for (y = 0; y < NY; y++)
	{
		for (x = 0; x < NX; x++)
		{
			indexMacr = idxScalarGlobal(x, y);

			rho[indexMacr] = RHO_0 + h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_RHO_INDEX, x / BLOCK_NX, y / BLOCK_NY)];
			ux[indexMacr] = h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UX_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;
			uy[indexMacr] = h_fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UY_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;
		}
	}

	// Sakthi-modifications
	// Creating master.p3d file
	// ======================================================================================================================
	std::string strInf = PATH_FILES;
	strInf += "/";
	strInf += ID_SIM;
	strInf += "/";
	strInf += "master.p3d"; // generate file name (with path)
	std::string prefix = "data";
	std::string suffix = ".f";
	std::ostringstream master_file;
	master_file << strInf;
	std::ofstream out(master_file.str());
	out << "{" << std::endl;
	;
	out << std::endl;
	;
	out << " \"auto-detect-format\": true," << std::endl;
	;
	out << std::endl;
	;
	out << " \"filenames\": [" << std::endl;
	;
	out << std::endl;
	;

	for (int iter = 1; iter < N_STEPS; iter++)
	{
		if (iter % MACR_SAVE == 0)
		{
			std::ostringstream filename_temp;
			filename_temp << prefix << (10000000 + iter) << suffix;
			std::string filename = filename_temp.str();

			if (out.is_open())
			{
				out << "{ \"time\" :  " << iter / MACR_SAVE << ", \"xyz\" : \"grid.x\", \"function\" : \"" << filename << "\" }," << std::endl;
				;
			}
			else
			{
				std::cerr << "Failed to open file: master.p3d" << std::endl;
			}
		}
	}
	out << std::endl;
	;
	out << "]" << std::endl;
	;
	out << "}" << std::endl;
	;
	out.close();
	// ======================================================================================================================

	// writing grid.x file
	// -------------------------------------------------------------------------------------------------------
	int nprocs = 1;
	if (nSteps == 0)
	{
		std::string strInf2 = PATH_FILES;
		strInf2 += "/";
		strInf2 += ID_SIM;
		strInf2 += "/";
		strInf2 += "grid.x"; // generate file name (with path)

		std::ofstream gridfile(strInf2, std::ios::binary);
		if (!gridfile)
		{
			std::cerr << "Error opening grid file" << std::endl;
		}

		// Write nprocs
		gridfile.write(reinterpret_cast<const char *>(&nprocs), sizeof(int));

		// Write (nx, ny) for each processor (Fortran loop: m = 1 to nprocs)
		for (int m = 1; m <= nprocs; ++m)
		{
			gridfile.write(reinterpret_cast<const char *>(&NX), sizeof(int));
			gridfile.write(reinterpret_cast<const char *>(&NY), sizeof(int));
		}

		// Write x and y arrays for each processor (m = 0 to nprocs - 1)
		for (int m = 0; m < nprocs; ++m)
		{
			// Fortran is column-major: loop j outer, i inner
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					float val = (double(i) - xc) / D; // already float
					gridfile.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}

			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
				{
					float val = (double(j) - yc) / D;
					gridfile.write(reinterpret_cast<const char *>(&val), sizeof(float));
				}
		}
	}

	// datafile

	std::ostringstream filename_temp;

	std::string strInf3 = PATH_FILES;
	strInf3 += "/";
	strInf3 += ID_SIM;
	strInf3 += "/";
	filename_temp << strInf3 << prefix << (10000000 + nSteps) << suffix;
	std::string filename = filename_temp.str();

	std::ofstream datafile(filename, std::ios::binary);
	if (!datafile)
	{
		std::cerr << "Error opening grid file" << std::endl;
	}

	// Write nprocs
	datafile.write(reinterpret_cast<const char *>(&nprocs), sizeof(int));

	// Write (nx, ny) for each processor (Fortran loop: m = 1 to nprocs)
	for (int l = 0; l < nprocs; ++l)
	{
		datafile.write(reinterpret_cast<const char *>(&NX), sizeof(int));
		datafile.write(reinterpret_cast<const char *>(&NY), sizeof(int));
		int nf = 3; // 3 fields: rho, ux, uy
		datafile.write(reinterpret_cast<const char *>(&nf), sizeof(int));
	}

	// Write x and y arrays for each processor (m = 0 to nprocs - 1)
	for (int m = 0; m < nprocs; ++m)
	{
		// Fortran is column-major: loop j outer, i inner
		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
			{
				indexMacr = idxScalarGlobal(i, j);
				float val = rho[indexMacr]; // already float
				datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
			}

		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
			{
				indexMacr = idxScalarGlobal(i, j);
				float val = ux[indexMacr];
				datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
			}

		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
			{
				indexMacr = idxScalarGlobal(i, j);
				float val = uy[indexMacr];
				datafile.write(reinterpret_cast<const char *>(&val), sizeof(float));
			}
	}
	datafile.close();

	std::string strFileRho, strFileUx, strFileUy;

	strFileRho = getVarFilename("rho", nSteps, ".bin");
	strFileUx = getVarFilename("ux", nSteps, ".bin");
	strFileUy = getVarFilename("uy", nSteps, ".bin");

	// saveVarBin(strFileRho, rho, MEM_SIZE_SCALAR);
	// saveVarBin(strFileUx, ux, MEM_SIZE_SCALAR);
	// saveVarBin(strFileUy, uy, MEM_SIZE_SCALAR);
}

std::string getVarFilename(
	const std::string varName,
	unsigned int step,
	const std::string ext)
{
	unsigned int n_zeros = 0, pot_10 = 10;
	unsigned int aux1 = 1000000; // 6 numbers on step
	// calculate number of zeros
	if (step != 0)
		for (n_zeros = 0; step * pot_10 < aux1; pot_10 *= 10)
			n_zeros++;
	else
		n_zeros = 5;

	// generates the file name as "PATH_FILES/id/id_varName000000.bin"
	std::string strFile = PATH_FILES;
	strFile += "/";
	strFile += ID_SIM;
	strFile += "/";
	strFile += ID_SIM;
	strFile += "_";
	strFile += varName;
	for (unsigned int i = 0; i < n_zeros; i++)
		strFile += "0";
	strFile += std::to_string(step);
	strFile += ext;

	return strFile;
}

void saveVarBin(
	std::string strFile,
	dfloat *var,
	size_t memSize)
{
	FILE *outFile = nullptr;

	outFile = fopen(strFile.c_str(), "wb");

	if (outFile != nullptr)
	{
		fwrite(var, memSize, 1, outFile);
		fclose(outFile);
	}
	else
	{
		printf("Error saving \"%s\" \nProbably wrong path!\n", strFile.c_str());
	}
}

std::string getSimInfoString(int step, dfloat MLUPS, int D, dfloat Dcy, size_t count, dfloat rho_infty)
{
#define BOOL_TO_STRING(b) ((b) ? "true" : "false")
	std::ostringstream strSimInfo("");

	strSimInfo << std::fixed;
	strSimInfo << std::setprecision(6);

	strSimInfo << "---------------------------- SIMULATION INFORMATION ----------------------------\n";
	strSimInfo << "      Simulation ID =  " << ID_SIM << "\n";
	strSimInfo << "       Velocity set = D2Q9\n";
	strSimInfo << "                 Re = " << RE << "\n";
	strSimInfo << "          Precision = float\n";
	strSimInfo << "                 NX = " << NX << "\n";
	strSimInfo << "                 NY = " << NY << "\n";
	strSimInfo << std::fixed << std::setprecision(6);
	/*strSimInfo << "                Tau: " << TAU << "\n";*/
	strSimInfo << "                 uo = " << U_MAX << "\n";
	strSimInfo << "          Macr_save = " << MACR_SAVE << "\n";
	strSimInfo << "             Nsteps = " << step << "\n";
	strSimInfo << "              MLUPS = " << MLUPS << "\n";
	strSimInfo << std::fixed << std::setprecision(0);
	strSimInfo << "                 BX = " << BLOCK_NX << "\n";
	strSimInfo << "                 BY = " << BLOCK_NY << "\n";
	strSimInfo << "                ncy = " << count << "\n";
	strSimInfo << "                  D = " << D << "\n";
	strSimInfo << std::fixed << std::setprecision(5);
	strSimInfo << "                Dcy = " << Dcy << "\n";
	strSimInfo << "          rho_infty = " << rho_infty << "\n";
	strSimInfo << "--------------------------------------------------------------------------------\n";

	return strSimInfo.str();
}

void folderSetup()
{
	// Windows
#if defined(_WIN32)
	std::string strPath;
	strPath = PATH_FILES;
	strPath += "\\\\"; // adds "\\"
	strPath += ID_SIM;
	std::string cmd = "md ";
	cmd += strPath;
	system(cmd.c_str());
	return;
#endif // !_WIN32

	// Unix
#if defined(__APPLE__) || defined(__MACH__) || defined(__linux__)
	std::string strPath;
	strPath = PATH_FILES;
	strPath += "/";
	strPath += ID_SIM;
	strPath += "/";
	strPath += "plots";
	std::string cmd = "mkdir -p ";
	cmd += strPath;
	const int i = system(cmd.c_str());
	static_cast<void>(i);
	return;
#endif // !Unix
	printf("I don't know how to setup folders for your operational system :(\n");
	return;
}

void saveSimInfo(int step, dfloat MLUPS, int D, dfloat Dcy, size_t count, dfloat rho_infty)
{
	std::string strInf = PATH_FILES;
	strInf += "/";
	strInf += ID_SIM;
	strInf += "/";
	strInf += ID_SIM;
	strInf += "_info.txt"; // generate file name (with path)
	FILE *outFile = nullptr;

	outFile = fopen(strInf.c_str(), "w");
	if (outFile != nullptr)
	{
		std::string strSimInfo = getSimInfoString(step, MLUPS, D, Dcy, count, rho_infty);
		fprintf(outFile, "%s\n", strSimInfo.c_str());
		fclose(outFile);
	}
	else
	{
		printf("Error saving \"%s\" \nProbably wrong path!\n", strInf.c_str());
	}
}
/**/

__host__ void saving_centerline_data(dfloat *ux, dfloat *uy)
{
	std::ostringstream source_path;
	source_path << PATH_FILES << "/" << ID_SIM << "/" << "centerline_" << ID_SIM << ".dat";
	std::ofstream out_file(source_path.str());

	out_file << "x_coord" << " ux " << " uy" << std::endl;

	for (size_t i = 0; i < L_back; ++i)
	{
		out_file << (i + L_front + D) << " " << ux[i] << " " << uy[i] << std::endl;
	}
}