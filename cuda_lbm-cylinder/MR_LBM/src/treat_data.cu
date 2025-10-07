#include "treat_data.cuh"

__host__ void calculate_pressure(cylinderProperties *h_cylinder_properties, unsigned int count, unsigned int step)
{
	if (step == STAT_BEGIN_TIME)
	{
		std::ostringstream filename_stream;
		filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "theta_" << ID_SIM << ".dat";
		std::string filename = filename_stream.str();

		std::ofstream data_file(filename.c_str(), std::ios::app);

		data_file << count;
		for (int i = 0; i < count; i++)
		{
			cylinderProperties property = h_cylinder_properties[i];

			data_file << std::setprecision(10) << " " << property.theta;
		}

		data_file << std::endl;
		data_file.close();
	}

	std::ostringstream filename_stream;
	filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "pressure_" << ID_SIM << ".dat";
	std::string filename = filename_stream.str();

	std::ofstream data_file(filename.c_str(), std::ios::app);

	data_file << step;

	for (int i = 0; i < count; i++)
	{
		cylinderProperties property = h_cylinder_properties[i];

		data_file << std::setprecision(10) << " " << property.ps;
	}
	data_file << std::endl;

	data_file.close();
}

__host__ void calculate_forces(cylinderProperties *h_cylinder_properties, unsigned int count, unsigned int step)
{
	dfloat f_x_net = 0.0;
	dfloat f_y_net = 0.0;

	for (int i = 0; i < count; i++)
	{
		f_x_net += h_cylinder_properties[i].Fx;
		f_y_net += h_cylinder_properties[i].Fy;
	}

	// write to a file
	std::ostringstream filename_stream;
	filename_stream << PATH_FILES << "/" << ID_SIM << "/" << "forces_" << ID_SIM << ".dat";
	std::string filename = filename_stream.str();

	std::ofstream data_file(filename.c_str(), std::ios::app);

	data_file << std::setprecision(10) << std::fixed << step << " " << f_x_net << " " << f_y_net << std::endl;
	data_file.close();
}

__host__ void calculate_inlet_density(dfloat *h_fMom, unsigned int step, dfloat* rho_infty)
{
	dfloat rho_inlet = 0;

	for (int y = 0; y < NY; y++)
	{
		rho_inlet += RHO_0 + h_fMom[idxMom(0 % BLOCK_NX, y % BLOCK_NY, M_RHO_INDEX, 0 / BLOCK_NX, y / BLOCK_NY)];
	}

	rho_inlet /= NY;

	*rho_infty = rho_inlet;

	std::ostringstream filename_rho;
	filename_rho << PATH_FILES << "/" << ID_SIM << "/" << "rho_inlet_" << ID_SIM << ".dat";
	std::string filename_r = filename_rho.str();

	std::ofstream rho_file(filename_r.c_str(), std::ios::app);

	rho_file << std::setprecision(10) << step << " " << rho_inlet << std::endl;
	rho_file.close();
}

__global__ void domain_avg(dfloat *fMom, dfloat *ux_mean, dfloat *uy_mean, unsigned int step)
{
	unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;

	unsigned int global_index = x + y * NX;

	int meanCounter = step - INI_MEAN_STEP;
	dfloat invCount = 1.0 / (meanCounter + 1.0);

	dfloat ux = fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UX_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;
	dfloat uy = fMom[idxMom(x % BLOCK_NX, y % BLOCK_NY, M_UY_INDEX, x / BLOCK_NX, y / BLOCK_NY)] / F_M_I_SCALE;

	ux_mean[global_index] = (ux_mean[global_index] * meanCounter + ux) * invCount;
	uy_mean[global_index] = (uy_mean[global_index] * meanCounter + uy) * invCount;
}

__global__ void velocity_on_centerline_average(dfloat *fMom, dfloat *ux_center, dfloat *uy_center, unsigned int step)
{
	const size_t offset = L_front + D;
	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	const size_t x = idx + offset;

	if(idx >= L_back) return;

	const size_t top_y_coord = NY / 2;
	const size_t bot_y_coord = top_y_coord - 1;

	const size_t global_x_thread = x % BLOCK_NX;
	const size_t global_top_y_thread = top_y_coord % BLOCK_NY;
	const size_t global_bot_y_thread = bot_y_coord % BLOCK_NY;

	const size_t global_x_block = x / BLOCK_NX;
	const size_t global_top_y_block = top_y_coord / BLOCK_NY;
	const size_t global_bot_y_block = bot_y_coord / BLOCK_NY;

	const dfloat ux_top = fMom[idxMom(global_x_thread, global_top_y_thread, M_UX_INDEX, global_x_block, global_top_y_block)] / F_M_I_SCALE;
	const dfloat ux_bot = fMom[idxMom(global_x_thread, global_bot_y_thread, M_UX_INDEX, global_x_block, global_bot_y_block)] / F_M_I_SCALE;

	const dfloat uy_top = fMom[idxMom(global_x_thread, global_top_y_thread, M_UY_INDEX, global_x_block, global_top_y_block)] / F_M_I_SCALE;
	const dfloat uy_bot = fMom[idxMom(global_x_thread, global_bot_y_thread, M_UY_INDEX, global_x_block, global_bot_y_block)] / F_M_I_SCALE;

	const dfloat ux = (ux_top + ux_bot) * 0.5;
	const dfloat uy = (uy_top + uy_bot) * 0.5;

	const size_t time_count = step - STAT_BEGIN_TIME;
	const dfloat inv_count = 1.0f / (1.0f + time_count);

	ux_center[idx] = (ux_center[idx] * time_count + ux) * inv_count;
	uy_center[idx] = (uy_center[idx] * time_count + uy) * inv_count;
}