#include "FluidSim.h"

void FluidSim::solvePoisson() {
	double dx2 = m_dx * m_dx;
	double residual = m_acc + 1; // initial residual
	double rho = 1.0; // density can be eliminated from the equation, so we just set it to 1.0

	Array2d& p = p_pressure->x();
	Array2d& u = p_velocity->x();	// x velocity
	Array2d& v = p_velocity->y();	// y velocity

	for (int it = 0; residual > m_acc && it < m_iter; ++it) {
		// Note that the boundaries are handles by the framework, so you iterations should be similar to:
		for (int y = 1; y < m_res_y - 1; ++y) {
			for (int x = 1; x < m_res_x - 1; ++x) {
				if (p_solidObstacle->x()(x, y) > 0.5) {
					if (p_solidObstacle->x()(x - 1, y) < 0.5) std::get<0>(m_boundaryPressure[{x, y}]) = p(x - 1, y);
					if (p_solidObstacle->x()(x + 1, y) < 0.5) std::get<1>(m_boundaryPressure[{x, y}]) = p(x + 1, y);
					if (p_solidObstacle->x()(x, y - 1) < 0.5) std::get<2>(m_boundaryPressure[{x, y}]) = p(x, y - 1);
					if (p_solidObstacle->x()(x, y + 1) < 0.5) std::get<3>(m_boundaryPressure[{x, y}]) = p(x, y + 1);
				} else {
					double b = -p_divergence->x()(x, y) / m_dt * rho;
					if (p_solidObstacle->x()(x + 1, y) > 0.5) p(x + 1, y) = std::get<0>(m_boundaryPressure[{x + 1, y}]);
					if (p_solidObstacle->x()(x - 1, y) > 0.5) p(x - 1, y) = std::get<1>(m_boundaryPressure[{x - 1, y}]);
					if (p_solidObstacle->x()(x, y + 1) > 0.5) p(x, y + 1) = std::get<2>(m_boundaryPressure[{x, y + 1}]);
					if (p_solidObstacle->x()(x, y - 1) > 0.5) p(x, y - 1) = std::get<3>(m_boundaryPressure[{x, y - 1}]);
					p(x, y) = (p(x + 1, y) + p(x - 1, y) + p(x, y + 1) + p(x, y - 1) + b * dx2) / 4.0;
				}
			}
		}

		// Compute the new residual, i.e. the sum of the squares of the individual residuals (squared L2-norm)
		residual = 0;
		for (int y = 1; y < m_res_y - 1; ++y) {
			for (int x = 1; x < m_res_x - 1; ++x) {
				double cellResidual;
				if (p_solidObstacle->x()(x, y) > 0.5) {
					if (p_solidObstacle->x()(x - 1, y) < 0.5) cellResidual = std::get<0>(m_boundaryPressure[{x, y}]) - p(x - 1, y);
					if (p_solidObstacle->x()(x + 1, y) < 0.5) cellResidual = std::get<1>(m_boundaryPressure[{x, y}]) - p(x + 1, y);
					if (p_solidObstacle->x()(x, y - 1) < 0.5) cellResidual = std::get<2>(m_boundaryPressure[{x, y}]) - p(x, y - 1);
					if (p_solidObstacle->x()(x, y + 1) < 0.5) cellResidual = std::get<3>(m_boundaryPressure[{x, y}]) - p(x, y + 1);
				} else {
					double b = -p_divergence->x()(x, y) / m_dt * rho;
					if (p_solidObstacle->x()(x + 1, y) > 0.5) p(x + 1, y) = std::get<0>(m_boundaryPressure[{x + 1, y}]);
					if (p_solidObstacle->x()(x - 1, y) > 0.5) p(x - 1, y) = std::get<1>(m_boundaryPressure[{x - 1, y}]);
					if (p_solidObstacle->x()(x, y + 1) > 0.5) p(x, y + 1) = std::get<2>(m_boundaryPressure[{x, y + 1}]);
					if (p_solidObstacle->x()(x, y - 1) > 0.5) p(x, y - 1) = std::get<3>(m_boundaryPressure[{x, y - 1}]); 
					cellResidual = p(x, y) - (p(x + 1, y) + p(x - 1, y) + p(x, y + 1) + p(x, y - 1) + b * dx2) / 4.0;
				}
				residual += cellResidual * cellResidual;
			}
		}

		// Get the L2-norm of the residual
		residual = sqrt(residual);

		// We assume the accuracy is meant for the average L2-norm per grid cell
		residual /= (m_res_x - 2) * (m_res_y - 2);

		// For your debugging, and ours, please add these prints after every iteration
		if (it == m_iter - 1)
			printf("Pressure solver: it=%d , res=%f \n", it, residual);
		if (residual < m_acc)
			printf("Pressure solver: it=%d , res=%f, converged \n", it, residual);
	}
}

void FluidSim::correctVelocity() {
	Array2d& p = p_pressure->x();
	Array2d& u = p_velocity->x();	// x velocity
	Array2d& v = p_velocity->y();	// y velocity
	double rho = 1.0; // density can be eliminated from the equation, so we just set it to 1.0

	// Note: velocity u_{i+1/2}
	for (int y = 1; y < m_res_y - 1; ++y)
		for (int x = 1; x < m_res_x; ++x) {
			// TODO: update u
			if (p_solidObstacle->x()(x, y) > 0.5 || p_solidObstacle->x()(x - 1, y) > 0.5) {
				u(x, y) = 0;
			} else {
				u(x, y) = u(x, y) - m_dt / rho * (p(x, y) - p(x - 1, y)) / m_dx;
			}
		}

	// Same for velocity v_{i+1/2}.
	for (int y = 1; y < m_res_y; ++y)
		for (int x = 1; x < m_res_x - 1; ++x) {
			// TODO: update v
			if (p_solidObstacle->x()(x, y) > 0.5 || p_solidObstacle->x()(x, y - 1) > 0.5) {
				v(x, y) = 0;
			} else {
				v(x, y) = v(x, y) - m_dt / rho * (p(x, y) - p(x, y - 1)) / m_dx;
			}
		}
}

void FluidSim::advectValues() {
	// Densities live on the grid centers, the velocities on the MAC grid
	// Separate their computation to avoid confusion

	Array2d& d = p_density->x();
	Array2d& u = p_velocity->x();
	Array2d& v = p_velocity->y();

	Array2d& d_ = p_density_tmp->x();
	Array2d& u_ = p_velocity_tmp->x();
	Array2d& v_ = p_velocity_tmp->y();

	// Densities, grid centers
	for (int y = 1; y < m_res_y - 1; ++y) {
		for (int x = 1; x < m_res_x - 1; ++x) {
			if (p_solidObstacle->x()(x, y) > 0.5) {
				d_(x, y) = 0;
				continue;
			}

			double last_x_velocity = (u(x, y) + u(x + 1, y)) * 0.5;
			double last_y_velocity = (v(x, y) + v(x, y + 1)) * 0.5;

			double last_x = x - m_dt * last_x_velocity / m_dx;
			double last_y = y - m_dt * last_y_velocity / m_dx;

			// outside walls
			if (last_x <= 1) last_x = 1 + m_acc;
			if (last_y <= 1) last_y = 1 + m_acc;
			if (last_x >= m_res_x - 2) last_x = m_res_x - 2 - m_acc;
			if (last_y >= m_res_y - 2) last_y = m_res_y - 2 - m_acc;

			// obstacles
			int x_low = (int)(last_x + 0.5);
			int y_low = (int)(last_y + 0.5);
			if (p_solidObstacle->x()(x_low, x_low) > 0.5) {
				if (p_solidObstacle->x()(x_low - 1, y_low) < 0.5) last_x = x_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low + 1, y_low) < 0.5) last_x = x_low + 1 + m_acc;
				if (p_solidObstacle->x()(x_low, y_low - 1) < 0.5) last_y = y_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low, y_low + 1) < 0.5) last_y = y_low + 1 + m_acc;
			}

			x_low = (int)last_x;
			y_low = (int)last_y;
			int x_high = x_low + 1;
			int y_high = y_low + 1;

			double x_weight = last_x - x_low;
			double y_weight = last_y - y_low;

			d_(x, y) = (1 - y_weight) * ((1 - x_weight) * d(x_low, y_low) + x_weight * d(x_high, y_low)) + y_weight * ((1 - x_weight) * d(x_low, y_high) + x_weight * d(x_high, y_high));
		}
	}

	// Velocities (u), MAC grid
	for (int y = 1; y < m_res_y - 1; ++y) {
		for (int x = 1; x < m_res_x; ++x) {
			if (p_solidObstacle->x()(x, y) > 0.5 || p_solidObstacle->x()(x - 1, y) > 0.5) {
				u_(x, y) = 0;
				continue;
			}

			double last_x_velocity = u(x, y);
			double last_y_velocity = (v(x - 1, y) + v(x - 1, y + 1) + v(x, y) + v(x, y + 1)) / 4.0;

			double last_x = x - m_dt * last_x_velocity / m_dx;
			double last_y = y - m_dt * last_y_velocity / m_dx;

			// outside walls
			if (last_x <= 1) last_x = 1 + m_acc;
			if (last_y <= 1) last_y = 1 + m_acc;
			if (last_x >= m_res_x - 1) last_x = m_res_x - 1 - m_acc;
			if (last_y >= m_res_y - 2) last_y = m_res_y - 2 - m_acc;

			// obstacles
			int x_low = (int)last_x;
			int y_low = (int)(last_y + 0.5);
			if (p_solidObstacle->x()(x_low, x_low) > 0.5) {
				if (p_solidObstacle->x()(x_low - 1, y_low) < 0.5) last_x = x_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low + 1, y_low) < 0.5) last_x = x_low + 1 + m_acc;
				if (p_solidObstacle->x()(x_low, y_low - 1) < 0.5) last_y = y_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low, y_low + 1) < 0.5) last_y = y_low + 1 + m_acc;
			}

			x_low = (int)last_x;
			y_low = (int)last_y;
			int x_high = x_low + 1;
			int y_high = y_low + 1;

			double x_weight = last_x - x_low;
			double y_weight = last_y - y_low;

			u_(x, y) = (1 - y_weight) * ((1 - x_weight) * u(x_low, y_low) + x_weight * u(x_high, y_low)) + y_weight * ((1 - x_weight) * u(x_low, y_high) + x_weight * u(x_high, y_high));
		}
	}

	// Velocities (v), MAC grid
	for (int y = 1; y < m_res_y; ++y) {
		for (int x = 1; x < m_res_x - 1; ++x) {
			if (p_solidObstacle->x()(x, y) > 0.5 || p_solidObstacle->x()(x, y - 1) > 0.5) {
				v_(x, y) = 0;
				continue;
			}

			double last_x_velocity = (u(x, y - 1) + u(x, y) + u(x + 1, y - 1) + u(x + 1, y)) / 4.0;
			double last_y_velocity = v(x, y);

			double last_x = x - m_dt * last_x_velocity / m_dx;
			double last_y = y - m_dt * last_y_velocity / m_dx;

			// outside walls
			if (last_x <= 1) last_x = 1 + m_acc;
			if (last_y <= 1) last_y = 1 + m_acc;
			if (last_x >= m_res_x - 2) last_x = m_res_x - 2 - m_acc;
			if (last_y >= m_res_y - 1) last_y = m_res_y - 1 - m_acc;

			// obstacles
			int x_low = (int)(last_x + 0.5);
			int y_low = (int)last_y;
			if (p_solidObstacle->x()(x_low, x_low) > 0.5) {
				if (p_solidObstacle->x()(x_low - 1, y_low) < 0.5) last_x = x_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low + 1, y_low) < 0.5) last_x = x_low + 1 + m_acc;
				if (p_solidObstacle->x()(x_low, y_low - 1) < 0.5) last_y = y_low - 1 - m_acc;
				if (p_solidObstacle->x()(x_low, y_low + 1) < 0.5) last_y = y_low + 1 + m_acc;
			}

			x_low = (int)last_x;
			y_low = (int)last_y;
			int x_high = x_low + 1;
			int y_high = y_low + 1;

			double x_weight = last_x - x_low;
			double y_weight = last_y - y_low;

			v_(x, y) = (1 - y_weight) * ((1 - x_weight) * v(x_low, y_low) + x_weight * v(x_high, y_low)) + y_weight * ((1 - x_weight) * v(x_low, y_high) + x_weight * v(x_high, y_high));
		}
	}

	// Copy the values in temp to the original buffers
	for (int y = 1; y < m_res_y - 1; ++y)
		for (int x = 1; x < m_res_x - 1; ++x)
			d(x, y) = d_(x, y);
	for (int y = 1; y < m_res_y - 1; ++y)
		for (int x = 1; x < m_res_x; ++x)
			u(x, y) = u_(x, y);
	for (int y = 1; y < m_res_y; ++y)
		for (int x = 1; x < m_res_x - 1; ++x)
			v(x, y) = v_(x, y);
}

