#include "SpringSim.h"

bool SpringSim::advance() {
    // perform time integration with different integrators

	// use p_cube, m_spring, m_dt, m_gravity

	Eigen::Vector3d v = p_cube->getLinearVelocity();
	Eigen::Vector3d p = p_cube->getPosition();

    // TODO
    /*
    dx/dt = v
    dv/dt = -(k(|x|-L)+γv)/m+g
    */

	// note that it is required to update both m_sptring.end and p_cube's position
    switch (m_method) {
        case 0: // analytical solution
        {
            double alpha = -m_spring.damping / (2 * m_mass);
            double beta = std::sqrt(4 * m_mass * m_spring.stiffness - m_spring.damping * m_spring.damping) / (2 * m_mass);
            double A = -m_mass * -m_gravity.y() / m_spring.stiffness;
            double B = -alpha * A / beta;
            double x = -A + std::exp(alpha * m_time) * (A * std::cos(beta * m_time) + B * std::sin(beta * m_time));
            double v = std::exp(alpha * m_time) * ((alpha * A + beta * B) * std::cos(beta * m_time) + (alpha * B - beta * A) * std::sin(beta * m_time));
            p_cube->setPosition(Eigen::Vector3d(0, -x + m_spring.start.y() - m_spring.length, 0));
            p_cube->setLinearVelocity(Eigen::Vector3d(0, -v, 0));
        }
            break;
        case 1: // explicit euler
        {
            Eigen::Vector3d a = -(m_spring.stiffness * ((p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * v) / m_mass + m_gravity;
            p_cube->setPosition(p + m_dt * v);
            p_cube->setLinearVelocity(v + m_dt * a);
        }
            break;

        case 2: // symplectic euler
        {
            Eigen::Vector3d a = -(m_spring.stiffness * ((p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * v) / m_mass + m_gravity;
            p_cube->setLinearVelocity(v + m_dt * a);
            p_cube->setPosition(p + m_dt * p_cube->getLinearVelocity());
        }
            break;

        case 3: // midpoint
        {
            Eigen::Vector3d a = -(m_spring.stiffness * ((p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * v) / m_mass + m_gravity;
            Eigen::Vector3d v_half = v + 0.5 * m_dt * a;
            p_cube->setPosition(p + m_dt * v_half);
            Eigen::Vector3d p_half = p + 0.5 * m_dt * v;
            Eigen::Vector3d a_half = -(m_spring.stiffness * ((p_half - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * v) / m_mass + m_gravity;
            p_cube->setLinearVelocity(v + m_dt * a_half);
        }
            break;
        
        case 4: // RK4
        {
            Eigen::Vector3d k1_v, k2_v, k3_v, k4_v, k1_p, k2_p, k3_p, k4_p;
            k1_v = -(m_spring.stiffness * ((p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * v) / m_mass + m_gravity;
            k1_p = v;
            k2_v = -(m_spring.stiffness * ((p + 0.5 * m_dt * k1_p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * (v + 0.5 * m_dt * k1_v)) / m_mass + m_gravity;
            k2_p = v + 0.5 * m_dt * k1_v;
            k3_v = -(m_spring.stiffness * ((p + 0.5 * m_dt * k2_p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * (v + 0.5 * m_dt * k2_v)) / m_mass + m_gravity;
            k3_p = v + 0.5 * m_dt * k2_v;
            k4_v = -(m_spring.stiffness * ((p + m_dt * k3_p - m_spring.start).norm() - m_spring.length) * (m_spring.end - m_spring.start).normalized() + m_spring.damping * (v + m_dt * k3_v)) / m_mass + m_gravity;
            k4_p = v + m_dt * k3_v;
            p_cube->setLinearVelocity(v + m_dt * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6);
            p_cube->setPosition(p + m_dt * (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6);
        }
            break;

        default:
            std::cerr << m_method << " is not a valid integrator method."
                        << std::endl;
    }

	// update spring end position
	m_spring.end = p_cube->getPosition();


    // advance m_time
    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) {
        m_trajectories.back().push_back(p_cube->getPosition());
    }

    return false;
}