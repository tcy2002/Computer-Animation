#include "CannonBallSim.h"


bool CannonBallSim::advance() {
    // perform time integration with different integrators

	// use p_ball, m_dt, m_gravity
    static Eigen::Vector3d v0 = p_ball->getLinearVelocity();
    static Eigen::Vector3d p0 = p_ball->getPosition();
	Eigen::Vector3d v = p_ball->getLinearVelocity();
	Eigen::Vector3d p = p_ball->getPosition();

    /*
    dx/dt = v
    dv/dt = g
    */

    // TODO
    switch (m_method) {
        case 0:
            // analytical solution
			// p(t) = v_0*t + 0.5*a*t^2
            p_ball->setPosition(p0 + v0 * m_time + 0.5 * m_gravity * m_time * m_time);
            p_ball->setLinearVelocity(v0 + m_time * m_gravity);
            break;

        case 1:
            // explicit euler
			// p' = p + dt*v
			// v' = v + dt*a
            p_ball->setPosition(p + m_dt * v);
            p_ball->setLinearVelocity(v + m_dt * m_gravity);
            break;

        case 2:
            // symplectic euler
			// v' = v + dt*a
			// p' = p + dt*v'
            p_ball->setLinearVelocity(v + m_dt * m_gravity);
            p_ball->setPosition(p + m_dt * p_ball->getLinearVelocity());
            break;

        case 3:
            // RK4
        {
            Eigen::Vector3d k1_v, k2_v, k3_v, k4_v, k1_p, k2_p, k3_p, k4_p;
            k1_v = m_gravity;
            k1_p = v;
            k2_v = m_gravity;
            k2_p = v + 0.5 * k1_v * m_dt;
            k3_v = m_gravity;
            k3_p = v + 0.5 * k2_v * m_dt;
            k4_v = m_gravity;
            k4_p = v + k3_v * m_dt;
            p_ball->setLinearVelocity(v + m_dt * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6);
            p_ball->setPosition(p + m_dt * (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6);

        }
            break;

        default:
            std::cerr << m_method << " is not a valid integrator method."
                        << std::endl;
    }

    // advance time
    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) {
        m_trajectories.back().push_back(p_ball->getPosition());
    }

    return false;
}