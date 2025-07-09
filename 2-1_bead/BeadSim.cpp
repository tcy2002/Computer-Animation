#include "BeadSim.h"


bool BeadSim::advance() {


    auto v = p_bead->getLinearVelocity();
    auto p = p_bead->getPosition();

    // TODO update position and velcity of p_bead
    // constraint C(x) = 0.5*(p.dot(p) - m_radius^2) = 0;
    double lamdda = ((-m_mass * m_gravity).dot(p) - m_mass * v.dot(v)) / (p.dot(p));
    
    // feedback
    auto f_n = -m_k * (p - m_radius * p.normalized());
    
    auto a = m_gravity + (lamdda * p + f_n) / m_mass;
    p_bead->setLinearVelocity(v + m_dt * a);
    p_bead->setPosition(p + m_dt * p_bead->getLinearVelocity());

    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) {
        m_trajectories.back().push_back(p_bead->getPosition());
    }
    return false;

}