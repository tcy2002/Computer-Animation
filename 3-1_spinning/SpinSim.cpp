#include "SpinSim.h"

bool SpinSim::advance() {
	Eigen::Vector3d w = p_body->getAngularVelocity();
	
	// TODO
	// update orientation
	switch (m_method) {
	case 0: {
		// matrix-based angular velocity
		double ang_speed = w.norm();
		if (ang_speed > 1e-5) {
			auto axis = w.normalized();
			auto delta_rot = Eigen::AngleAxisd(ang_speed * m_dt, axis).toRotationMatrix();
			auto rot = delta_rot * p_body->getRotationMatrix();
			p_body->setRotation(rot);
		}
		break;
	}

	case 1: {
		// quaternion-based angular velocity
		auto q = p_body->getRotation();
		auto dq = Eigen::Quaterniond(0, w.x(), w.y(), w.z()) * q;
		auto new_q = Eigen::Quaterniond(
			q.w() + 0.5 * m_dt * dq.w(), 
			q.x() + 0.5 * m_dt * dq.x(), 
			q.y() + 0.5 * m_dt * dq.y(), 
			q.z() + 0.5 * m_dt * dq.z()
		);
		new_q.normalize();
		p_body->setRotation(new_q);
		break;
	}
	default:
            std::cerr << m_method << " is not a valid rotation method."
                        << std::endl;
	}

	// advance time
	m_time += m_dt;
	m_step++;

	return false;
}