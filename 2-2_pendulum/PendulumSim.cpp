#include "PendulumSim.h"

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 2, 6> Matrix26d;
typedef Eigen::Matrix<double, 6, 6> Matrix66d;
typedef Eigen::Matrix<double, 2, 2> Matrix22d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef Eigen::Matrix<double, 3, 1> Vector3d;

bool PendulumSim::advance() {
    auto pos1 = p1->getPosition();
    auto pos2 = p2->getPosition();
    auto vel1 = p1->getLinearVelocity();
    auto vel2 = p2->getLinearVelocity();

    static double l1_2 = pos1.squaredNorm();
    static double l2_2 = (pos1 - pos2).squaredNorm();

    // TODO update positions and velocities of particle p1, p2
    // c1 = 0.5 * (x1.dot(x1) - l1^2) = 0
    // c2 = 0.5 * ((xi-x2).dot(xi-x2) - l2^2) = 0
    Matrix26d J;
    J << pos1(0), pos1(1), pos1(2), 0, 0, 0,
    (pos1(0)-pos2(0)), (pos1(1)-pos2(1)), (pos1(2)-pos2(2)), -(pos1(0)-pos2(0)), -(pos1(1)-pos2(1)), -(pos1(2)-pos2(2));
    Matrix26d J1;
    J1 << vel1(0), vel1(1), vel1(2), 0, 0, 0,
    (vel1(0)-vel2(0)), (vel1(1)-vel2(1)), (vel1(2)-vel2(2)), -(vel1(0)-vel2(0)), -(vel1(1)-vel2(1)), -(vel1(2)-vel2(2));

    Matrix66d W;
    double inv_mass1 = 1.0 / p1->getMass(), inv_mass2 = 1.0 / p2->getMass();
    W.setZero();
    W.diagonal() << inv_mass1, inv_mass1, inv_mass1, inv_mass2, inv_mass2, inv_mass2;
    Vector6d q1;
    q1 << vel1(0), vel1(1), vel1(2), vel2(0), vel2(1), vel2(2);
    Vector6d f;
    auto g1 = p1->getMass() * m_gravity, g2 = p2->getMass() * m_gravity;
    f << g1(0), g1(1), g1(2), g2(0), g2(1), g2(2);

    // feedback
    double k_s = 0.1, k_d = 0.1;
    double c1 = 0.5 * (pos1.dot(pos1) - l1_2);
    double c2 = 0.5 * ((pos1 - pos2).dot(pos1 - pos2) - l2_2);
    double c1_1 = pos1.dot(vel1);
    double c2_1 = (pos1 - pos2).dot(vel1 - vel2);

    Matrix22d A = J * W * J.transpose();
    Vector2d b = -J1 * q1 - J * W * f - k_s * Vector2d(c1, c2) - k_d * Vector2d(c1_1, c2_1);

    Vector2d lambda = A.ldlt().solve(b); //use ldlt to solve the linear system
    Vector6d F = J.transpose() * lambda + f;
    Vector3d a1 = F.head(3) * inv_mass1, a2 = F.tail(3) * inv_mass2;
    p1->setLinearVelocity(vel1 + m_dt * a1);
    p2->setLinearVelocity(vel2 + m_dt * a2);
    p1->setPosition(pos1 + m_dt * p1->getLinearVelocity());
    p2->setPosition(pos2 + m_dt * p2->getLinearVelocity());

    // advance m_time
    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) {
        m_trajectories[0].push_back(p1->getPosition());
        m_trajectories[1].push_back(p2->getPosition());
    }

    return false;
}