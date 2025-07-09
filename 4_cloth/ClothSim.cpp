#include "ClothSim.h"
#include "igl/Timer.h"

using Eigen::Vector2i;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Matrix<double, 3, 3> Matrix33;

// add matrix into triplet array with index [i*3, i*3+1, i*3+2]*[j*3, j*3+1, j*3+2]
static void fillInTriplet(int i, int j, std::vector<Triplet>& triplets, const Matrix33& matrix) {
    for(int a = 0; a < 3; ++a)
        for(int b = 0; b < 3; ++b)
            triplets.push_back(Triplet(i*3+a, j*3+b, matrix(a, b)));
}

// add v into vector b in index [i*3, i*3+1, i*3+2]
static void addInVector(int i, Eigen::VectorXd& b, const Eigen::Vector3d& v) {
    for(int a = 0; a < 3; ++a)
        b[i*3+a] += v[a];
}

// To test if a sparse matrix A is symmetric
static bool isSymmetrical(Eigen::SparseMatrix<double>& A) {
	bool ret = true;
	Eigen::SparseMatrix<double> res = Eigen::SparseMatrix<double>(A.transpose()) - A;
	for (unsigned int k = 0; k < res.outerSize(); k++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(res, k); it; ++it)
		{
			if (std::fabs(it.value()) > 1e-6) {
				std::cout<< "row:{} col:{} value: {}" << it.row() << it.col() << A.coeffRef(it.row(), it.col());
				return false;
			}
		}
	}
	return ret;
}

// useful functions:
// 1. particleCoordinateToIdx()
// 2. isValidCoor()
bool ClothSim::advance() {
    igl::Timer t;
    t.start();

    static auto explicit_euler = [&]() {
        std::vector<Eigen::Vector3d> f;
        f.resize(particles.size(), Eigen::Vector3d::Zero());

        static double L[3] = {dx, 2 * dx, std::sqrt(2) * dx};
        static double K[3] = {m_spring.k_struct, m_spring.k_bend, m_spring.k_shear};

        // TODO compute force
        for(int idx = 0; idx < particles.size(); ++idx) {
            Eigen::Vector2i a_coord = particleIdxToCoordinate(idx);

            // gravity
            f[idx] += m_gravity * m_mass;

            // structure/bending/shearing forces
            static std::vector<std::vector<std::pair<int, int>>> d = {
                {{-1, 0}, {1, 0}, {0, -1}, {0, 1}},
                {{-2, 0}, {2, 0}, {0, -2}, {0, 2}},
                {{-1, 1}, {-1, -1}, {1, -1}, {1, 1}}
            };
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                    Eigen::Vector2i b_coord(a_coord.x() + d[i][j].first, a_coord.y() + d[i][j].second);
                    int b_idx = particleCoordinateToIdx(b_coord);
                    if (isValidCoor(Eigen::Vector2i(b_coord))) {
                        Eigen::Vector3d p_ab = particles[idx].position - particles[b_idx].position;
                        Eigen::Vector3d p_ab_hat = p_ab.normalized();
                        double p_ab_norm = p_ab.norm();
                        Eigen::Vector3d v_ab = particles[idx].velocity - particles[b_idx].velocity;
                        Eigen::Vector3d force = K[i] * (L[i] - p_ab_norm) * p_ab_hat;
                        Eigen::Vector3d force_d = -m_spring.damping * v_ab.dot(p_ab_hat) * p_ab_hat;
                        f[idx] += force + force_d;
                    }
                }
            }

            // friction forces
            Eigen::Vector3d v = particles[idx].velocity;
            if (particles[idx].position.y() < 0 && (std::abs(v.x()) > 1e-6 || std::abs(v.z()) > 1e-6)) {
                Eigen::Vector3d v = particles[idx].velocity;
                Eigen::Vector3d fric = friction * ((1 + restitution) * v.y() / m_dt + std::abs(m_gravity.y())) * m_mass * -Eigen::Vector3d(v.x(), 0, v.z()).normalized();
                f[idx] += fric;
            }
        }

        // TODO constrain fixed particles
        for(int i = 0; i < fixedParticleIdx.size(); ++i) {
            f[fixedParticleIdx[i]] = Eigen::Vector3d::Zero();
        }

        // TODO update velocity and position of particles
        for(int idx = 0; idx < particles.size(); ++idx) {
            if (isFixedParticle[idx]) continue;
            Eigen::Vector3d a = f[idx] / m_mass;
            Eigen::Vector3d old_p = particles[idx].position;
            Eigen::Vector3d old_v = particles[idx].velocity;
            particles[idx].velocity = old_v + a * m_dt;
            particles[idx].position = old_p + old_v * m_dt;
            if (particles[idx].position.y() < 0) {
                particles[idx].position.y() = -1e-6;
                particles[idx].velocity.y() = -particles[idx].velocity.y() * restitution;
            }
        }
    };

    static auto implicit_euler = [&]() {
        std::vector<Triplet> triplets;
        Eigen::VectorXd b(particles.size()*3); b.setZero();
        SparseMatrix A(particles.size()*3, particles.size()*3);

        static double L[3] = {dx, 2 * dx, std::sqrt(2) * dx};
        static double K[3] = {m_spring.k_struct, m_spring.k_bend, m_spring.k_shear};

        // TODO compute triplets arrays and right-hand-size b
        // NOTE: remember to process fixed particles, use isFixedParticle to test if a particle is fixed
        for(int idx = 0; idx < particles.size(); ++idx) {
            Eigen::Vector2i a_coord = particleIdxToCoordinate(idx);

            Eigen::Matrix3d m = Eigen::Matrix3d::Identity() * m_mass;
            if (isFixedParticle[idx]) {
                m = Eigen::Matrix3d::Identity() * 1e100;
            }
            Eigen::Matrix3d dt_df_dv = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d dt2_df_dp = Eigen::Matrix3d::Zero();

            Eigen::Vector3d f_n = Eigen::Vector3d::Zero();
            Eigen::Vector3d dt_df_dp_v = Eigen::Vector3d::Zero();

            // gravity
            if (!isFixedParticle[idx]) {
                f_n += m_gravity * m_mass;
            }

            // structure/bending/shearing forces
            static std::vector<std::vector<std::pair<int, int>>> d = {
                {{-1, 0}, {1, 0}, {0, -1}, {0, 1}},
                {{-2, 0}, {2, 0}, {0, -2}, {0, 2}},
                {{-1, 1}, {-1, -1}, {1, -1}, {1, 1}}
            };
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                    Eigen::Vector2i b_coord(a_coord.x() + d[i][j].first, a_coord.y() + d[i][j].second);
                    int b_idx = particleCoordinateToIdx(b_coord);
                    if (isValidCoor(Eigen::Vector2i(b_coord))) {
                        Eigen::Vector3d p_ab = particles[idx].position - particles[b_idx].position;
                        Eigen::Vector3d p_ab_hat = p_ab.normalized();
                        double p_ab_norm = p_ab.norm();
                        Eigen::Vector3d v_ab = particles[idx].velocity - particles[b_idx].velocity;
                        Eigen::Vector3d force = K[i] * (L[i] - p_ab_norm) * p_ab_hat;
                        Eigen::Vector3d force_d = -m_spring.damping * v_ab.dot(p_ab_hat) * p_ab_hat;

                        Eigen::Matrix3d p_ab_pa = (Eigen::Matrix3d::Identity() - p_ab_hat * p_ab_hat.transpose()) / p_ab_norm;
                        Eigen::Matrix3d df_dp = K[i] * (-p_ab_hat * p_ab_hat.transpose() + (L[i] - p_ab_norm) * p_ab_pa);
                        df_dp += m_spring.damping * (-v_ab.dot(p_ab_hat) * p_ab_pa + p_ab_pa * -v_ab * p_ab_hat.transpose());

                        Eigen::Matrix3d dt_df_dva = m_dt * -m_spring.damping * p_ab_hat * p_ab_hat.transpose();
                        Eigen::Matrix3d dt2_df_dpa = m_dt * m_dt * df_dp;
                        dt_df_dv += dt_df_dva;
                        dt2_df_dp += dt2_df_dpa;

                        f_n += force + force_d;
                        dt_df_dp_v += m_dt * (df_dp * particles[idx].velocity - df_dp * particles[b_idx].velocity);

                        // add non-diagonal triplets
                        fillInTriplet(idx, b_idx, triplets, dt_df_dva + dt2_df_dpa);
                    }
                }
            }

            // friction forces
            Eigen::Vector3d v = particles[idx].velocity;
            if (particles[idx].position.y() < 0 && (std::abs(v.x()) > 1e-6 || std::abs(v.z()) > 1e-6)) {
                double fric_n = friction * ((1 + restitution) * v.y() / m_dt + std::abs(m_gravity.y())) * m_mass;
                Eigen::Vector3d fric_d = -Eigen::Vector3d(v.x(), 0, v.z()).normalized();
                f_n += fric_n * fric_d;
                double x2 = v.x() * v.x(), z2 = v.z() * v.z(), xz = v.x() * v.z();
                double d2 = x2 + z2;
                double d = std::sqrt(d2);
                Eigen::Matrix3d deri;
                deri << (d - x2 / d) / d2, 0, -xz / (d * d2),
                        0, 0, 0,
                        -xz / (d * d2), 0, (d - z2 / d) / d2;
                dt_df_dv += m_dt * (Eigen::Vector3d(0, -friction * (1 + restitution) / m_dt * m_mass, 0) * fric_d.transpose() + fric_n * -deri);
            }

            // add diagonal triplets
            fillInTriplet(idx, idx, triplets, m - dt_df_dv - dt2_df_dp);

            // add b
            addInVector(idx, b, m_dt * (f_n + dt_df_dp_v));
        }
        
        A.setFromTriplets(triplets.begin(), triplets.end());
        // NOTE: this is for debug purpose, make sure your A matrix is symmetric
        //if(!isSymmetrical(A)) {std::cout<<"wrong"<<std::endl; return true;}
        cgSolver.compute(A);
        Eigen::VectorXd x = cgSolver.solve(b);
        // std::cout<< "error:" << cgSolver.error() << ", iter:" << cgSolver.iterations() << std::endl;

        // TODO with x computed, update position and velocity of particles
        for(int i = 0; i < particles.size(); ++i) {
            if (isFixedParticle[i]) continue;
            Eigen::Vector3d old_p = particles[i].position;
            Eigen::Vector3d old_v = particles[i].velocity;
            particles[i].velocity = old_v + x.segment<3>(i*3);
            particles[i].position = old_p + particles[i].velocity * m_dt;
            if (particles[i].position.y() < 0) {
                particles[i].position.y() = -1e-6;
                particles[i].velocity.y() = -particles[i].velocity.y() * restitution;
            }
        }
    };

    static auto semi_implicit_euler = [&]() {
        // TODO
        std::vector<Eigen::Vector3d> f;
        f.resize(particles.size(), Eigen::Vector3d::Zero());

        static double L[3] = {dx, 2 * dx, std::sqrt(2) * dx};
        static double K[3] = {m_spring.k_struct, m_spring.k_bend, m_spring.k_shear};

        // TODO compute force
        for(int idx = 0; idx < particles.size(); ++idx) {
            Eigen::Vector2i a_coord = particleIdxToCoordinate(idx);

            // gravity
            f[idx] += m_gravity * m_mass;

            // structure/bending/shearing forces
            static std::vector<std::vector<std::pair<int, int>>> d = {
                {{-1, 0}, {1, 0}, {0, -1}, {0, 1}},
                {{-2, 0}, {2, 0}, {0, -2}, {0, 2}},
                {{-1, 1}, {-1, -1}, {1, -1}, {1, 1}}
            };
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 4; j++) {
                    Eigen::Vector2i b_coord(a_coord.x() + d[i][j].first, a_coord.y() + d[i][j].second);
                    int b_idx = particleCoordinateToIdx(b_coord);
                    if (isValidCoor(Eigen::Vector2i(b_coord))) {
                        Eigen::Vector3d p_ab = particles[idx].position - particles[b_idx].position;
                        Eigen::Vector3d p_ab_hat = p_ab.normalized();
                        double p_ab_norm = p_ab.norm();
                        Eigen::Vector3d v_ab = particles[idx].velocity - particles[b_idx].velocity;
                        Eigen::Vector3d force = K[i] * (L[i] - p_ab_norm) * p_ab_hat;
                        Eigen::Vector3d force_d = -m_spring.damping * v_ab.dot(p_ab_hat) * p_ab_hat;
                        f[idx] += force + force_d;
                    }
                }
            }

            // friction forces
            Eigen::Vector3d v = particles[idx].velocity;
            if (particles[idx].position.y() < 0 && (std::abs(v.x()) > 1e-6 || std::abs(v.z()) > 1e-6)) {
                Eigen::Vector3d fric = friction * ((1 + restitution) * v.y() / m_dt + std::abs(m_gravity.y())) * m_mass * -Eigen::Vector3d(v.x(), 0, v.z()).normalized();
                f[idx] += fric;
            }
        }

        // TODO constrain fixed particles
        for(int i = 0; i < fixedParticleIdx.size(); ++i) {
            f[fixedParticleIdx[i]] = Eigen::Vector3d::Zero();
        }

        // TODO update velocity and position of particles
        for(int idx = 0; idx < particles.size(); ++idx) {
            if (isFixedParticle[idx]) continue;
            Eigen::Vector3d a = f[idx] / m_mass;
            Eigen::Vector3d old_p = particles[idx].position;
            Eigen::Vector3d old_v = particles[idx].velocity;
            particles[idx].velocity = old_v + a * m_dt;
            particles[idx].position = old_p + particles[idx].velocity * m_dt;
            if (particles[idx].position.y() < 0) {
                particles[idx].position.y() = -1e-6;
                particles[idx].velocity.y() = -particles[idx].velocity.y() * restitution;
            }
        }
    };

    if(m_method == 0) { // euler
        explicit_euler();
    } else if (m_method == 1) { // implicit euler
        implicit_euler();
    } else if (m_method == 2) { // semi-implicit euler
        semi_implicit_euler();
    } else { // combination of a half step of explicit Euler and a half step of implicit Euler
        static bool i = true;
        if (i) {
            explicit_euler();
        } else {
            implicit_euler();
        }
        i = !i;
    }
    
    // advance m_time
    m_time += m_dt;
    m_step++;

    // log
    if ((m_step % m_log_frequency) == 0) { 
        std::cout<<t.getElapsedTimeInMicroSec()<<std::endl;
    }
    
    return false;
}