#include "DummySim.h"

using namespace std;
/*
 * Example simulation that changes the colors of a cube.
 */
DummySim::DummySim() : Simulation() {
	init();
}

static void addSphereMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, double radius, Eigen::Vector3d center = Eigen::Vector3d::Zero()) {
	int n = 20;
	int m = 20;
	int oldVSize = V.rows();
	int oldFSize = F.rows();
	V.conservativeResize(oldVSize + n * m, 3);
	F.conservativeResize(oldFSize + 2 * n * m, 3);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			double theta = 2 * M_PI * i / n;
			double phi = M_PI * j / (m - 1);
			V.row(oldVSize + i * m + j) << radius * sin(phi) * cos(theta) + center.x(), radius * -cos(phi) + center.y(), radius * sin(phi) * sin(theta) + center.z();
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m - 1; j++) {
			int i1 = (i + 1) % n;
			int j1 = j + 1;
			F.row(oldFSize + 2 * (i * m + j)) << oldVSize + i * m + j, oldVSize + i * m + j1, oldVSize + i1 * m + j;
			F.row(oldFSize + 2 * (i * m + j) + 1) << oldVSize + i1 * m + j, oldVSize + i * m + j1, oldVSize + i1 * m + j1;
		}
	}
}

void DummySim::init() {
	// create a sphere
	addSphereMesh(m_V, m_F, 5.0, Eigen::Vector3d(0, 5, 0));

	m_C.resize(m_V.rows(), 3);
	m_C.setZero();
	m_C.col(0).setOnes();

	m_dt = 0.01;

	reset();
}

void DummySim::resetMembers() {
	m_C.setZero();
	m_C.col(0).setOnes();
}

void DummySim::updateRenderGeometry() {
	m_renderV = m_V;
	m_renderF = m_F;
	m_renderC = m_C;
}

bool DummySim::advance() {
	// TODO: update m_V (and m_F maybe) here
	static double y = 5;
	static Eigen::MatrixXd m_V_start = m_V;
	for (int i = 0; i < m_V.rows(); i++) {
		m_V(i, 1) = (m_V_start(i, 1) - y) * (0.25 * cos(8 * m_time) + 0.75) + y - 10 * sin(2 * m_time);
		m_V(i, 0) = m_V_start(i, 0) * (0.25 * cos(5 * m_time) + 0.75);
		m_V(i, 2) = m_V_start(i, 2) * (0.25 * cos(5 * m_time) + 0.75); 
	}

	// advance step
	m_step++;
	m_time += m_dt;
	return false;
}

void DummySim::renderRenderGeometry(
	igl::opengl::glfw::Viewer &viewer) {
	viewer.data().set_mesh(m_renderV, m_renderF);
	viewer.data().set_colors(m_renderC);
}