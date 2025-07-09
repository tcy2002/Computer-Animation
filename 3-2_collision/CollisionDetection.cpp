#include "CollisionDetection.h"
#include <map>
#include <unordered_map>
#include <unordered_set>

void CollisionDetection::computeBroadPhase(int broadPhaseMethod) {
    // compute possible collisions
    m_overlappingBodys.clear();
    
	switch (broadPhaseMethod) {
	case 0: { // none
		for (size_t i = 0; i < m_objects.size(); i++) {
			for (size_t j = i + 1; j < m_objects.size(); j++) {
				m_overlappingBodys.push_back(std::make_pair(i, j));
			}
		}
		break;
	}
     
	case 1: {  // AABB
        // compute bounding boxes
        std::vector<AABB> aabbs(m_objects.size());
        for (size_t i = 0; i < aabbs.size(); i++) {
            aabbs[i].computeAABB(m_objects[i]);
        }
        for (size_t i = 0; i < m_objects.size(); i++) {
            for (size_t j = i + 1; j < m_objects.size(); j++) {
                // add pair of objects to possible collision if their
                // bounding boxes overlap
                if (aabbs[i].testCollision(aabbs[j])) {
                    m_overlappingBodys.push_back(std::make_pair(i, j));
                }
            }
        }
        break;
    }
	
	case 2: {
		// TODO: implement other broad phase algorithm
        // sweep and prune
        static int target_axis = 0;

        // compute bounding boxes
        std::map<RigidObject*, AABB> aabbs;
        for (size_t i = 0; i < m_sortedObjects.size(); i++) {
            aabbs[m_sortedObjects[i]].computeAABB(*m_sortedObjects[i]);
        }

        // sort on the target axis
        std::sort(m_sortedObjects.begin(), m_sortedObjects.end(), [&](const auto& a, const auto& b) {
            return aabbs[a].getMinCoord()[target_axis] < aabbs[b].getMinCoord()[target_axis];
        });

        Eigen::Vector3d s = Eigen::Vector3d::Zero(), s2 = Eigen::Vector3d::Zero();
        for (size_t i = 0; i < m_sortedObjects.size(); i++) {
            // update sum and sum of squares to calculate mean and variance
            AABB& aabb = aabbs[m_sortedObjects[i]];
            Eigen::Vector3d center = (aabb.getMinCoord() + aabb.getMaxCoord()) / 2;
            s += center;
            s2 += center.cwiseProduct(center);

            // test collision pairs
            for (size_t j = i + 1; j < m_sortedObjects.size(); j++) {
                auto& aabb2 = aabbs[m_sortedObjects[j]];
                if (aabb2.getMinCoord()[target_axis] > aabb.getMaxCoord()[target_axis]) break;
                if (aabb2.testCollision(aabb)) m_overlappingBodys.push_back(std::make_pair(j, i));
            }
        }
		break;

        // update axis sorted to be the one with the largest variance
        auto v = s2 - s.cwiseProduct(s) / (double)m_sortedObjects.size();
        target_axis = 0;
        if (v.y() > v.x()) target_axis = 1;
        if (v.z() > v[target_axis]) target_axis = 2;
	}
	}
}

bool isApproxVector3d(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return std::abs(a.x() - b.x()) < 1e-5 && std::abs(a.y() - b.y()) < 1e-5 && std::abs(a.z() - b.z()) < 1e-5;
}

namespace std {
template <>
struct hash<Eigen::Vector3d> {
    size_t operator()(const Eigen::Vector3d& v) const {
        return std::hash<int>()((int)(v.x() * 1e4)) ^ std::hash<int>()((int)(v.y() * 1e4)) ^ std::hash<int>()((int)(v.z() * 1e4));
    }
};
template <>
struct equal_to<Eigen::Vector3d> {
    bool operator()(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const {
        return isApproxVector3d(a, b);
    }
};
}

// test if an axis is a separating axis
static bool isSeperatingAxis(const Eigen::Vector3d& axis, const std::vector<Eigen::Vector3d>& Va, std::vector<std::vector<size_t>>& Fa, const std::vector<Eigen::Vector3d>& Vb, std::vector<std::vector<size_t>>& Fb, double& depth) {
    // project all vertices of both meshes onto the axis
    double min_a = std::numeric_limits<double>::infinity();
    double max_a = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < Va.size(); i++) {
        double proj = axis.dot(Va[i]);
        min_a = std::min(min_a, proj);
        max_a = std::max(max_a, proj);
    }

    double min_b = std::numeric_limits<double>::infinity();
    double max_b = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < Vb.size(); i++) {
        double proj = axis.dot(Vb[i]);
        min_b = std::min(min_b, proj);
        max_b = std::max(max_b, proj);
    }

    // check for overlap
    if (max_a < min_b || max_b < min_a) {
        return true; // separating axis found
    }

    // compute depth of intersection
    depth = std::min(max_a - min_b, max_b - min_a);
    return false; // no separating axis found
}

// find closest face to given axis (the normal of the face is almost the axis)
static std::vector<size_t> findClosestFace(const Eigen::Vector3d& sep_axis, const std::vector<Eigen::Vector3d>& V, std::vector<std::vector<size_t>>& F) {
    size_t closest_face;
    double min_dist_n = 2;
    for (size_t i = 0; i < F.size(); i++) {
        Eigen::Vector3d a = V[F[i][0]];
        Eigen::Vector3d b = V[F[i][1]];
        Eigen::Vector3d c = V[F[i][2]];
        Eigen::Vector3d n = (b - a).cross(c - a).normalized();
        double dist_n = std::abs(1 - sep_axis.dot(n));
        if (dist_n < min_dist_n) {
            min_dist_n = dist_n;
            closest_face = i;
        }
    }
    return F[closest_face];
}

static void clipFace(const std::vector<Eigen::Vector3d>& p_in, std::vector<Eigen::Vector3d>& p_out, const Eigen::Vector3d& normal, double eq) {
    size_t num_verts = p_in.size();
    if (num_verts < 2) return;

    Eigen::Vector3d first_vertex = p_in.back();
    double ds = normal.dot(first_vertex) + eq;

    for (int ve = 0; ve < num_verts; ve++) {
        Eigen::Vector3d second_vertex = p_in[ve];
        double de = normal.dot(second_vertex) + eq;
        if (ds < 0) {
            if (de < 0) {
                p_out.push_back(second_vertex);
            } else {
                double t = ds / (ds - de);
                p_out.push_back(first_vertex + t * (second_vertex - first_vertex));
            }
        } else {
            if (de < 0) {
                double t = ds / (ds - de);
                p_out.push_back(first_vertex + t * (second_vertex - first_vertex));
                p_out.push_back(second_vertex);
            }
        }
        first_vertex = second_vertex;
        ds = de;
    }
}

// intersect face against face, add contact points to result
static void intersectFaceFace(std::vector<Contact>& result, const Eigen::Vector3d& sep_axis, RigidObject* rba, const std::vector<Eigen::Vector3d>& Va, const std::vector<size_t>& fa, RigidObject* rbb, const std::vector<Eigen::Vector3d>& Vb, const std::vector<size_t>& fb) {
    Eigen::Vector3d v1 = Va[fa[0]], v2 = Va[fa[1]], v3 = Va[fa[2]];
    Eigen::Vector3d n_face = (v2 - v1).cross(v3 - v1).normalized();
    std::vector<Eigen::Vector3d> a_verts, b_verts;
    for (int i = 0; i < fb.size(); i++) {
        b_verts.push_back(Vb[fb[i]]);
    }

    std::vector<Eigen::Vector3d>* p_in = &b_verts;
    std::vector<Eigen::Vector3d>* p_out = &a_verts;
    for (size_t e0 = 0; e0 < fa.size(); e0++) {
        Eigen::Vector3d a = Va[fa[e0]];
        Eigen::Vector3d b = Va[fa[(e0 + 1) % fa.size()]];
        Eigen::Vector3d edge0 = a - b;
        Eigen::Vector3d plane_normal = -edge0.cross(n_face).normalized();
        double plane_eq = -a.dot(plane_normal);

        clipFace(*p_in, *p_out, plane_normal, plane_eq);
        std::swap(p_in, p_out);
        p_out->clear();
    }

    // only add points that are inside the face
    Eigen::Vector3d v_face = Va[fa[0]];
    for (size_t i = 0; i < p_in->size(); i++) {
        Eigen::Vector3d pt = (*p_in)[i];
        double depth = n_face.dot(pt - v_face);
        if (depth < 0) {
            Contact contact;
            contact.p = pt;
            contact.n = sep_axis;
            contact.a = rba;
            contact.b = rbb;
            result.push_back(contact);
        }
    }
}

// restruct the mesh to hull
void convertTriangleMesh2Hull(Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Eigen::Vector3d>& V_hull, std::vector<std::vector<size_t>>& F_hull) {
    std::unordered_map<Eigen::Vector3d, std::unordered_set<Eigen::Vector3d>> faces;
    for (size_t i = 0; i < F.rows(); i++) {
        Eigen::Vector3d v1 = V.row(F(i, 0)), v2 = V.row(F(i, 1)), v3 = V.row(F(i, 2));
        Eigen::Vector3d n = (v2 - v1).cross(v3 - v1).normalized();
        if (faces.find(n) == faces.end()) {
            faces[n] = {};
        }
        faces[n].insert(v1);
        faces[n].insert(v2);
        faces[n].insert(v3);
    }

    for (auto& face : faces) {
        std::vector<Eigen::Vector3d> verts(face.second.begin(), face.second.end());
        Eigen::Vector3d first = verts.front();
        std::sort(verts.begin(), verts.end(), [&](const auto& a, const auto& b) {
            if (isApproxVector3d(a, first)) return true;
            else if (isApproxVector3d(b, first)) return false;
            return (a - first).cross(b - first).dot(face.first) > 0;
        });

        F_hull.emplace_back();
        size_t first_idx = V_hull.size();
        for (size_t i = 0; i < verts.size(); i++) {
            V_hull.push_back(verts[i]);
            F_hull.back().push_back(first_idx + i);
        }
    }
}

void CollisionDetection::computeNarrowPhase(int narrowPhaseMethod) {
    switch (narrowPhaseMethod) {
    case 0: {
        // exhaustive
        // iterate through all pairs of possible collisions
        for (auto overlap : m_overlappingBodys) {
            std::vector<Contact> temp_contacts[2];
            // compute intersection of a with b first and intersectino
            // of b with a and store results in temp_contacts
            for (int switcher = 0; switcher < 2; switcher++) {
                RigidObject* a =
                    &m_objects[(!switcher) ? overlap.first
                                            : overlap.second];
                RigidObject* b =
                    &m_objects[(!switcher) ? overlap.second
                                            : overlap.first];

                Eigen::MatrixXd Va, Vb;
                Eigen::MatrixXi Fa, Fb;
                a->getMesh(Va, Fa);
                b->getMesh(Vb, Fb);

                // iterate through all faces of first object
                for (int face = 0; face < Fa.rows(); face++) {
                    // iterate through all edges of given face
                    for (size_t j = 0; j < 3; j++) {
                        int start = Fa(face, j);
                        int end = Fa(face, (j + 1) % 3);

                        // check if there is a collision
                        ContactType ct = isColliding(
                            Va.row(start), Va.row(end), Vb, Fb);

                        // find collision and check for duplicates
                        switch (ct) {
                            case ContactType::VERTEXFACE: {
                                auto ret = m_penetratingVertices.insert(
                                    Fa(face, j));
                                // if not already in set
                                if (ret.second) {
                                    Contact temp_col =
                                        findVertexFaceCollision(
                                            Va.row(Fa(face, j)), Vb,
                                            Fb);
                                    temp_col.a = a;
                                    temp_col.b = b;
                                    temp_col.type =
                                        ContactType::VERTEXFACE;
                                    temp_contacts[switcher].push_back(
                                        temp_col);
                                }
                                break;
                            }
                            case ContactType::EDGEEDGE: {
                                int orderedStart = std::min(start, end);
                                int orderedEnd = std::max(start, end);
                                auto ret = m_penetratingEdges.insert(
                                    std::make_pair(orderedStart,
                                                    orderedEnd));
                                // if not already in set
                                if (ret.second) {
                                    Contact temp_col =
                                        findEdgeEdgeCollision(
                                            Va.row(orderedStart),
                                            Va.row(orderedEnd), Vb, Fb);
                                    temp_col.a = a;
                                    temp_col.b = b;
                                    temp_col.type =
                                        ContactType::EDGEEDGE;
                                    temp_contacts[switcher].push_back(
                                        temp_col);
                                }
                                break;
                            }
                            case ContactType::NONE: {
                                break;
                            }
                        }
                    }
                }
                m_penetratingVertices.clear();
                m_penetratingEdges.clear();
            }

            // look for vertexFace
            bool found = false;
            for (int i = 0; i < 2; i++) {
                for (auto cont : temp_contacts[i]) {
                    if (cont.type == ContactType::VERTEXFACE) {
                        m_contacts.push_back(cont);
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (found) {
                continue;
            }

            // take single edgeedge if possible
            if (temp_contacts[0].size() > 0 &&
                temp_contacts[0].size() < temp_contacts[1].size()) {
                for (int i = 0; i < temp_contacts[0].size(); i++) {
                    m_contacts.push_back(temp_contacts[0][i]);
                }
            } else if (temp_contacts[1].size() > 0 &&
                        temp_contacts[0].size() >
                            temp_contacts[1].size()) {
                for (int i = 0; i < temp_contacts[1].size(); i++) {
                    m_contacts.push_back(temp_contacts[1][i]);
                }
            } else if (temp_contacts[0].size() > 0) {
                for (int i = 0; i < temp_contacts[0].size(); i++) {
                    m_contacts.push_back(temp_contacts[0][i]);
                }
            } else if (temp_contacts[1].size() > 0) {
                for (int i = 0; i < temp_contacts[1].size(); i++) {
                    m_contacts.push_back(temp_contacts[1][i]);
                }
            }
        }
        break;
    }

	case 1: {
		// TODO: implement other narrow phase algorithm
        // SAT
        for (auto overlap : m_overlappingBodys) {
            auto* a = m_sortedObjects[overlap.first];
            auto* b = m_sortedObjects[overlap.second];
            Eigen::MatrixXd Va, Vb;
            Eigen::MatrixXi Fa, Fb;
            a->getMesh(Va, Fa);
            b->getMesh(Vb, Fb);
            std::vector<Eigen::Vector3d> Va_hull, Vb_hull;
            std::vector<std::vector<size_t>> Fa_hull, Fb_hull;
            convertTriangleMesh2Hull(Va, Fa, Va_hull, Fa_hull);
            convertTriangleMesh2Hull(Vb, Fb, Vb_hull, Fb_hull);
            
            Eigen::Vector3d rel = a->getPosition() - b->getPosition();

            /// first: find potential separating axis
            double d_min = std::numeric_limits<double>::infinity();
            Eigen::Vector3d sep_axis; // always pointing from B to A
            bool flag = false;

            // test normal from convex mesh A
            for (size_t i = 0; i < Fa_hull.size(); i++) {
                Eigen::Vector3d v1 = Va_hull[Fa_hull[i][0]], v2 = Va_hull[Fa_hull[i][1]], v3 = Va_hull[Fa_hull[i][2]];
                auto normal = (v2 - v1).cross(v3 - v1).normalized();
                if (normal.dot(rel) > 0) continue;
                double depth;
                if (isSeperatingAxis(normal, Va_hull, Fa_hull, Vb_hull, Fb_hull, depth)) {
                    // separating axis found
                    flag = true;
                    break;
                }
                if (depth < d_min) {
                    d_min = depth;
                    sep_axis = -normal;
                }
            }
            if (flag) continue;

            // test normal from convex mesh B
            for (size_t i = 0; i < Fb_hull.size(); i++) {
                Eigen::Vector3d v1 = Vb_hull[Fb_hull[i][0]], v2 = Vb_hull[Fb_hull[i][1]], v3 = Vb_hull[Fb_hull[i][2]];
                auto normal = (v2 - v1).cross(v3 - v1).normalized();
                if (normal.dot(rel) < 0) continue;
                double depth;
                if (isSeperatingAxis(normal, Va_hull, Fa_hull, Vb_hull, Fb_hull, depth)) {
                    // separating axis found
                    flag = true;
                    break;
                }
                if (depth < d_min) {
                    d_min = depth;
                    sep_axis = normal;
                }
            }
            if (flag) continue;

            // test edge cross products
            std::unordered_set<Eigen::Vector3d> edges_a, edges_b;
            for (size_t i = 0; i < Fa_hull.size(); i++) {
                for (size_t j = 0; j < Fa_hull[i].size(); j++) {
                    auto& v1 = Va_hull[Fa_hull[i][j]];
                    auto& v2 = Va_hull[Fa_hull[i][(j + 1) % Fa_hull[i].size()]];
                    Eigen::Vector3d edge = (v1 - v2).normalized();
                    if (edges_a.find(edge) == edges_a.end() && edges_a.find(-edge) == edges_a.end()) {
                        edges_a.insert(edge);
                    }
                }
            }
            for (size_t i = 0; i < Fb_hull.size(); i++) {
                for (size_t j = 0; j < Fb_hull[i].size(); j++) {
                    auto& v1 = Vb_hull[Fb_hull[i][j]];
                    auto& v2 = Vb_hull[Fb_hull[i][(j + 1) % Fb_hull[i].size()]];
                    Eigen::Vector3d edge = (v1 - v2).normalized();
                    if (edges_b.find(edge) == edges_b.end() && edges_b.find(-edge) == edges_b.end()) {
                        edges_b.insert(edge);
                    }
                }
            }
            for (auto edge_a : edges_a) {
                for (auto edge_b : edges_b) {
                    auto normal = edge_a.cross(edge_b);
                    if (normal.norm() < 1e-5) continue;
                    normal.normalize();
                    if (normal.dot(rel) < 0) normal = -normal;
                    double depth;
                    if (isSeperatingAxis(normal, Va_hull, Fa_hull, Vb_hull, Fb_hull, depth)) {
                        // separating axis found
                        flag = true;
                        break;
                    }
                    if (depth < d_min) {
                        d_min = depth;
                        sep_axis = normal;
                    }
                }
            }
            if (flag) continue;

            /// second: find closest face on each mesh
            auto fa = findClosestFace(-sep_axis, Va_hull, Fa_hull);
            auto fb = findClosestFace(sep_axis, Vb_hull, Fb_hull);

            /// last: clip face against face and add contact points
            intersectFaceFace(m_contacts, sep_axis, a, Va_hull, fa, b, Vb_hull, fb);
        }
        
		break;
	}
    }
}

void CollisionDetection::applyImpulse(double eps) {
    // compute impulse for all contacts
    for (auto contact : m_contacts) {
        Eigen::Vector3d vrel_vec = contact.a->getVelocity(contact.p) -
                                    contact.b->getVelocity(contact.p);
        double vrel_n = contact.n.dot(vrel_vec);
        Eigen::Vector3d contact_t = (contact.n.cross(vrel_vec).cross(contact.n)).normalized();
        double vrel_t = contact_t.dot(vrel_vec);
        
		// TODO: compute impulse and update the following momentums
		Eigen::Vector3d r_a = contact.p - contact.a->getPosition();
        Eigen::Vector3d r_b = contact.p - contact.b->getPosition();

        double j_a_n = contact.n.dot((contact.a->getInertiaInv() * r_a.cross(contact.n)).cross(r_a));
        double j_b_n = contact.n.dot((contact.b->getInertiaInv() * r_b.cross(contact.n)).cross(r_b));
        double j_n = -(1 + eps) * vrel_n / (contact.a->getMassInv() + contact.b->getMassInv() + j_a_n + j_b_n);

        if (vrel_n < 0) {
            contact.a->setLinearMomentum(contact.a->getLinearMomentum() + j_n * contact.n);
            contact.b->setLinearMomentum(contact.b->getLinearMomentum() + j_n * -contact.n);
            contact.a->setAngularMomentum(contact.a->getAngularMomentum() + r_a.cross(j_n * contact.n));
            contact.b->setAngularMomentum(contact.b->getAngularMomentum() + r_b.cross(j_n * -contact.n));
        }
        if (vrel_vec.squaredNorm() - vrel_n * vrel_n > 1e-5) {
            double j_a_t = contact_t.dot((contact.a->getInertiaInv() * r_a.cross(contact_t)).cross(r_a));
            double j_b_t = contact_t.dot((contact.b->getInertiaInv() * r_b.cross(contact_t)).cross(r_b));
            double j_t = -(1 + eps) * vrel_t / (contact.a->getMassInv() + contact.b->getMassInv() + j_a_t + j_b_t);
            j_t = std::clamp(j_t, -m_friction * j_n, m_friction * j_n);

            contact.a->setLinearMomentum(contact.a->getLinearMomentum() + j_t * contact_t);
            contact.b->setLinearMomentum(contact.b->getLinearMomentum() + j_t * -contact_t);
            contact.a->setAngularMomentum(contact.a->getAngularMomentum() + r_a.cross(j_t * contact_t));
            contact.b->setAngularMomentum(contact.b->getAngularMomentum() + r_b.cross(j_t * -contact_t));
        }
    }
}
