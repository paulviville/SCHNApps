#include <schnapps/plugins/tubular_mesh/vessels_builder.h>

#include <cgogn/geometry/types/eigen.h>
#include <cgogn/geometry/functions/orientation.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/geometry/functions/basics.h>

#include <unordered_map>

namespace schnapps
{

namespace vessels_building{
    bool build_hexmesh(UGraph& ug, CMap3& cmap3, CMap2& cmap2){
        add_radius(ug, 1.0);
//        ug_stats(ug);

        bool valid;
        /// collecting pre existing graph attributes
        UG_Attributes ug_attribs;
        valid = ug_get_attributes(ug, ug_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to load attributes";
            return false;
        }
        cgogn_log_info("success: ") << "graph attributes collected";

        /// preprocessing, cleaning, and subdividing the graph
        valid = ug_preprocessing(ug, ug_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to preprocess";
            return false;
        }
        cgogn_log_info("success: ") << "graph preprocessed";

        /// collecting simplified graph information: branching points and branches
        Graph simplified_graph = ug_analysis(ug);
        if(!simplified_graph.valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to analyze";
            return false;
        }
        cgogn_log_info("success: ") << "graph analyzed";
        cgogn_log_info("") << simplified_graph.intersections.size() << " intersections; "
                                           << simplified_graph.branches.size() << " branches";

        CMap2_Attributes m2_attribs;
        ug_attribs.m2_CC = ug.add_attribute<Dart, UGVertex>("m2_CC");
        ug_attribs.m2_interface = ug.add_attribute<Dart, UGDart>("m2_interface");
        /// building connection interfaces
        valid = build_connection_interfaces(ug, ug_attribs, cmap2);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to build connection interfaces";
            return false;
        }
        cgogn_log_info("success: ") << "connection interfaces built";

//        cmap2.create_embedding<M2Vertex::ORBIT>();
        m2_attribs.position = cmap2.add_attribute<VEC3, M2Vertex>("position");
        m2_attribs.center = cmap2.add_attribute<VEC3, M2Volume>("center");
        m2_attribs.f_branch = cmap2.add_attribute<Dart, M2Face>("f_branches");
        m2_attribs.f_point = cmap2.add_attribute<VEC3, M2Face>("f_points");
        ug_attribs.frames = ug.add_attribute<Frame, UGDart>("frames");

        valid = complete_intersections(ug, ug_attribs, simplified_graph, cmap2, m2_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to complete connection interfaces";
            return false;
        }
        cgogn_log_info("success: ") << "connection interfaces completed";

        valid = propagate_frames(ug, ug_attribs, simplified_graph, cmap2);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to propagate frames";
            return false;
        }
        cgogn_log_info("success: ") << "frames propagated";


        valid = set_interfaces_geometry(ug, ug_attribs, cmap2, m2_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to build interface geometry";
            return false;
        }
        cgogn_log_info("success: ") << "interface geometry built";

        m2_attribs.edge_pos = cmap2.add_attribute<VEC3, M2Edge>("edge_position");
        valid = set_edge_geometry(ug, ug_attribs, cmap2, m2_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to set edge geometry";
            return false;
        }
        cgogn_log_info("success: ") << "interface edge geometry built";

        m2_attribs.connections = cmap2.add_attribute<Dart, M2Dart>("connections");
        valid = build_branch_sections(ug, ug_attribs, cmap2, m2_attribs, cmap3);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to build branch sections";
            return false;
        }
        cgogn_log_info("success: ") << "branch sections built";

        valid = sew_sections(cmap2, m2_attribs, cmap3);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to sew sections";
            return false;
        }
        cgogn_log_info("success: ") << "sections sown";

        valid = set_m3_geometry(cmap2, m2_attribs, cmap3);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to set m3 geometry";
            return false;
        }
        cgogn_log_info("success: ") << "m3 geometry set";


        /// removing attributes meant only for construction
        ug.remove_attribute(ug_attribs.m2_interface);
        ug.remove_attribute(ug_attribs.frames);

        cmap2.remove_attribute(m2_attribs.center);
        cmap2.remove_attribute(m2_attribs.f_branch);
        cmap2.remove_attribute(m2_attribs.f_point);
        cmap2.remove_attribute(m2_attribs.connections);
        return true;
    }

    void add_radius(UGraph& ug, Scalar radius){
        UGraph::VertexAttribute<Scalar> radii = ug.get_attribute<Scalar, UGVertex>("radius");
        if(!radii.is_valid()){
            radii = ug.add_attribute<Scalar, UGVertex>("radius");
            ug.foreach_cell([&](UGVertex ugv){
                uint32 valence = ug.nb_darts_of_orbit(ugv);
                Scalar rad = radius;
                if(valence > 2){
                    rad = radius * (1 + 0.15*(valence-2));
                }
                radii[ugv] = rad;
            });
        }
    }

    bool ug_stats(const UGraph &ug){
        bool valid = true;
        std::unordered_map<uint32, uint32> vertex_valence;
        uint32 nb_isolated = 0;
        uint32 nb_intersections = 0;
        ug.foreach_cell([&](UGVertex ugv){
            uint32 valence = ug.nb_darts_of_orbit(ugv);
            Dart d0, d1;
            switch(valence){
                case 1:
                    if (ug.is_boundary(ug.alpha0(ugv.dart))){
                        nb_isolated++;
                        valid = false;
                    }
                    break;
                case 2:
                    break;
                default:
                    vertex_valence[valence]++;
                    nb_intersections++;
                    break;
            }
        });
        cgogn_log_info("Graph stats: ") << "valid = " << valid;
        cgogn_log_info("nb vertices: ") << ug.nb_cells<UGVertex>();
        cgogn_log_info("nb isolated vertices: ") << nb_isolated;
        cgogn_log_info("nb intersections: ") << nb_intersections;
        for(const auto& n : vertex_valence) cgogn_log_info("valence: ") << n.first << " \tnb: " << n.second;
        UGraph::VertexAttribute<Scalar> radii = ug.get_attribute<Scalar, UGVertex>("radius");
        if (!radii.is_valid())
            valid = false;
        else{
            ug.foreach_cell([&](UGVertex ugv){
                cgogn_log_info("vertex: ") << ug.embedding(ugv.dart, UGVertex::ORBIT) << "  rad: " << radii[ugv];
            });
        }
        return valid;
    }

    bool ug_get_attributes(UGraph& ug, UG_Attributes& ug_attribs){
        bool valid = true;
        ug_attribs.radii = ug.get_attribute<Scalar, UGVertex>("radius");
        ug_attribs.position = ug.get_attribute<VEC3, UGVertex>("position");
        if(!ug_attribs.radii.is_valid()){
            cgogn_log_info("ug_get_attribute error: ") << "radii attribute invalid";
            valid &= false;
        }
        if(!ug_attribs.position.is_valid()){
            cgogn_log_info("ug_get_attribute error: ") << "position attribute invalid";
            valid &= false;
        }

        return valid;
    }

    bool ug_preprocessing(UGraph &ug, UG_Attributes& ug_attribs){
        bool valid;
        valid = ug_clean(ug, ug_attribs);
        if(!valid){
            cgogn_log_info("ug_preprocessing error: ") << "clean up failed";
            return false;
        }
        valid = ug_sudivide(ug, ug_attribs);
        if(!valid){
            cgogn_log_info("ug_preprocessing error: ") << "subdivision failed";
            return false;
        }

        return true;
    }

    bool ug_clean(UGraph &ug, UG_Attributes& ug_attribs){
        bool valid;
        /// checks if the graph is empty
        valid = (ug.nb_cells<UGVertex>() > 0);
        if(!valid){
            cgogn_log_info("ug_clean error: ") << "0 vertices in graph";
            return false;
        }

        ug.foreach_cell([&](UGVertex ugv){
            uint32 valence = ug.nb_darts_of_orbit(ugv);
            /// if a vertex is isolated it is destroyed
            if(valence == 1 && ug.is_boundary(ug.alpha0(ugv.dart)))
                    ug.remove_vertex(ugv);
        });

        /// if the graph contains no edges it's invalid
        valid = (ug.nb_cells<UGEdge>() > 0);
        if(!valid){
            cgogn_log_info("ug_clean error: ") << "0 edges in graph";
            return false;
        }

        // add radius+position check ?
        return true;
    }

    bool ug_sudivide(UGraph& ug, UG_Attributes& ug_attribs){
        UGDartMarker ug_dart_marker(ug);
        UGraph::VertexAttribute<VEC3>& ug_pos = ug_attribs.position;
        UGraph::VertexAttribute<Scalar>& ug_radii = ug_attribs.radii;

        UGCache edge_cache(ug);
        edge_cache.build<UGEdge>();
        while(edge_cache.size<UGEdge>()){
            ug.foreach_cell([&](UGEdge uge){
                VEC3 p1 = ug_attribs.position[uge.dart];
                VEC3 p2 = ug_attribs.position[ug.alpha0(uge.dart)];
                Scalar radius = (ug_radii[uge.dart] + ug_radii[ug.alpha0(uge.dart)])/2;
                if((p1 - p2).norm() > 2 * radius){
                    UGVertex ugv = ug.cut_edge(uge);
                    ug_attribs.position[ugv] = (p1 + p2) / 2;
                    ug_radii[ugv] = radius;
                    ug_dart_marker.mark_orbit(UGEdge(ugv.dart));
                    ug_dart_marker.mark_orbit(UGEdge(ug.alpha1(ugv.dart)));
                }
            }, edge_cache);

            edge_cache.clear<UGEdge>();
            edge_cache.build<UGEdge>([&](UGEdge uge){
                return ug_dart_marker.is_marked(uge.dart);
            });
            ug_dart_marker.unmark_all();
        }
        return true;
    }

    Graph ug_analysis(const UGraph& ug){
        Graph graph;
        graph.valid = true;

        /// looking for all branches extremities
        std::vector<UGVertex> extremities;
        ug.foreach_cell([&](UGVertex ugv) -> bool{
            switch(ug.nb_darts_of_orbit(ugv)){
                case 1:
                    extremities.push_back(ugv);
                    break;
                case 2:
                    break;
                default:
                    extremities.push_back(ugv);
                    graph.intersections.push_back(ugv);
                    break;
            }
            return true;
        });

        if(!graph.valid){
            cgogn_log_info("error ug_analysis: ") << "error in extremities gathering";
            return graph;
        }

        /// building branches
        UGDartMarker ug_marker(ug);
        for(UGVertex ugv : extremities){
            ug.foreach_dart_of_orbit(ugv, [&](Dart d0){
                if(!ug_marker.is_marked(d0)){
                    ug_marker.mark(d0);
                    Dart d1 = ug.alpha0(d0);
                    ug_marker.mark(d1);
                    while(ug.nb_darts_of_orbit(UGVertex(d1)) == 2){
                        d1 = ug.alpha1(d1);
                        ug_marker.mark(d1);
                        d1 = ug.alpha0(d1);
                        ug_marker.mark(d1);
                    }
                    graph.branches.push_back(Branch(d0, d1));
                }
            });
        }

        /// all encompassing simplified graph
        ug.foreach_dart([&](Dart d) -> bool {
            if(!ug_marker.is_marked(d)){
                graph.valid = false;
                cgogn_log_info("error ug_analysis: ") << "graph not fully covered";
                return false;
            }
            return true;
        });
        return graph;
    }

    VEC3 project_on_sphere(VEC3 P, Scalar R, VEC3 C){
        return C + (P - C).normalized() * R;
    }

    bool build_connection_interfaces(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2){
        bool valid;
        ug.foreach_cell([&](UGVertex ugv) -> bool {
            switch (ug.nb_darts_of_orbit(ugv)) {
            case 1:
                valid = build_interface_1(ug, ug_attribs, cmap2, ugv);
                return valid;
            case 2:
                valid = build_interface_2(ug, ug_attribs, cmap2, ugv);
                return valid;
            default:
                valid = build_interface_n(ug, ug_attribs, cmap2, ugv);
                return valid;
            }
        });

        return valid;
    }

    bool build_interface_1(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv){
        Dart f = cmap2.add_face(4).dart;
        ug_attribs.m2_CC[ugv.dart] = f;
        ug_attribs.m2_interface[ugv.dart] = f;
        return true;
    }

    bool build_interface_2(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv){
        M2Builder m2builder(cmap2);
        std::vector<Dart> Fd = {m2builder.add_face_topo_fp(4), m2builder.add_face_topo_fp(4)};
        m2builder.phi2_sew(Fd[0], Fd[1]);
        m2builder.phi2_sew(cmap2.phi<1>(Fd[0]), cmap2.phi<111>(Fd[1]));
        m2builder.phi2_sew(cmap2.phi<11>(Fd[0]), cmap2.phi<11>(Fd[1]));
        m2builder.phi2_sew(cmap2.phi<111>(Fd[0]), cmap2.phi<1>(Fd[1]));
        ug_attribs.m2_CC[ugv] = Fd[0];
        Dart d = ugv.dart;
        ug_attribs.m2_interface[d] = Fd[0];
        ug_attribs.m2_interface[ug.alpha1(d)] = cmap2.phi<12>(Fd[0]);
        return true;
    }

    bool build_interface_n(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv){
        M2Builder m2builder(cmap2);
        std::vector<Dart> Fd = {m2builder.add_face_topo_fp(4), m2builder.add_face_topo_fp(4), m2builder.add_face_topo_fp(4)};

        m2builder.phi2_sew(Fd[0], cmap2.phi<1>(Fd[1]));
        m2builder.phi2_sew(cmap2.phi<1>(Fd[0]), Fd[2]);
        m2builder.phi2_sew(cmap2.phi<11>(Fd[0]), cmap2.phi<111>(Fd[2]));
        m2builder.phi2_sew(cmap2.phi<111>(Fd[0]), cmap2.phi<11>(Fd[1]));
        m2builder.phi2_sew(cmap2.phi<1>(Fd[2]), Fd[1]);
        m2builder.phi2_sew(cmap2.phi<11>(Fd[2]), cmap2.phi<111>(Fd[1]));

        ug_attribs.m2_CC[ugv] = Fd[0];

        return true;
    }

    bool complete_intersections(UGraph& ug, UG_Attributes& ug_attribs, const Graph& graph, CMap2& cmap2, CMap2_Attributes& m2_attribs){
        for(UGVertex ugv : graph.intersections){
            if(!complete_intersection(ug, ug_attribs, cmap2, m2_attribs, ugv)){
                cgogn_log_info("error complete_intersections") << "failed to complete intersection";
                return false;
            }
        }

        for(UGVertex ugv : graph.intersections){
            if(!create_intersection_frames(ug, ug_attribs, cmap2, m2_attribs, ugv)){
                cgogn_log_info("error complete_intersections") << "failed to create intersection frames";
                return false;
            }
        }
        return true;
    }

    bool complete_intersection(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv){
        switch (ug.nb_darts_of_orbit(ugv)) {
        case 3:
            return complete_intersection_3(ug, ug_attribs, cmap2, m2_attribs, ugv);
        default:
            return complete_intersection_n(ug, ug_attribs, cmap2, m2_attribs, ugv);
        }
        return true;
    }

    bool complete_intersection_3(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv){
        VEC3 center = ug_attribs.position[ugv];
        Scalar radius = ug_attribs.radii[ugv];
        Dart cc = ug_attribs.m2_CC[ugv];
//        m2_attribs.center[cc] = center;

        std::vector<Dart> Fd = {cc, cmap2.phi<2111>(cc), cmap2.phi<12>(cc)};
        std::vector<Dart> Pdart = {ugv.dart, ug.alpha1(ugv.dart), ug.alpha1(ug.alpha1(ugv.dart))};
        std::vector<VEC3> Ppos = {
            project_on_sphere(ug_attribs.position[ug.alpha0(Pdart[0])], radius, center),
            project_on_sphere(ug_attribs.position[ug.alpha0(Pdart[1])], radius, center),
            project_on_sphere(ug_attribs.position[ug.alpha0(Pdart[2])], radius, center)
        };
        VEC3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
        std::vector<VEC3> Q = {center + V * radius, center - V * radius};
        std::vector<VEC3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
                 center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
                 center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

        M[0] = mean_dir(center, radius, M[0], {Ppos[0], Ppos[1]});
        M[1] = mean_dir(center, radius, M[1], {Ppos[2], Ppos[1]});
        M[2] = mean_dir(center, radius, M[2], {Ppos[0], Ppos[2]});

        m2_attribs.position[Fd[0]] = M[0];
        m2_attribs.position[Fd[1]] = M[1];
        m2_attribs.position[Fd[2]] = M[2];
        m2_attribs.position[cmap2.phi<1>(Fd[0])] = Q[0];
        m2_attribs.position[cmap2.phi<111>(Fd[0])] = Q[1];

        ug_attribs.m2_interface[Pdart[0]] = Fd[0];
        ug_attribs.m2_interface[Pdart[1]] = Fd[1];
        ug_attribs.m2_interface[Pdart[2]] = Fd[2];
        return true;
    }

    bool complete_intersection_n(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv){
        using cgogn::geometry::intersection_ray_triangle;
        VEC3 center = ug_attribs.position[ugv];
        Scalar radius = ug_attribs.radii[ugv];
        Dart cc = ug_attribs.m2_CC[ugv];
//        m2_attribs.center[cc] = center;

        std::vector<Dart> Fd = {cc, cmap2.phi<2111>(cc), cmap2.phi<12>(cc)};

        std::vector<Dart> Pdart;
        std::vector<VEC3> Ppos;
        ug.foreach_adjacent_vertex_through_edge(ugv, [&](UGVertex ugv){
            Ppos.push_back(project_on_sphere(ug_attribs.position[ugv], radius, center));
            Pdart.push_back(ug.alpha0(ugv.dart));
        });

        VEC3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
        std::vector<VEC3> Q = {center + V * radius, center - V * radius};
        std::vector<VEC3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
                 center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
                 center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

        M[0] = mean_dir(center, radius, M[0], {Ppos[0], Ppos[1]});
        M[1] = mean_dir(center, radius, M[1], {Ppos[2], Ppos[1]});
        M[2] = mean_dir(center, radius, M[2], {Ppos[0], Ppos[2]});

        m2_attribs.position[Fd[0]] = M[0];
        m2_attribs.position[Fd[1]] = M[1];
        m2_attribs.position[Fd[2]] = M[2];
        m2_attribs.position[cmap2.phi<1>(Fd[0])] = Q[0];
        m2_attribs.position[cmap2.phi<111>(Fd[0])] = Q[1];

        m2_attribs.f_point[Fd[0]] = Ppos[0];
        m2_attribs.f_point[Fd[1]] = Ppos[1];
        m2_attribs.f_point[Fd[2]] = Ppos[2];

        m2_attribs.f_branch[Fd[0]] = Pdart[0];
        m2_attribs.f_branch[Fd[1]] = Pdart[1];
        m2_attribs.f_branch[Fd[2]] = Pdart[2];

        for(uint i = 3; i < Ppos.size(); i++){
            VEC3 P0 = Ppos[i];
            Dart F0 = Pdart[i];
            M2Face face2cut;
            std::vector<VEC3> Quadp; //positions
            std::vector<Dart> Quadd; //darts

            bool face_found = false;
            cmap2.foreach_cell([&](M2Face f) -> bool {
                Dart fd = convex_quad(cmap2, m2_attribs, f.dart);
                Quadd = {fd, cmap2.phi<1>(fd),
                         cmap2.phi<11>(fd),
                         cmap2.phi<111>(fd)};
                Quadp = {m2_attribs.position[Quadd[0]],
                         m2_attribs.position[Quadd[1]],
                         m2_attribs.position[Quadd[2]],
                         m2_attribs.position[Quadd[3]]};
                if((face_found = (intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[1], Quadp[2])
                                 || intersection_ray_triangle(center, P0 - center, Quadp[0], Quadp[2], Quadp[3])))){
                    face2cut = f;
                    return false;
                }
                return true;
            });
            if(!face_found){
                cgogn_log_info("error complete_intersection_n: ") << "no face found";
                return false;
            }

            VEC3 P1 = m2_attribs.f_point[M2Face(face2cut)];
            Dart F1 = m2_attribs.f_branch[M2Face(face2cut)];

            Dart cut0, cut1;
            VEC3 AC = (Quadp[2] - Quadp[0]).normalized();
            VEC3 BD = (Quadp[3] - Quadp[1]).normalized();
            VEC3 P0P1 = (P1 - P0).normalized();
            if(abs(AC.dot(P0P1)) < abs(BD.dot(P0P1))){
                cut0 = Quadd[0]; cut1 = Quadd[2];
            } else {
                cut0 = Quadd[1]; cut1 = Quadd[3];
            }
            cmap2.cut_face(cut0, cut1);
            M2Vertex v = cmap2.cut_edge(M2Edge(cmap2.phi_1(cut0)));
            m2_attribs.position[v] = project_on_sphere((P0 + P1) * Scalar(0.5), radius, center);

            Dart newFace, oldFace;
            VEC3 out0 = m2_attribs.position[v] - center;
            VEC3 out1 = ((P0 - m2_attribs.position[cut0]).normalized().cross((P1 - m2_attribs.position[cut0]).normalized()));
            if(out1.dot(out0) >= 0){
                newFace = cut0; oldFace = cut1;
            }
            else{
                newFace = cut1; oldFace = cut0;
            }

            m2_attribs.f_branch[M2Face(newFace)] = F0;
            m2_attribs.f_branch[M2Face(oldFace)] = F1;
            m2_attribs.f_point[M2Face(newFace)] = P0;
            m2_attribs.f_point[M2Face(oldFace)] = P1;

            cmap2.foreach_incident_vertex(M2Face(newFace), [&](M2Vertex m2v){
                std::vector<VEC3> points;
                cmap2.foreach_incident_face(m2v, [&](M2Face m2f){
                    points.push_back(m2_attribs.f_point[m2f]);
                });
                m2_attribs.position[m2v] = mean_dir(center, radius, m2_attribs.position[m2v], points);
            });
        }

        cmap2.foreach_incident_face(M2Volume(cc), [&](M2Face m2f){
            ug_attribs.m2_interface[m2_attribs.f_branch[m2f]] = m2f.dart;
        });
        return true;
    }

    bool create_intersection_frames(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv){
        VEC3 center = ug_attribs.position[ugv];
//        Dart cc = ug_attribs.m2_CC[ugv];
        ug.foreach_dart_of_orbit(ugv, [&](Dart d){
            Dart d0 = ug_attribs.m2_interface[d];
            Dart d1 = cmap2.phi<11>(d0);
            VEC3 R, S, T, diag, temp;
            T = (ug_attribs.position[ug.alpha0(d)] - center).normalized();
            diag =(m2_attribs.position[d1] - m2_attribs.position[d0]).normalized();
//            S = T.cross(diag).normalized();
//            R = S.cross(T).normalized();
            R = diag.cross(T).normalized();
            S = T.cross(R).normalized();
            ug_attribs.frames[d].col(0) = R;
            ug_attribs.frames[d].col(1) = S;
            ug_attribs.frames[d].col(2) = T;
        });
        return true;
    }

    Dart convex_quad(CMap2& cmap2, CMap2_Attributes& m2_attribs, Dart f){
        Dart res;
        VEC3 A = m2_attribs.position[M2Vertex(f)];
        VEC3 B = m2_attribs.position[M2Vertex(cmap2.phi<1>(f))];
        VEC3 C = m2_attribs.position[M2Vertex(cmap2.phi<11>(f))];
        VEC3 D = m2_attribs.position[M2Vertex(cmap2.phi<111>(f))];
        VEC3 AB = B - A;
        VEC3 AC = C - A;
        VEC3 AD = D - A;
        VEC3 N0 = AC.cross(AD);
        VEC3 N1 = AB.cross(AC);
        VEC3 N = N0.cross(N1).normalized();

        if(N.dot(AC) < 0){
            res = cmap2.phi<1>(f);
        }
        else{
            res = f;
        }
        return res;
    }

    VEC3 mean_dir(VEC3 center, Scalar radius, VEC3 point, std::vector<VEC3> points){
        using Quat = Eigen::Quaterniond;
        uint32 valence = points.size();

        std::vector<VEC3> directions;
        for(VEC3 p : points){
            directions.push_back((p - center).normalized());
        }
        VEC3 avg_dir = (point - center).normalized();

        std::vector<Quat> rotations;
        rotations.reserve(valence);

//        for(uint i = 0; i < valence; i++){
            for(VEC3 dir : directions){
                Quat q = Quat::FromTwoVectors(avg_dir, dir);
                q.normalize();
                rotations.push_back(q);
            }

            Eigen::MatrixXd m(4, valence);
            for(uint32 j = 0; j < valence; ++j){
                const Quat& q = rotations[j];
                m.col(j) = VEC4(q.w(), q.x(), q.y(), q.z());
            }

            Eigen::MatrixXd mm = m * m.transpose();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mm);
            VEC4 r = es.eigenvectors().col(3);

            Quat mean_rot(r[0], r[1], r[2], r[3]);
            mean_rot.normalize();
            avg_dir = mean_rot._transformVector(avg_dir);
            rotations.clear();
//        }

        return avg_dir * radius + center;
    }

    bool propagate_frames(UGraph& ug, UG_Attributes& ug_attribs, const Graph& graph, CMap2& cmap2){
        for(Branch branch : graph.branches){
            if(ug.nb_darts_of_orbit(UGVertex(branch.first)) > 1){
                propagate_frame_n_1(ug, ug_attribs, branch.first);
                if(ug.nb_darts_of_orbit(UGVertex(branch.second)) == 1){
                    propagate_frame_n_1(ug, ug_attribs, branch.second);
                }
                else{
                    propagate_frame_n_n(ug, ug_attribs, branch.first, cmap2);
                }
            }
            else{
                propagate_frame_n_1(ug, ug_attribs, branch.second);
            }
        }
        return true;
    }

    bool propagate_frame_n_1(UGraph& ug, UG_Attributes& ug_attribs, Dart d){
        Frame U0 = ug_attribs.frames[d];
        Dart d0 = d;
        Dart d1 = ug.alpha0(d0);
        uint32 valence = ug.nb_darts_of_orbit(UGVertex(d1));
        Frame U, U_;
        U = U0;

        while(valence == 2){
            VEC3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
            Scalar c1, c2;

            Ri = U.col(0);
            Ti = U.col(2);
            v1 = ug_attribs.position[d1] - ug_attribs.position[d0];
            c1 = v1.dot(v1);
            RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
            TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
            // placeholder Ti1
            v3 = (ug_attribs.position[ug.alpha0(ug.alpha1(d1))] - ug_attribs.position[d1]);
            Ti1 = (v1.normalized()+v3.normalized()).normalized();
            v2 = Ti1 - TiL;
            c2 = v2.dot(v2);
            Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
            Si1 = Ti1.cross(Ri1);
            U.col(0) = Ri1;
            U.col(1) = Si1;
            U.col(2) = Ti1;
            /// placeholder


            U_.col(0) = U.col(0);
            U_.col(1) = -U.col(1);
            U_.col(2) = -U.col(2);

            ug_attribs.frames[d1] = U_;
            d0 = ug.alpha1(d1);
            d1 = ug.alpha0(d0);
            ug_attribs.frames[d0] = U;
            valence = ug.nb_darts_of_orbit(UGVertex(d1));
        }

        if(valence == 1){
            VEC3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
            Scalar c1, c2;
            Ri = U.col(0);
            Ti = U.col(2);
            v1 = ug_attribs.position[d1] - ug_attribs.position[d0];
            c1 = v1.dot(v1);
            RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
            TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
            // placeholder Ti1
//            v3 = (ug_attribs.position[ug.alpha0(ug.alpha1(d1))] - ug_attribs.position[d1]);
            Ti1 = v1.normalized();
            v2 = Ti1 - TiL;
            c2 = v2.dot(v2);
            Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
            Si1 = Ti1.cross(Ri1);
            U.col(0) = Ri1;
            U.col(1) = Si1;
            U.col(2) = Ti1;

            U_.col(0) = U.col(0);
            U_.col(1) = -U.col(1);
            U_.col(2) = -U.col(2);
            ug_attribs.frames[d1] = U_;
        }

        return true;
    }

    bool propagate_frame_n_n(UGraph& ug, UG_Attributes& ug_attribs, Dart d, CMap2& cmap2){
        Frame U0 = ug_attribs.frames[d];
        Dart d0 = d;
        Dart d1 = ug.alpha0(d0);
        uint32 valence = ug.nb_darts_of_orbit(UGVertex(d1));
        Frame U, U_;
        U = U0;

        uint nb_e = 0;

        VEC3 v1, v2, v3, Ri, Ti, RiL, TiL, Ri1, Ti1, Si1;
        Scalar c1, c2;

        while(valence == 2){
            nb_e++;

            Ri = U.col(0);
            Ti = U.col(2);
            v1 = ug_attribs.position[d1] - ug_attribs.position[d0];
            c1 = v1.dot(v1);
            RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
            TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;

            v3 = (ug_attribs.position[ug.alpha0(ug.alpha1(d1))] - ug_attribs.position[d1]);
            Ti1 = (v1.normalized()+v3.normalized()).normalized();
            v2 = Ti1 - TiL;
            c2 = v2.dot(v2);
            Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
            Si1 = Ti1.cross(Ri1);
            U.col(0) = Ri1;
            U.col(1) = Si1;
            U.col(2) = Ti1;

            d0 = ug.alpha1(d1);
            d1 = ug.alpha0(d0);
            ug_attribs.frames[d0] = U;
            valence = ug.nb_darts_of_orbit(UGVertex(d1));
        }

        Ri = U.col(0);
        Ti = U.col(2);
        v1 = ug_attribs.position[d1] - ug_attribs.position[d0];
        c1 = v1.dot(v1);
        RiL = Ri - (2/c1)*(v1.dot(Ri))*v1;
        TiL = Ti - (2/c1)*(v1.dot(Ti))*v1;
        Ti1 = v1.normalized();
        v2 = Ti1 - TiL;
        c2 = v2.dot(v2);
        Ri1 = RiL -(2/c2)*(v2.dot(RiL)) * v2;
        Si1 = Ti1.cross(Ri1);
        U.col(0) = Ri1;
        U.col(1) = Si1;
        U.col(2) = Ti1;

        U_.col(0) = U.col(0);
        U_.col(1) = -U.col(1);
        U_.col(2) = -U.col(2);

        /// shift end fram to minimize twist
        VEC3 X = (U_.col(0) + U_.col(1)).normalized();
        Frame UE = ug_attribs.frames[d1];
        VEC3 RE = UE.col(0), SE = UE.col(1);
        bool A = (RE.dot(X) >= 0);
        bool B = (SE.dot(X) >= 0);
        uint nb_shifts = 0;
        if(!A && B) nb_shifts = 1;
        else if(!A && !B) nb_shifts = 2;
        else if(A && !B) nb_shifts = 3;
        if(nb_shifts){
            UE = shift_frame(UE, nb_shifts);
            ug_attribs.frames[d1] = UE;
            ug_attribs.m2_interface[d1] = shift_interface(cmap2, ug_attribs.m2_interface[d1], nb_shifts);
        }

        if(nb_e){
            Scalar angle = std::acos(UE.col(0).dot(U_.col(0)));
            cgogn_log_info("angle: ") << angle;
            Scalar angle_step = angle / nb_e;
            Dart d0 = d;
            Dart d1 = ug.alpha0(d0);
            uint valence = ug.nb_darts_of_orbit(UGVertex(d1));
            uint step = 0;

            while(valence == 2) {
                step++;
                d0 = ug.alpha1(d1);
                U = ug_attribs.frames[d0];
                AngleAxisd rot (angle_step * step, U.col(2));
                U.col(0) = rot * U.col(0);
                U.col(1) = U.col(2).cross(U.col(0));

                U_.col(0) = U.col(0);
                U_.col(1) = -U.col(1);
                U_.col(2) = -U.col(2);

                ug_attribs.frames[d1] = U_;
                ug_attribs.frames[d0] = U;

                d1 = ug.alpha0(d0);
                valence = ug.nb_darts_of_orbit(UGVertex(d1));
            }
        }
        return true;
    }

    bool set_interfaces_geometry(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs){
        ug.foreach_cell([&](UGVertex ugv) -> bool {
            VEC3 center = ug_attribs.position[ugv];
            Scalar radius = ug_attribs.radii[ugv];
            Dart f0 = ug_attribs.m2_interface[ugv.dart];

            m2_attribs.center[M2Volume(f0)] = center;

            if(ug.nb_darts_of_orbit(ugv) > 2){
                return true;
            }

            Frame F = ug_attribs.frames[ugv.dart];

            m2_attribs.position[f0] = center - F.col(1) * radius;
            m2_attribs.position[cmap2.phi<1>(f0)] = center + F.col(0) * radius;
            m2_attribs.position[cmap2.phi<11>(f0)] = center + F.col(1) * radius;
            m2_attribs.position[cmap2.phi<111>(f0)] = center - F.col(0) * radius;

            return true;
        });
        return true;
    }

    bool set_edge_geometry(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs){
        ug.foreach_cell([&](UGVertex ugv){
            VEC3 center = ug_attribs.position[ugv];
            Scalar radius = ug_attribs.radii[ugv];
            Dart cc = ug_attribs.m2_CC[ugv];

            cmap2.foreach_incident_edge(M2Volume(cc), [&](M2Edge m2e){
                Dart d0 = m2e.dart;
                Dart d1 = cmap2.phi2(d0);

                VEC3 edge_pos = project_on_sphere((m2_attribs.position[d0] + m2_attribs.position[d1])/2, radius, center);
                m2_attribs.edge_pos[m2e] = edge_pos;
            });
        });
        return true;
    }

    bool build_branch_sections(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3){
        ug.foreach_cell([&](UGEdge uge){
            Dart e0 = uge.dart;
            Dart e1 = ug.alpha0(e0);
            std::vector<Dart> F = {ug_attribs.m2_interface[e0], ug_attribs.m2_interface[e1]};
            std::vector<Dart> F0 = {F[0], cmap2.phi<1>(F[0]), cmap2.phi<11>(F[0]), cmap2.phi<111>(F[0])};
            std::vector<Dart> F1 = {F[1], cmap2.phi<1>(F[1]), cmap2.phi<11>(F[1]), cmap2.phi<111>(F[1])};

            Dart D = new_section(cmap3);
            std::vector<Dart> D0 = {D, cmap3.phi<2321>(D), cmap3.phi<23212321>(D), cmap3.phi<111232>(D)};
            std::vector<Dart> D1 = {cmap3.phi<2112>(D0[0]), cmap3.phi<2112>(D0[1]), cmap3.phi<2112>(D0[2]), cmap3.phi<2112>(D0[3])};

            m2_attribs.connections[F0[0]] = cmap3.phi<1>(D0[0]);
            m2_attribs.connections[F0[1]] = cmap3.phi<1>(D0[1]);
            m2_attribs.connections[F0[2]] = cmap3.phi<1>(D0[2]);
            m2_attribs.connections[F0[3]] = cmap3.phi<1>(D0[3]);

            m2_attribs.connections[F1[0]] = cmap3.phi<11>(D1[1]);
            m2_attribs.connections[F1[1]] = cmap3.phi<11>(D1[0]);
            m2_attribs.connections[F1[2]] = cmap3.phi<11>(D1[3]);
            m2_attribs.connections[F1[3]] = cmap3.phi<11>(D1[2]);
        });
        return true;
    }

    Dart new_section(CMap3& cmap3){
        M3Builder m3builder(cmap3);
        std::vector<Dart> D = {m3builder.add_prism_topo_fp(4u),
                              m3builder.add_prism_topo_fp(4u),
                              m3builder.add_prism_topo_fp(4u),
                              m3builder.add_prism_topo_fp(4u)
                              };

        m3builder.sew_volumes_fp(cmap3.phi<2>(D[0]), cmap3.phi<1112>(D[1]));
        m3builder.sew_volumes_fp(cmap3.phi<2>(D[1]), cmap3.phi<1112>(D[2]));
        m3builder.sew_volumes_fp(cmap3.phi<2>(D[2]), cmap3.phi<1112>(D[3]));
        m3builder.sew_volumes_fp(cmap3.phi<2>(D[3]), cmap3.phi<1112>(D[0]));
        return D[0];
    }

    bool sew_sections(CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3){
        M3Builder m3builder(cmap3);
        bool valid = true;
        cmap2.foreach_cell([&](M2Edge m2e) -> bool {
            Dart m2d_0 = m2e.dart;
            Dart m2d_1 = cmap2.phi2(m2d_0);

            if(cmap2.is_boundary(m2d_0)||cmap2.is_boundary(m2d_1)){
                return true;
            }

            m3builder.sew_volumes_fp(m2_attribs.connections[m2d_0], cmap3.phi<1>(m2_attribs.connections[m2d_1]));

            return true;
        });
        m3builder.close_map();
        return valid;
    }

    bool set_m3_geometry(CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3){
        CMap3::VertexAttribute<VEC3> position = cmap3.add_attribute<VEC3, M3Vertex>("position");
        M3DartMarker m3Marker(cmap3);

        cmap2.foreach_cell([&](M2Volume m2w){
            Dart m3d = cmap3.phi_1(m2_attribs.connections[m2w.dart]);
            position[m3d] = m2_attribs.center[m2w.dart];
        });

        cmap2.foreach_cell([&](M2Edge m2e){
            Dart d = cmap3.phi1(m2_attribs.connections[m2e.dart]);
            position[d] = m2_attribs.edge_pos[m2e];
        });

        cmap2.foreach_dart([&](Dart m2d){
            if(!cmap2.is_boundary(m2d)){
                Dart m3d = m2_attribs.connections[m2d];
                if(!m3Marker.is_marked(m3d)){
                    position[m3d] = m2_attribs.position[cmap2.phi1(m2d)];
                    m3Marker.mark_orbit(M3Vertex(m3d));
                }
            }
        });
        return true;
    }

    Dart shift_interface(CMap2& cmap2, Dart m2f, uint32 nb_shifts){
        Dart d = m2f;
        for(uint i = 0; i < nb_shifts; ++i)
            d = cmap2.phi1(d);
        return d;
    }

    Frame shift_frame(Frame frame, uint nb_shifts){
        Frame f = frame;
        for(uint i = 0; i < nb_shifts; ++i){
            VEC3 R, S, T;
            R = f.col(1);
            S = -f.col(0);
            T = f.col(2);
            f.col(0) = R;
            f.col(1) = S;
            f.col(2) = T;
        }
        return f;
    }
}




//Vessels_Builder::Vessels_Builder() {}

//void Vessels_Builder::set_skeleton(UGraph* ug){
//    ug_ = ug;
//    ug_marker_ = new UGDartMarker(*ug_);

//    UGTangents_ = ug_->add_attribute<VEC3, UGVertex>("tangents");
//    UGposition_ = ug_->template get_attribute<VEC3, UGVertex>("position");

//}

//void Vessels_Builder::set_cmap3(CMap3* cmap3){
//    cmap3_ = cmap3;
//    m3builder_ = new M3Builder(*cmap3_);
//    if (! cmap3_->is_embedded<M3Vertex>())
//            m3builder_->create_embedding<M3Vertex::ORBIT>();
//    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
//    M3position_ = ca.get_chunk_array<VEC3>("position");
//    if (M3position_ == nullptr)
//        M3position_ = ca.add_chunk_array<VEC3>("position");
//}

//bool Vessels_Builder::compute_cmap3(){
////    subdivide_graph();
//    UGNormals_ = ug_->add_attribute<VEC3, UndirectedGraph::CDart>("normals");
//    UGConnections_ = ug_->template add_attribute<Dart, UndirectedGraph::CDart>("connections");
//    analyse_graph();
//    build_cmap3();
////    clean_up();
//    return true;
//}


////using std::unordered_map;
//void Vessels_Builder::graph_stats(){
//    std::unordered_map<uint32, uint32> vertex_valence;
//    uint32 nb_isolated = 0;
//    ug_->foreach_cell([&](UGVertex ugv){
//        uint32 valence = ug_->nb_darts_of_orbit(ugv);
//        Dart d0, d1;
//        switch(valence){
//            case 1:
//                if (ug_->is_boundary(ug_->alpha0(ugv.dart)))
//                    nb_isolated++;
//                break;
//            case 2:
//                break;
//            case 5:
//            d0 = ug_->alpha1(ugv.dart);
//            cgogn_log_info("ugv: ") << ugv.dart;
//            while(d0 != ugv.dart){
//                cgogn_log_info("dart: ") << d0;
//                d1 = ug_->alpha1(d0);
//                ug_->disconnect_vertices(UGEdge(d0));
//                d0 = d1;
//            }
//            ug_->disconnect_vertices(UGEdge(ugv.dart));
//               ug_->remove_vertex(ugv);
//            break;
//            default:
//                vertex_valence[valence]++;
//                break;
//        }
//    });

//    cgogn_log_info("GRAPH STATS: ") << " ";
//    cgogn_log_info("isolated vertices") << nb_isolated;
//    for(const auto& n : vertex_valence){
//        cgogn_log_info("valence: ") << n.first << " \tnb: " << n.second;
//    }
//}

//void Vessels_Builder::analyse_graph(){
//    find_ends();
//    find_branches();
//    cgogn_log_info("graph analysis:\n") << "nb inter: " << intersections_.size() << "   nb branches: " << branches_.size();
//    uint32 nb_branches = branches_.size();
////    check_isolated_loops();
////    cgogn_log_info("") << "\t\tnb loops: " << branches_.size() - nb_branches;

//}

//void Vessels_Builder::find_ends(){
//    ug_->foreach_cell([&](UGVertex v){
//        switch(ug_->nb_darts_of_orbit(v)){
//            case 1: // end
//                extremities_.push_back(v);
//            case 2: // joint
//            break;
//            case 3: // intersection 3+
//            default:
//                extremities_.push_back(v);
//                intersections_.push_back(v);
//            break;
//        }
//    });
////    cgogn_log_info("extremities: ") << extremities_.size();
//}

//void Vessels_Builder::find_branches(){
//    for(UGVertex v : extremities_){
//        ug_->foreach_dart_of_orbit(v, [&](Dart d0){
////            Dart d0 = v.dart;
//            if(!ug_marker_->is_marked(d0)){
//                ug_marker_->mark(d0);
//                Dart d1 = ug_->alpha0(d0);
//                ug_marker_->mark(d1);
//                while(!ug_->is_boundary(d1) && ug_->nb_darts_of_orbit(UGVertex(d1)) == 2){
//                    d1 = ug_->alpha1(d1);
//                    ug_marker_->mark(d1);
//                    d1 = ug_->alpha0(d1);
//                    ug_marker_->mark(d1);
//                }
//                if(!ug_->is_boundary(d1)) // pas un sommet isolé
//                    branches_.push_back(std::pair<Dart, Dart>(d0, d1));
////                cgogn_log_info("extremities: ") << d0 << "-" << ug_->nb_darts_of_orbit(UGVertex(d0)) << "-" << ug_->embedding(d0, UGVertex::ORBIT)
////                        << " / " << d1 << "-" << ug_->nb_darts_of_orbit(UGVertex(d1)) << "-" << ug_->embedding(d1, UGVertex::ORBIT) ;
//            }
//        });
//    }
//}

//void Vessels_Builder::check_isolated_loops(){
//    ug_->foreach_dart([&](Dart d0){
//        if(!ug_marker_->is_marked(d0)){
//            ug_marker_->mark(d0);
//            Dart d1 = ug_->alpha0(d0);
//            ug_marker_->mark(d1);
//            while(!(ug_->alpha1(d1) == d0) && !ug_->is_boundary(d1) && ug_->nb_darts_of_orbit(UGVertex(d1)) == 2){
//                d1 = ug_->alpha1(d1);
//                ug_marker_->mark(d1);
//                d1 = ug_->alpha0(d1);
//                ug_marker_->mark(d1);
//            }
//            if(!ug_->is_boundary(d1))
//                branches_.push_back(std::pair<Dart, Dart>(d0, d1));
//        }
//    });
//}

//void Vessels_Builder::subdivide_graph(){
////    ug_->nb_cells();
//    Scalar radius(0.5);
//    ug_marker_->unmark_all();
//    bool over = false;
//    while(!over){
//        over = true;
//        ug_->foreach_cell([&](UGEdge e){
//            if(!ug_marker_->is_marked(e.dart)){
//                Dart v0 = e.dart;
//                Dart v1 = ug_->alpha0(e.dart);
//                VEC3 edge_vec = UGposition_[v0] - UGposition_[v1];
////                cgogn_log_info("edge: ") << edge_vec.norm();
//                if(edge_vec.norm() > 2 * radius){
//                    over &= false;
//                    UGVertex vert = ug_->add_vertex();
//                    UGposition_[vert] = (UGposition_[v0] + UGposition_[v1]) * Scalar(0.5) ;
////                    cgogn_log_info("length longer than radius: ") << 1;
//                    ug_->connect_vertices(UGVertex(v0), vert);
//                    ug_->connect_vertices(UGVertex(v1), vert);
//                    ug_->disconnect_vertices(e);
//                }
//                else{
//                    ug_marker_->mark_orbit(e);
//                }
//            }
//        });
//    }
//    ug_marker_->unmark_all();
//}

//void Vessels_Builder::build_cmap3(){
//    cgogn_log_info("Building cmap3");
//    build_intersections();
//    compute_tangents();
//    build_branches();
//    m3builder_->close_map();
////    cmap3_->check_map_integrity();
//}

//void Vessels_Builder::build_intersections(){
//    intersection_m2builder_.set_cmap2(cmap2_);
//    intersection_m2builder_.set_ugraph(ug_);

//    for(UGVertex ugv : intersections_){
//        if(ug_->nb_darts_of_orbit(UGVertex(ugv)) == 3)
//            build_intersection3(ugv);
//        else
////        if(ug_->nb_darts_of_orbit(UGVertex(ugv)) == 5)
//            build_intersectionN(ugv);
//    }
//}

//void Vessels_Builder::build_intersection3(UGVertex ugv){
//    Scalar radius = Scalar(1.5f);
//    VEC3 Ctr = UGposition_[ugv];
//    std::vector<VEC3> P, PV, C, Q, M, E, F;
//    std::vector<Dart> PD;

////    cgogn_log_info("INTERSECTION: ") << ug_->embedding(ugv.dart, UGVertex::ORBIT);
//    ug_->foreach_adjacent_vertex_through_edge(ugv, [&](UGVertex ugv2){
//        VEC3 v = (UGposition_[ugv2] - Ctr).normalized();
//        ug_->for
////        cgogn_log_info("d") << ug_->embedding(ugv2.dart, UGVertex::ORBIT) ;
//        P.push_back(Ctr + v * Scalar(radius));
//        PD.push_back(ug_->alpha0(ugv2.dart));
//    });

//    VEC3 V = -(P[1] - P[0]).cross(P[2] - P[0]).normalized();
//    Q.push_back(Ctr + V * Scalar(radius));
//    Q.push_back(Ctr - V * Scalar(radius));

//    PV.push_back((P[1]-P[0]).normalized());
//    PV.push_back((P[2]-P[1]).normalized());
//    PV.push_back((P[0]-P[2]).normalized());
//    M.push_back(Ctr - PV[0].cross(V) * radius);
//    M.push_back(Ctr - PV[1].cross(V) * radius);
//    M.push_back(Ctr - PV[2].cross(V) * radius);

//    C.push_back((M[0]+M[1]+M[2]+Q[0])*Scalar(0.25f));
//    C.push_back((M[0]+M[1]+M[2]+Q[1])*Scalar(0.25f));

//    E.push_back((M[0] + M[1])*Scalar(0.5f));
//    E.push_back((M[1] + M[2])*Scalar(0.5f));
//    E.push_back((M[2] + M[0])*Scalar(0.5f));
//    E.push_back((M[0] + Q[0])*Scalar(0.5f));
//    E.push_back((M[1] + Q[0])*Scalar(0.5f));
//    E.push_back((M[2] + Q[0])*Scalar(0.5f));
//    E.push_back((M[0] + Q[1])*Scalar(0.5f));
//    E.push_back((M[1] + Q[1])*Scalar(0.5f));
//    E.push_back((M[2] + Q[1])*Scalar(0.5f));

//    F.push_back((M[0] + M[1] + Q[0]) / Scalar(3.0f));
//    F.push_back((M[1] + M[2] + Q[0]) / Scalar(3.0f));
//    F.push_back((M[2] + M[0] + Q[0]) / Scalar(3.0f));
//    F.push_back((M[0] + M[1] + Q[1]) / Scalar(3.0f));
//    F.push_back((M[1] + M[2] + Q[1]) / Scalar(3.0f));
//    F.push_back((M[2] + M[0] + Q[1]) / Scalar(3.0f));

//    /// Projection sur la shère de centre Ctr de l'embranchement
////    for(uint32 i = 0; i < E.size(); i++){
////        E[i] = (E[i] - Ctr).normalized();
////        E[i] = Ctr + E[i] * radius;
////    }

////    for(uint32 i = 0; i < F.size(); i++){
////        F[i] = (F[i] - Ctr).normalized();
////        F[i] = Ctr + F[i] * radius;
////    }

//    /// building CMAP3
//    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
//    /// Insertion des points dans la 3 carte
//    std::vector<uint32> Ci, Qi, Mi, Ei, Fi;
//    uint32 Ctri;

//    uint32 l = ca.insert_lines<1>();
//    (*M3position_)[l] = Ctr;
//    Ctri = l;

//    for(VEC3 v: C){
//        uint32 l = ca.insert_lines<1>();
//        (*M3position_)[l] = v;
//        Ci.push_back(l);
//    }
//    C.clear();

//    for(VEC3 v: Q){
//        uint32 l = ca.insert_lines<1>();
//        (*M3position_)[l] = v;
//        Qi.push_back(l);
//    }
//    Q.clear();

//    for(VEC3 v: M){
//        uint32 l = ca.insert_lines<1>();
//        (*M3position_)[l] = v;
//        Mi.push_back(l);
//    }
//    M.clear();

//    for(VEC3 v: E){
//        uint32 l = ca.insert_lines<1>();
//        (*M3position_)[l] = v;
//        Ei.push_back(l);
//    }
//    E.clear();

//    for(VEC3 v: F){
//        uint32 l = ca.insert_lines<1>();
//        (*M3position_)[l] = v;
//        Fi.push_back(l);
//    }
//    F.clear();

//    /// Construction des hexaèdres de l'embranchement
//    std::vector<Dart> H;
//    H.push_back(add_hexa({Mi[0], Ei[2], Fi[2], Ei[3], Ei[0], Ctri, Ci[0], Fi[0]}));
//    H.push_back(add_hexa({Ei[0], Ctri, Ci[0], Fi[0], Mi[1], Ei[1], Fi[1], Ei[4]}));
//    H.push_back(add_hexa({Ctri, Ei[2], Fi[2], Ci[0], Ei[1], Mi[2], Ei[5], Fi[1]}));
//    H.push_back(add_hexa({Ei[3], Fi[2], Ei[5], Qi[0], Fi[0], Ci[0], Fi[1], Ei[4]}));
//    H.push_back(add_hexa({Ei[2], Mi[0], Ei[6], Fi[5], Ctri, Ei[0], Fi[3], Ci[1]}));
//    H.push_back(add_hexa({Ei[0], Ctri, Ei[1], Mi[1], Fi[3], Ci[1], Fi[4], Ei[7]}));
//    H.push_back(add_hexa({Ctri, Ei[2], Mi[2], Ei[1], Ci[1], Fi[5], Ei[8], Fi[4]}));
//    H.push_back(add_hexa({Fi[5], Ei[6], Qi[1], Ei[8], Ci[1], Fi[3], Ei[7], Fi[4]}));

//    /// couture topologique des hexaèdres de l'embranchement
//    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[0]), H[1]);
//    m3builder_->sew_volumes_fp(cmap3_->phi<121>(H[0]), H[2]);
//    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[0]), cmap3_->phi<2>(H[3]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[0]), cmap3_->phi<2>(H[4]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[1]), cmap3_->phi<1112>(H[2]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[1]), cmap3_->phi<2112>(H[3]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[1]), H[5]);
//    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[2]), cmap3_->phi<121>(H[3]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[2]), H[6]);
//    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[4]), cmap3_->phi<2>(H[5]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<212>(H[4]), cmap3_->phi<2>(H[6]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[4]), cmap3_->phi<2>(H[7]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[5]), cmap3_->phi<1112>(H[6]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[5]), cmap3_->phi<2112>(H[7]));
//    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[6]), cmap3_->phi<212>(H[7]));

//    /// Mise en place de connections à l'embranchement dans l'attribut
//    UGConnections_[PD[0]] = cmap3_->phi<123211>(H[0]);
//    UGConnections_[PD[1]] = cmap3_->phi<2121>(H[0]);
//    UGConnections_[PD[2]] = cmap3_->phi<21121>(H[1]);

//    /// stockage de la normal
////    UGNormals_[ugv] = V;

//    UGNormals_[PD[1]] = V;
//    UGNormals_[PD[2]] = V;
//    UGNormals_[PD[0]] = V;

//}

//void Vessels_Builder::build_intersectionN(UGVertex ugv){
////    CMap2* map2 = new CMap2();
//    intersection_m2builder_.build_intersection(ugv);
//    CMap2::VertexAttribute<VEC3> M2Position = cmap2_->get_attribute<VEC3, M2Vertex>("position");
//    CMap2::FaceAttribute<Dart> M2FaceBranch = cmap2_->get_attribute<Dart, M2Face>("face_branches");
//    CMap2::CDartAttribute<Dart> M2BlocConnections = cmap2_->add_attribute<Dart, CMap2::CDart>("connections");
//    CMap2::VertexAttribute<uint32> M2_EmbVert = cmap2_->add_attribute<uint32, M2Vertex>("emb_vert");
//    CMap2::VertexAttribute<uint32> M2_EmbMid2Ctr = cmap2_->add_attribute<uint32, M2Vertex>("emb_mid2ctr");
//    CMap2::EdgeAttribute<uint32> M2_EmbMidEdge = cmap2_->add_attribute<uint32, M2Edge>("emb_midedge");
//    CMap2::EdgeAttribute<uint32> M2_EmbFace2Ctr = cmap2_->add_attribute<uint32, M2Edge>("emb_face2ctr");

//    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();

//    uint32 Ctri = ca.insert_lines<1>();
//    (*M3position_)[Ctri] = UGposition_[ugv];

//    cmap2_->foreach_cell([&](M2Vertex m2v){
//        M2_EmbMid2Ctr[m2v] = ca.insert_lines<1>();
//        (*M3position_)[M2_EmbMid2Ctr[m2v]] = (UGposition_[ugv] + M2Position[m2v])/2;
//        M2_EmbVert[m2v] = ca.insert_lines<1>();
//        (*M3position_)[M2_EmbVert[m2v]] = M2Position[m2v];
//    });

//    cmap2_->foreach_cell([&](M2Edge m2e){
//        Dart v0 = m2e.dart;
//        Dart v1 = cmap2_->phi2(v0);
//        M2_EmbMidEdge[m2e] = ca.insert_lines<1>();
//        (*M3position_)[M2_EmbMidEdge[m2e]] = (M2Position[v0] + M2Position[v1])/2;
//        M2_EmbFace2Ctr[m2e] = ca.insert_lines<1>();
//        (*M3position_)[M2_EmbFace2Ctr[m2e]] = (UGposition_[ugv] + M2Position[v0] + M2Position[v1])/3;
////        cgogn_log_info("embmid2ctr: ") << M2_EmbMidEdge[m2e];
//    });

//    cmap2_->foreach_cell([&](M2Face f){
//        Dart fd = intersection_m2builder_.convex_quad(f.dart);
//        Scalar radius = Scalar(1.5f);

//        std::vector<uint32> Ci, Qi, Mi, Ei, Fi;
//        Ci = {ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>()};
//        Qi = {M2_EmbVert[cmap2_->phi_1(fd)], M2_EmbVert[cmap2_->phi1(fd)]};
//        Mi = {M2_EmbVert[cmap2_->phi<11>(fd)], M2_EmbVert[fd], Ctri};
//        Ei = {ca.insert_lines<1>(), M2_EmbMid2Ctr[fd], M2_EmbMid2Ctr[cmap2_->phi<11>(fd)],
//             M2_EmbMidEdge[cmap2_->phi<11>(fd)], M2_EmbMidEdge[cmap2_->phi_1(fd)], M2_EmbMid2Ctr[cmap2_->phi_1(fd)],
//             M2_EmbMidEdge[cmap2_->phi1(fd)], M2_EmbMidEdge[fd], M2_EmbMid2Ctr[cmap2_->phi1(fd)]};
//        Fi = {ca.insert_lines<1>(), M2_EmbFace2Ctr[cmap2_->phi_1(fd)], M2_EmbFace2Ctr[cmap2_->phi<11>(fd)],
//             ca.insert_lines<1>(), M2_EmbFace2Ctr[fd], M2_EmbFace2Ctr[cmap2_->phi1(fd)]};

//        (*M3position_)[Ei[0]] = (M2Position[fd] + M2Position[cmap2_->phi<11>(fd)])/2;
//        (*M3position_)[Fi[0]] = (M2Position[fd] + M2Position[cmap2_->phi_1(fd)] + M2Position[cmap2_->phi<11>(fd)])/3;
//        (*M3position_)[Fi[3]] = (M2Position[fd] + M2Position[cmap2_->phi1(fd)] + M2Position[cmap2_->phi<11>(fd)])/3;
//        (*M3position_)[Ci[0]] = (M2Position[fd] + M2Position[cmap2_->phi<11>(fd)] + UGposition_[ugv])/3;
//        (*M3position_)[Ci[1]] = (M2Position[fd] + M2Position[cmap2_->phi_1(fd)] + M2Position[cmap2_->phi<11>(fd)] + UGposition_[ugv])/4;
//        (*M3position_)[Ci[2]] = (M2Position[fd] + M2Position[cmap2_->phi1(fd)] + M2Position[cmap2_->phi<11>(fd)] + UGposition_[ugv])/4;

//        std::vector<Dart> H;
//        H.push_back(add_hexa({Mi[0], Ei[2], Fi[2], Ei[3], Ei[0], Ci[0], Ci[1], Fi[0]}));
//        H.push_back(add_hexa({Ei[0], Ci[0], Ci[1], Fi[0], Mi[1], Ei[1], Fi[1], Ei[4]}));
//        H.push_back(add_hexa({Ci[0], Ei[2], Fi[2], Ci[1], Ei[1], Mi[2], Ei[5], Fi[1]}));
//        H.push_back(add_hexa({Ei[3], Fi[2], Ei[5], Qi[0], Fi[0], Ci[1], Fi[1], Ei[4]}));
//        H.push_back(add_hexa({Ei[2], Mi[0], Ei[6], Fi[5], Ci[0], Ei[0], Fi[3], Ci[2]}));
//        H.push_back(add_hexa({Ei[0], Ci[0], Ei[1], Mi[1], Fi[3], Ci[2], Fi[4], Ei[7]}));
//        H.push_back(add_hexa({Ci[0], Ei[2], Mi[2], Ei[1], Ci[2], Fi[5], Ei[8], Fi[4]}));
//        H.push_back(add_hexa({Fi[5], Ei[6], Qi[1], Ei[8], Ci[2], Fi[3], Ei[7], Fi[4]}));

//        m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[0]), H[1]);
//        m3builder_->sew_volumes_fp(cmap3_->phi<121>(H[0]), H[2]);
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[0]), cmap3_->phi<2>(H[3]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[0]), cmap3_->phi<2>(H[4]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[1]), cmap3_->phi<1112>(H[2]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[1]), cmap3_->phi<2112>(H[3]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[1]), H[5]);
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[2]), cmap3_->phi<121>(H[3]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[2]), H[6]);
//        m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[4]), cmap3_->phi<2>(H[5]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<212>(H[4]), cmap3_->phi<2>(H[6]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[4]), cmap3_->phi<2>(H[7]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[5]), cmap3_->phi<1112>(H[6]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[5]), cmap3_->phi<2112>(H[7]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[6]), cmap3_->phi<212>(H[7]));

////        UGConnections_[M2FaceBranch[f]] = cmap3_->phi<2121>(H[0]);
//        UGConnections_[M2FaceBranch[f]] = cmap3_->phi<11121>(H[5]);
//        UGNormals_[M2FaceBranch[f]] = ((*M3position_)[Ci[2]] - (*M3position_)[Ci[1]]).normalized();

//        M2BlocConnections[fd] = cmap3_->phi<112111>(H[5]);
//        M2BlocConnections[cmap2_->phi<1>(fd)] = cmap3_->phi<1>(H[7]);
//        M2BlocConnections[cmap2_->phi<11>(fd)] = cmap3_->phi<111>(H[0]);
//        M2BlocConnections[cmap2_->phi<111>(fd)] = cmap3_->phi<112111>(H[3]);

//    });

//    uint32 n = 0;
//    cmap2_->foreach_cell([&](M2Edge e){
//        n++;
//        Dart d0 = M2BlocConnections[e.dart];
//        Dart d1 = M2BlocConnections[cmap2_->phi2(e.dart)];
//        m3builder_->sew_volumes_fp(cmap3_->phi<111232>(d0), cmap3_->phi<111>(d1));
//        m3builder_->sew_volumes_fp(cmap3_->phi<111>(d0), cmap3_->phi<111232>(d1));
//        m3builder_->sew_volumes_fp(cmap3_->phi<11232>(d0), cmap3_->phi<11232111>(d1));
//    });



//}

//void Vessels_Builder::build_branches(){
////    uint i = 0;
////    for(std::pair<Dart, Dart> branch : branches_){
////        cgogn_log_info("branches: ") << i++ << " " << ug_->embedding(branch.first, UGVertex::ORBIT) << " " << ug_->embedding(branch.second, UGVertex::ORBIT);
////    }

//    ug_marker_->unmark_all();
//    for(std::pair<Dart, Dart> branch : branches_){
//        if(ug_->nb_darts_of_orbit(UGVertex(branch.first)) >= 3){
//            build_branch(branch.first);
//        }
//        else{
//            if(ug_->nb_darts_of_orbit(UGVertex(branch.second)) > 1)
//            build_branch(branch.second);
//        }
//    }
//}

//void Vessels_Builder::compute_tangents(){
//    ug_marker_->unmark_all();
//    for(std::pair<Dart, Dart> branch : branches_){
//        if(ug_->nb_darts_of_orbit(UGVertex(branch.first)) >= 3){
//            compute_tangent_3_1(branch.first);
//        }
//        else{
//            compute_tangent_3_1(branch.second);
//        }
//    }
//}

//void Vessels_Builder::compute_tangent_3_1(Dart d){
////    cgogn_log_info("branch 3 1");
//    VEC3 N = UGNormals_[d];
//    Dart d0 = d;
//    Dart d1 = ug_->alpha0(d);
//    VEC3 V0 = (UGposition_[d1] - UGposition_[d0]).normalized();
//    uint32 valence = ug_->nb_darts_of_orbit(UGVertex(d1));
//    while(valence == 2){
////        cgogn_log_info("darts: ") << d0 << " " << d1;
//        UGNormals_[d1] = N;
//        d0 = ug_->alpha1(d1);
//        d1 = ug_->alpha0(d0);
//        VEC3 V1 = (UGposition_[d1] - UGposition_[d0]).normalized();
//        UGNormals_[d0] = N;
//        VEC3 T = (V1 + V0).normalized();
//        UGTangents_[d0] = T;
//        V0 = V1;
//        valence = ug_->nb_darts_of_orbit(UGVertex(d1));
//    }

//    if(valence == 1){
//        UGNormals_[d1] = N;
//        UGTangents_[d1] = V0;
//    }
//}

//void Vessels_Builder::build_branch(Dart d){
////    cgogn_log_info("looping1");
//    Scalar radius = Scalar(1.5f);
//    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();

//    Dart d0 = d;
//    Dart d1 = ug_->alpha0(d0);
//    uint32 valence = ug_->nb_darts_of_orbit(UGVertex(d1));
//    uint32 vertex = ug_->embedding(d, UGVertex::ORBIT);
//    while(true){
////        cgogn_log_info("valence: ") << valence << " " << d0 << " " << d1 ;
//        VEC3 T = UGTangents_[d1];
//        VEC3 N = UGNormals_[d1];
//        VEC3 M = (N.cross(T)).normalized();

//        Dart D = UGConnections_[d0];
//        if(D.is_nil()){
//            cgogn_log_info("Branch error: connection is nil") << d0 << " " << D;
//            return;
//        }

//        // entrance connection darts and
//        std::vector<Dart> D0 = {D, cmap3_->phi<1232>(D), cmap3_->phi<11232>(D), cmap3_->phi<111232>(D), cmap3_->phi<1112321232>(D), cmap3_->phi<11123211232>(D)};
//        std::vector<uint32> E0 = {
//            cmap3_->embedding(D0[2], M3Vertex::ORBIT),
//            cmap3_->embedding(D0[1], M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<1>(D0[0]), M3Vertex::ORBIT),
//            cmap3_->embedding(D0[0], M3Vertex::ORBIT),
//            cmap3_->embedding(D0[5], M3Vertex::ORBIT),
//            cmap3_->embedding(D0[4], M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<111>(D0[5]), M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<11>(D0[5]), M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<11>(D0[4]), M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<11>(D0[2]), M3Vertex::ORBIT),
//            cmap3_->embedding(cmap3_->phi<11>(D0[1]), M3Vertex::ORBIT)
//        };

//        std::vector<Dart> D1 = {add_hexa(), add_hexa(), add_hexa(), add_hexa(), add_hexa(), add_hexa()};
//        m3builder_->sew_volumes_fp(cmap3_->phi<1112>(D1[0]), cmap3_->phi<2>(D1[1]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(D1[0]), cmap3_->phi<2>(D1[2]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[0]), cmap3_->phi<2>(D1[3]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[1]), cmap3_->phi<1112>(D1[2]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[2]), cmap3_->phi<1112>(D1[4]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<1112>(D1[3]), cmap3_->phi<2>(D1[4]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<112>(D1[3]), cmap3_->phi<2>(D1[5]));
//        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[4]), cmap3_->phi<1112>(D1[5]));

//        m3builder_->sew_volumes_fp(D0[0], D1[0]);
//        m3builder_->sew_volumes_fp(D0[1], D1[1]);
//        m3builder_->sew_volumes_fp(D0[2], D1[2]);
//        m3builder_->sew_volumes_fp(D0[3], D1[3]);
//        m3builder_->sew_volumes_fp(D0[4], D1[4]);
//        m3builder_->sew_volumes_fp(D0[5], D1[5]);

//        std::vector<uint32> E1;
//        /// gestion des extrémités et articulations
//        if(valence <= 2){
//            E1 = {ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(),
//                  ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>()
//           };
//            (*M3position_)[E1[0]] = UGposition_[d1];
//            (*M3position_)[E1[1]] = UGposition_[d1] + N * radius * Scalar(0.5f);
//            (*M3position_)[E1[2]] = UGposition_[d1] + (N + M).normalized() * radius;
//            (*M3position_)[E1[3]] = UGposition_[d1] + M * radius;
//            (*M3position_)[E1[4]] = UGposition_[d1] - (N - M).normalized() * radius;
//            (*M3position_)[E1[5]] = UGposition_[d1] - N * radius * Scalar(0.5f);
//            (*M3position_)[E1[6]] = UGposition_[d1] - N * radius;
//            (*M3position_)[E1[7]] = UGposition_[d1] - (N + M).normalized() * radius;
//            (*M3position_)[E1[8]] = UGposition_[d1] - M * radius;
//            (*M3position_)[E1[9]] = UGposition_[d1] + (N - M).normalized() * radius;
//            (*M3position_)[E1[10]] = UGposition_[d1] + N * radius;
//        }
//        /// gestion des embranchements
//        if(valence >= 3){
//            Dart E = UGConnections_[d1];
////            cgogn_log_info("Branching point") << E ;
//            if(E.is_nil()){
//                cgogn_log_info("Branch error: connection is nil") << E ;
//            }

//            /// Inversion des connections si les normales sont trop différentes
////            if(UGNormals_[UGVertex(d0)].dot(UGNormals_[UGVertex(d1)]) >= 0)
//            if(UGNormals_[d0].dot(UGNormals_[d1]) >= 0)
//                E = cmap3_->phi<1123211>(E);
//            else
//                E = cmap3_->phi<111232111>(E);

////            std::vector<Dart> D2 = {E, cmap3_->phi<111232>(E), cmap3_->phi<11232>(E), cmap3_->phi<1232>(E), cmap3_->phi<1232111232>(E), cmap3_->phi<123211232>(E)};
//            std::vector<Dart> D2 = {E, cmap3_->phi<111232>(E), cmap3_->phi<11232>(E), cmap3_->phi<1232>(E), cmap3_->phi<1232111232>(E), cmap3_->phi<123211232>(E)};
//            E1 = {
//                cmap3_->embedding(D2[3], M3Vertex::ORBIT),
//                cmap3_->embedding(D2[2], M3Vertex::ORBIT),
//                cmap3_->embedding(D2[0], M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<1>(D2[0]), M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<2>(D2[5]), M3Vertex::ORBIT),
//                cmap3_->embedding(D2[5], M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<11>(D2[5]), M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<12>(D2[4]), M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<11>(D2[2]), M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<11>(D2[1]), M3Vertex::ORBIT),
//                cmap3_->embedding(cmap3_->phi<111>(D2[1]), M3Vertex::ORBIT),
//            };
//            m3builder_->sew_volumes_fp(D2[0], cmap3_->phi<2112>(D1[0]));

//            m3builder_->sew_volumes_fp(D2[1], cmap3_->phi<2112>(D1[1]));
//            m3builder_->sew_volumes_fp(D2[2], cmap3_->phi<2112>(D1[2]));
//            m3builder_->sew_volumes_fp(D2[3], cmap3_->phi<2112>(D1[3]));
//            m3builder_->sew_volumes_fp(D2[4], cmap3_->phi<2112>(D1[4]));
//            m3builder_->sew_volumes_fp(D2[5], cmap3_->phi<2112>(D1[5]));
//            uint32 vertex2 = ug_->embedding(d1, UGVertex::ORBIT);
////            cgogn_log_info("valence: ") << valence <<" " << vertex << " " << vertex2;
//        }

//        embed_hexa(D1[0], {E0[2], E0[3], E0[0], E0[1], E1[2], E1[3], E1[0], E1[1]});
//        embed_hexa(D1[1], {E0[2], E0[1], E0[9], E0[10], E1[2], E1[1], E1[9], E1[10]});
//        embed_hexa(D1[2], {E0[1], E0[0], E0[8], E0[9], E1[1], E1[0], E1[8], E1[9]});
//        embed_hexa(D1[3], {E0[0], E0[3], E0[4], E0[5], E1[0], E1[3], E1[4], E1[5]});
//        embed_hexa(D1[4], {E0[0], E0[5], E0[7], E0[8], E1[0], E1[5], E1[7], E1[8]});
//        embed_hexa(D1[5], {E0[5], E0[4], E0[6], E0[7], E1[5], E1[4], E1[6], E1[7]});

//        if(valence == 1) break;
//        if(valence == 2){
//            d0 = ug_->alpha1(d1);
//            d1 = ug_->alpha0(d0);
//        }
//        if(valence >= 3) break;


//        UGConnections_[d0] = cmap3_->phi<2112>(D1[0]);
//        valence = ug_->nb_darts_of_orbit(UGVertex(d1));
//    }

//}

//cgogn::Dart Vessels_Builder::add_hexa(){
//    return m3builder_->add_prism_topo_fp(4u);
//}

//cgogn::Dart Vessels_Builder::add_hexa(std::vector<uint32> verts){
//    Dart d = m3builder_->add_prism_topo_fp(4u);
//    const std::array<Dart, 8> vertices_of_hexa = {
//        d,
//        cmap3_->phi1(d),
//        cmap3_->phi1(cmap3_->phi1(d)),
//        cmap3_->phi_1(d),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi_1(d))))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(d)))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(d))))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(d))))))
//    };

//    uint32 index = 0;
//    for (Dart dv : vertices_of_hexa)
//        m3builder_->template set_orbit_embedding<M3Vertex>(M3Vertex2(dv), verts[index++]);
//    return d;
//}

//void Vessels_Builder::embed_hexa(Dart d, std::vector<uint32> verts){
//    const std::array<Dart, 8> vertices_of_hexa = {
//        d,
//        cmap3_->phi1(d),
//        cmap3_->phi1(cmap3_->phi1(d)),
//        cmap3_->phi_1(d),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi_1(d))))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(d)))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(d))))),
//        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(d))))))
//    };

//    uint32 index = 0;
//    for (Dart dv : vertices_of_hexa)
//        m3builder_->template set_orbit_embedding<M3Vertex>(M3Vertex2(dv), verts[index++]);
//}

//void Vessels_Builder::clean_up(){
//    ug_->remove_attribute(UGConnections_);
//    ug_->remove_attribute(UGTangents_);
//    ug_->remove_attribute(UGNormals_);
//}


//Intersection_M2Builder::Intersection_M2Builder() : radius_(1.5), cmap2_(nullptr), ug_(nullptr){

//}

//void Intersection_M2Builder::set_cmap2(CMap2 *cmap2){
//    cmap2_ = cmap2;
//    cmap2_->clear_and_remove_attributes();

//    m2builder_ = new M2Builder(*cmap2_);
//}

//void Intersection_M2Builder::set_ugraph(UGraph *ug){
//    ug_ = ug;
//    UGposition_ = ug_->template get_attribute<VEC3, UGVertex>("position");
//}

//void Intersection_M2Builder::build_intersection(UGVertex ugv){
//    cmap2_->clear_and_remove_attributes();

//    Ppos_.clear();
//    Pdart_.clear();

//    ugv_ = ugv;
//    center_ = UGposition_[ugv];
//    radius_ = Scalar(1.5); /// PLACEHOLDER

//    get_intersection_data();
//    build_all();
//}

//bool Intersection_M2Builder::in_quad(Dart face, VEC3 P){
//    using cgogn::geometry::intersection_ray_triangle;

//    VEC3 A = M2Position_[face];
//    VEC3 B = M2Position_[cmap2_->phi<1>(face)];
//    VEC3 C = M2Position_[cmap2_->phi<11>(face)];
//    VEC3 D = M2Position_[cmap2_->phi<111>(face)];
//    VEC3 DIR = P - center_;

//    return (intersection_ray_triangle(center_, DIR, A, B, C)
//            || intersection_ray_triangle(center_, DIR, A, C, D));
//}

//void Intersection_M2Builder::get_intersection_data(){
//    ug_->foreach_adjacent_vertex_through_edge(ugv_, [&](UGVertex ugv){
//        Ppos_.push_back(project_on_sphere(UGposition_[ugv]));
//        Pdart_.push_back(ug_->alpha0(ugv.dart));
//    });
//}

//void Intersection_M2Builder::build_core(){
//    std::vector<VEC3> Qp, Mp;
//    std::vector<Dart> Fd;

//    Fd = {m2builder_->add_face_topo_fp(4), m2builder_->add_face_topo_fp(4), m2builder_->add_face_topo_fp(4)};

//    m2builder_->phi2_sew(Fd[0], cmap2_->phi<1>(Fd[2]));
//    m2builder_->phi2_sew(cmap2_->phi<1>(Fd[0]), Fd[1]);
//    m2builder_->phi2_sew(cmap2_->phi<11>(Fd[0]), cmap2_->phi_1(Fd[1]));
//    m2builder_->phi2_sew(cmap2_->phi_1(Fd[0]), cmap2_->phi<11>(Fd[2]));
//    m2builder_->phi2_sew(cmap2_->phi<1>(Fd[1]), Fd[2]);
//    m2builder_->phi2_sew(cmap2_->phi<11>(Fd[1]), cmap2_->phi_1(Fd[2]));

//    M2Position_ = cmap2_->add_attribute<VEC3, M2Vertex>("position");
//    M2FacePoints_ = cmap2_->add_attribute<VEC3, M2Face>("face_points");
//    M2FaceBranch_ = cmap2_->add_attribute<Dart, M2Face>("face_branches");

//    VEC3 V = (Ppos_[1] - Ppos_[0]).cross(Ppos_[2] - Ppos_[0]).normalized();
//    Qp = {center_ + V * radius_, center_ - V * radius_};
//    Mp = {center_ + (Ppos_[1]-Ppos_[0]).normalized().cross(V) * radius_,
//         center_ + (Ppos_[2]-Ppos_[1]).normalized().cross(V) * radius_,
//         center_ + (Ppos_[0]-Ppos_[2]).normalized().cross(V) * radius_};

//    M2Position_[M2Vertex(Fd[0])] = Mp[2];
//    M2Position_[M2Vertex(Fd[1])] = Mp[0];
//    M2Position_[M2Vertex(Fd[2])] = Mp[1];
//    M2Position_[M2Vertex(cmap2_->phi<1>(Fd[0]))] = Qp[1];
//    M2Position_[M2Vertex(cmap2_->phi_1(Fd[0]))] = Qp[0];

//    M2FacePoints_[M2Face(Fd[0])] = Ppos_[0];
//    M2FacePoints_[M2Face(Fd[1])] = Ppos_[1];
//    M2FacePoints_[M2Face(Fd[2])] = Ppos_[2];

//    M2FaceBranch_[M2Face(Fd[0])] = Pdart_[0];
//    M2FaceBranch_[M2Face(Fd[1])] = Pdart_[1];
//    M2FaceBranch_[M2Face(Fd[2])] = Pdart_[2];
//}

//void Intersection_M2Builder::build_all(){
////    cgogn_log_info("valence N") << Ppos_.size();
//    build_core();
////    move_points();
//    for(uint i = 3; i < Ppos_.size(); i++){
//        add_point(i);
//    }
//    move_points();
////    cmap2_->create_embedding<M2Edge::ORBIT>();
//}

//void Intersection_M2Builder::add_point(uint32 pt_nb){
//    VEC3 P0 = Ppos_[pt_nb];
//    Dart F0 = Pdart_[pt_nb];
//    M2Face face2cut;
//    std::vector<VEC3> Quadp; //positions
//    std::vector<Dart> Quadd; //darts

//    bool face_found = false;
//    cmap2_->foreach_cell([&](M2Face f) -> bool {
//        if((face_found = in_quad(f.dart, P0))){
//            Quadd = {f.dart, cmap2_->phi<1>(f.dart), cmap2_->phi<11>(f.dart), cmap2_->phi<111>(f.dart)};
//            Quadp = {M2Position_[Quadd[0]], M2Position_[Quadd[1]], M2Position_[Quadd[2]], M2Position_[Quadd[3]]};
//            face2cut = f; return false;
//        }
//        return true;
//    });
//    if(!face_found)
//        cgogn_log_info("no face found") << face_found;

//    VEC3 P1 = M2FacePoints_[M2Face(face2cut)];
//    Dart F1 = M2FaceBranch_[M2Face(face2cut)];

//    Dart cut0, cut1;
//    VEC3 AC = (Quadp[2] - Quadp[0]).normalized();
//    VEC3 BD = (Quadp[3] - Quadp[1]).normalized();
//    VEC3 P0P1 = (P1 - P0).normalized();
//    if(abs(AC.dot(P0P1)) < abs(BD.dot(P0P1))){
//        cut0 = Quadd[0]; cut1 = Quadd[2];
//    } else {
//        cut0 = Quadd[1]; cut1 = Quadd[3];
//    }
//    cmap2_->cut_face(cut0, cut1);
//    M2Vertex v = cmap2_->cut_edge(M2Edge(cmap2_->phi_1(cut0)));
//    M2Position_[v] = project_on_sphere((P0 + P1) * Scalar(0.5));

//    Dart newFace, oldFace;
//    VEC3 out0 = M2Position_[v] - center_;
//    VEC3 out1 = ((P0 - M2Position_[cut0]).normalized().cross((P1 - M2Position_[cut0]).normalized()));
//    if(out1.dot(out0) >= 0){
//        newFace = cut0; oldFace = cut1;
//    }
//    else{
//        newFace = cut1; oldFace = cut0;
//    }

//    M2FaceBranch_[M2Face(newFace)] = F0;
//    M2FaceBranch_[M2Face(oldFace)] = F1;
//    M2FacePoints_[M2Face(newFace)] = P0;
//    M2FacePoints_[M2Face(oldFace)] = P1;

//    move_points(newFace);
//}

//VEC3 Intersection_M2Builder::project_on_sphere(VEC3 P){
//    return center_ + (P - center_).normalized() * radius_;
//}

//VEC3 Intersection_M2Builder::mean_dir(M2Vertex vert){
//    using Quat = Eigen::Quaterniond;

//    uint32 valence = cmap2_->nb_darts_of_orbit(M2Vertex(vert));

//    std::vector<Quat> rotations;
//    rotations.reserve(valence);

//    VEC3 start = (M2Position_[vert] - center_).normalized();
//    cmap2_->foreach_dart_of_orbit(vert, [&](Dart face){
//        VEC3 dir = (M2FacePoints_[M2Face(face)] - center_).normalized();
//        Quat q = Quat::FromTwoVectors(start, dir);
//        q.normalize();
//        rotations.push_back(q);
//    });

//    Eigen::MatrixXd m(4, valence);
//    for(uint32 i = 0; i < valence; ++i){
//        const Quat& q = rotations[i];
//        m.col(i) = VEC4(q.w(), q.x(), q.y(), q.z());
//    }

//    Eigen::MatrixXd mm = m * m.transpose();
//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mm);
//    VEC4 r = es.eigenvectors().col(3);

//    Quat mean_rot(r[0], r[1], r[2], r[3]);
//    mean_rot.normalize();
//    return mean_rot._transformVector(start) * radius_ + center_;

////    VEC3 mean_dir = (M2Position_[vert] - center_).normalized();
////    for(uint32 i = 0; i < 50; ++i){
////    Eigen::Matrix3d rot;
////    rot.setIdentity();
////    cmap2_->foreach_dart_of_orbit(vert, [&](Dart face){
////        VEC3 v = (M2FacePoints_[M2Face(face)] - center_).normalized();
////        double angle = cgogn::geometry::angle(mean_dir, v);
////        angle = -1.0 * (M_PI - angle) / 100.0;
////        Eigen::AngleAxisd aa({ angle, mean_dir.cross(v) });
////        rot *= aa.toRotationMatrix();
////    });
////    mean_dir = rot * mean_dir;
////    mean_dir.normalize();
////    }
////    return -mean_dir * radius_ + center_;
//}

////VEC3 Intersection_M2Builder::mean_dir2(M2Vertex vert){

////}

//void Intersection_M2Builder::move_point(Dart vert){
////    uint32 valence = cmap2_->nb_darts_of_orbit(M2Vertex(vert));
////    if(valence == 2){
////        VEC3 V0 = M2FacePoints_[vert];
////        VEC3 V1 = M2FacePoints_[cmap2_->phi<2>(vert)];
////        M2Position_[vert] = project_on_sphere((V0 + V1) / Scalar(2));
////        return;
////    }
////    if(valence == 3){
////        VEC3 V0 = M2FacePoints_[vert];
////        VEC3 V1 = M2FacePoints_[cmap2_->phi<2>(vert)];
////        VEC3 V2 = M2FacePoints_[cmap2_->phi<212>(vert)];
////        VEC3 N = ((V2 - V0).normalized()).cross((V1 - V0).normalized()).normalized();
////        M2Position_[vert] = center_ + N * radius_;
////        return;
////    }
////    if(valence > 3){
//////        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic > coord(3, valence);
//////        uint32 i = 0;
//////        cmap2_->foreach_dart_of_orbit(M2Vertex(vert), [&](Dart face){
//////            coord.col(i++) = M2FacePoints_[face];
//////        });

//////        VEC3 centroid (coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());
//////        coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);
//////        auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
//////        VEC3 N = svd.matrixU().rightCols<1>();
//////        cgogn_log_info("result ") << N[0] << " " << N[1] << " " << N[2];
//////        if(N.dot(centroid - center_) < 0) N = -N;

//////        M2Position_[vert] = center_ + N * radius_;
////        for(uint32 i = 0; i < 10; ++i){
////        M2Position_[vert] = mean_dir(M2Vertex(vert));
////        }
////        return;
////    }
////    for(uint32 i = 0; i < 30; ++i){
//    M2Position_[vert] = mean_dir(M2Vertex(vert));
////    }
//}

//cgogn::Dart Intersection_M2Builder::convex_quad(Dart f){
//    Dart res;
//    VEC3 A = M2Position_[M2Vertex(f)];
//    VEC3 B = M2Position_[M2Vertex(cmap2_->phi<1>(f))];
//    VEC3 C = M2Position_[M2Vertex(cmap2_->phi<11>(f))];
//    VEC3 D = M2Position_[M2Vertex(cmap2_->phi<111>(f))];
//    VEC3 AB = B - A;
//    VEC3 AC = C - A;
//    VEC3 AD = D - A;
//    VEC3 N0 = AC.cross(AD);
//    VEC3 N1 = AB.cross(AC);
//    VEC3 N = N0.cross(N1).normalized();

//    if(N.dot(AC) < 0){
//        res = cmap2_->phi<1>(f);
//    }
//    else{
//        res = f;
//    }
//    return res;
//}

//void Intersection_M2Builder::move_points(Dart face){
//    cmap2_->foreach_dart_of_orbit(M2Face(face), [&](Dart vert){
//        move_point(vert);
//    });
//}

//void Intersection_M2Builder::move_points(){
//    cmap2_->foreach_cell([&](M2Face face){
////        convex_quad(face.dart);
//    });
//}

//void Intersection_M2Builder::clear(){

//}

} // namespace schnapps
