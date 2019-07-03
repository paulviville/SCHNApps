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
//        add_radius(ug, 1.0);
//        ug_stats(ug);

//        VEC3 v0 = VEC3(1,1,1);
//        VEC3 v1 = VEC3(2,2,2);
//        VEC3 v2 = VEC3(3,3,3);
//        Matrix3d m;
//        m << v0, v1, v2;
//        cgogn_log_info("m:") << m.col(0)[0] << " " << m.col(0)[1] << " " << m.col(0)[2];
//        cgogn_log_info("m:") << m.col(1)[0] << " " << m.col(1)[1] << " " << m.col(1)[2];
//        cgogn_log_info("m:") << m.col(2)[0] << " " << m.col(2)[1] << " " << m.col(2)[2];
//        m << v0, v1, v0;
//        cgogn_log_info("m:") << m.col(0)[0] << " " << m.col(0)[1] << " " << m.col(0)[2];
//        cgogn_log_info("m:") << m.col(1)[0] << " " << m.col(1)[1] << " " << m.col(1)[2];
//        cgogn_log_info("m:") << m.col(2)[0] << " " << m.col(2)[1] << " " << m.col(2)[2];


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

        /// building connection interfaces
        valid = build_connection_interfaces(ug, ug_attribs, cmap2);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to build connection interfaces";
            return false;
        }
        cgogn_log_info("success: ") << "connection interfaces built";

        CMap2_Attributes m2_attribs;
        m2_add_attributes(cmap2, m2_attribs);

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

        valid = set_edge_geometry(ug, ug_attribs, cmap2, m2_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to set edge geometry";
            return false;
        }
        cgogn_log_info("success: ") << "interface edge geometry built";

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

        m2_attribs_clean_up(cmap2, m2_attribs);
        ug_attribs_clean_up(ug, ug_attribs);
//        subdivide(cmap3);
        ug.foreach_cell([&](UGEdge uge) -> bool{
//           add_layer_edge(ug, ug_attribs, cmap3, uge);
//           add_layer_edge(ug, ug_attribs, cmap3, uge);
//           cgogn_log_info("integrity: ") << cmap3.check_map_integrity();
//           add_layer_edge(ug, ug_attribs, cmap3, UGEdge(ug.alpha1(uge.dart)));
//           add_layer_edge(ug, ug_attribs, cmap3, UGEdge(ug.alpha1(ug.alpha0(uge.dart))));
//           add_layer_edge(ug, ug_attribs, cmap3, UGEdge(ug.alpha1(ug.alpha0(uge.dart))));
//           add_layer_edge(ug, ug_attribs, cmap3, uge);
           return false;
//           return true;
        });

//        uint boundary_darts = 0;
//        uint not_boundary = 0;
//        cmap3.foreach_dart([&](Dart d) -> bool{
//            bool boundary = cmap3.is_boundary(d);
//            if(boundary) boundary_darts++;
//            else not_boundary++;
//            cgogn_log_info("d: ") << boundary << " " << d << " // phi1: " << cmap3.phi1(d) << " // phi2: " << cmap3.phi2(d)<< " // phi3: " << cmap3.phi3(d);
//            return true;
//        });
//        cgogn_log_info("boundary faces:") << (boundary_darts/ 4) << " " << boundary_darts;
//        cgogn_log_info("not boundary:") << not_boundary;

//        subdivide(cmap3);
        CMap3::VertexAttribute<VEC3> m3pos = cmap3.get_attribute<VEC3, M3Vertex>("position");
        CMap3_Quality_Attributes m3_QA_Attribs;
        quality_build_frames(cmap3, m3pos, m3_QA_Attribs);
        quality_scaled_jacobian(cmap3, m3_QA_Attribs);

        return true;
    }

    bool build_hexmesh2(UGraph& ug, CMap3& cmap3, CMap2& cmap2){
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
        valid = ug_clean(ug, ug_attribs);
        if(!valid){
            cgogn_log_info("ug_preprocessing error: ") << "clean up failed";
            return false;
        }

        /// collecting simplified graph information: branching points and branches
        Graph simplified_graph = ug_analysis(ug);
        if(!simplified_graph.valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to analyze";
            return false;
        }
        cgogn_log_info("success: ") << "graph analyzed";
        cgogn_log_info("") << simplified_graph.intersections.size() << " intersections; "
                                           << simplified_graph.branches.size() << " branches";

        /// building connection interfaces
        valid = build_intersection_interfaces(ug, ug_attribs, cmap2);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to build connection interfaces";
            return false;
        }
        cgogn_log_info("success: ") << "connection interfaces built";

        CMap2_Attributes m2_attribs;
        m2_add_attributes(cmap2, m2_attribs);

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

        valid = set_edge_geometry(ug, ug_attribs, cmap2, m2_attribs);
        if(!valid){
            cgogn_log_info("error build_hexmesh: ") << "failed to set edge geometry";
            return false;
        }
        cgogn_log_info("success: ") << "interface edge geometry built";

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

        m2_attribs_clean_up(cmap2, m2_attribs);
        ug_attribs_clean_up(ug, ug_attribs);
        return true;
    }


    void add_radius(UGraph& ug, Scalar radius){
        UGraph::VertexAttribute<Scalar> radii = ug.get_attribute<Scalar, UGVertex>("radii");
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
        ug_attribs.radii = ug.get_attribute<Scalar, UGVertex>("radii");
        ug_attribs.position = ug.get_attribute<VEC3, UGVertex>("position");
        if(!ug_attribs.radii.is_valid()){
            cgogn_log_info("ug_get_attribute error: ") << "radii attribute invalid";
            valid &= false;
        }
        if(!ug_attribs.position.is_valid()){
            cgogn_log_info("ug_get_attribute error: ") << "position attribute invalid";
            valid &= false;
        }
        ug_attribs.m2_CC = ug.add_attribute<Dart, UGVertex>("m2_CC");
        ug_attribs.m2_interface = ug.add_attribute<Dart, UGDart>("m2_interface");
        ug_attribs.frames = ug.add_attribute<Frame, UGDart>("frames");
        ug_attribs.sections = ug.add_attribute<Dart, UGDart>("sections");
//        ug_attribs.layers = ug.add_attribute<uint, UGEdge>("layers");
        return valid;
    }

    bool m2_add_attributes(CMap2& cmap2, CMap2_Attributes& m2_attribs){
        m2_attribs.position = cmap2.add_attribute<VEC3, M2Vertex>("position");
        m2_attribs.center = cmap2.add_attribute<VEC3, M2Volume>("center");
        m2_attribs.f_branch = cmap2.add_attribute<Dart, M2Face>("f_branches");
        m2_attribs.f_point = cmap2.add_attribute<VEC3, M2Face>("f_points");
        m2_attribs.edge_pos = cmap2.add_attribute<VEC3, M2Edge>("edge_position");
        m2_attribs.connections = cmap2.add_attribute<Dart, M2Dart>("connections");
        return true;
    }

    void m2_attribs_clean_up(CMap2& cmap2, CMap2_Attributes& m2_attribs){
        cmap2.remove_attribute(m2_attribs.center);
        cmap2.remove_attribute(m2_attribs.f_branch);
        cmap2.remove_attribute(m2_attribs.f_point);
        cmap2.remove_attribute(m2_attribs.connections);
    }

    void ug_attribs_clean_up(UGraph& ug, UG_Attributes& ug_attribs){
        ug.remove_attribute(ug_attribs.m2_interface);
        ug.remove_attribute(ug_attribs.frames);
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
        UGraph::VertexAttribute<VEC3>& ug_pos = ug_attribs.position;
        UGraph::VertexAttribute<Scalar>& ug_radii = ug_attribs.radii;

            ug.foreach_cell([&](UGEdge uge){

                VEC3 P0 = ug_pos[uge.dart];
                VEC3 Pn = ug_pos[ug.alpha0(uge.dart)];
                Scalar R0 = ug_radii[uge.dart];
                Scalar Rn = ug_radii[ug.alpha0(uge.dart)];

                Scalar avg_radius = (R0 + Rn)/2;
                VEC3 edge = Pn - P0;
                Scalar D = edge.norm();
                edge = edge.normalized();

                /// Placeholder for n
                uint n = D / avg_radius;
                if(n > 1){
                    Scalar y = (Rn - R0)/D;
                    Scalar ratio = pow(Rn/R0, 1.0/Scalar(n));

                    Scalar sum_ratio = 0;
                    Scalar alpha = 1;
                    std::vector<Scalar> ratios_sums;
                    for(uint i = 0; i < n; ++i){
                        sum_ratio += alpha;
                        ratios_sums.push_back(sum_ratio);
                        alpha *= ratio;
                    }
                    ratios_sums.pop_back(); // last is useless because P0 + sum_ratio * edge = Pn

                    VEC3 D0 = D / sum_ratio * edge;
                    std::vector<VEC3> Pi;
                    for(Scalar r : ratios_sums){
                        Pi.push_back(P0 + r * D0);
                    }

                    Dart d  = uge.dart;
                    alpha = ratio * R0;
                    for(VEC3 P : Pi){
                        d = ug.cut_edge(UGEdge(d)).dart;
                        ug_pos[d] = P;
                        ug_radii[d] = alpha;
                        alpha *= ratio;
                        d = ug.alpha1(d);
                    }
                }
            });
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
//        if(cmap2.is_embedded<M2Vertex>()){
//            m2builder.new_orbit_embedding(M2Vertex(Fd[0]));
//            m2builder.new_orbit_embedding(M2Vertex(cmap2.phi<1>(Fd[0])));
//            m2builder.new_orbit_embedding(M2Vertex(cmap2.phi<11>(Fd[0])));
//            m2builder.new_orbit_embedding(M2Vertex(cmap2.phi<111>(Fd[0])));
//        }

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
        UGraph::VertexAttribute<VEC3>& ug_pos = ug_attribs.position;
        CMap2::VertexAttribute<VEC3>& m2_pos = m2_attribs.position;

        VEC3 center = ug_pos[ugv];
        Scalar radius = ug_attribs.radii[ugv];
        Dart cc = ug_attribs.m2_CC[ugv];

        std::vector<Dart> Fd = {cc, cmap2.phi<2111>(cc), cmap2.phi<12>(cc)};
        std::vector<Dart> Pdart = {ugv.dart, ug.alpha1(ugv.dart), ug.alpha1(ug.alpha1(ugv.dart))};
        std::vector<VEC3> Ppos = {
            project_on_sphere(ug_pos[ug.alpha0(Pdart[0])], radius, center),
            project_on_sphere(ug_pos[ug.alpha0(Pdart[1])], radius, center),
            project_on_sphere(ug_pos[ug.alpha0(Pdart[2])], radius, center)
        };
        VEC3 V = (Ppos[1] - Ppos[0]).cross(Ppos[2] - Ppos[0]).normalized();
        std::vector<VEC3> Q = {center + V * radius, center - V * radius};
        std::vector<VEC3> M = {center + (Ppos[1] - Ppos[0]).normalized().cross(V) * radius,
                 center + (Ppos[2] - Ppos[1]).normalized().cross(V) * radius,
                 center + (Ppos[0] - Ppos[2]).normalized().cross(V) * radius};

//        M[0] = mean_dir(center, radius, M[0], {Ppos[0], Ppos[1]});
//        M[1] = mean_dir(center, radius, M[1], {Ppos[2], Ppos[1]});
//        M[2] = mean_dir(center, radius, M[2], {Ppos[0], Ppos[2]});

        m2_pos[Fd[0]] = M[0];
        m2_pos[Fd[1]] = M[1];
        m2_pos[Fd[2]] = M[2];
        m2_pos[cmap2.phi<1>(Fd[0])] = Q[0];
        m2_pos[cmap2.phi<111>(Fd[0])] = Q[1];

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
        ug.foreach_dart_of_orbit(ugv, [&](Dart d){
            Dart d0 = ug_attribs.m2_interface[d];
            Dart d1 = cmap2.phi<11>(d0);
            VEC3 R, S, T, diag, temp;
            T = (ug_attribs.position[ug.alpha0(d)] - center).normalized();
            diag =(m2_attribs.position[d1] - m2_attribs.position[d0]).normalized();
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
            Scalar cos = UE.col(0).dot(U_.col(0));
            Scalar angle = cos > 1 ? std::acos(1) : std::acos(cos);
            Scalar angle_step = angle / nb_e;

            cgogn_log_info("angle_step: ") << (UE.col(0).dot(U_.col(0)) > 1) << " " << angle << " " << angle_step << " " << nb_e;
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
//                VEC3 edge_pos = center + ((m2_attribs.position[d0] - center) + (m2_attribs.position[d1] - center).normalized() * radius);
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

            ug_attribs.sections[uge.dart] = D;
            ug_attribs.sections[ug.alpha0(uge.dart)] = cmap3.phi<23112>(D);
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

    void subdivide(CMap3& cmap3){
        CMap3::VertexAttribute<VEC3> m3pos = cmap3.get_attribute<VEC3, M3Vertex>("position");
        M3DartMarker vertexMarker(cmap3);
        M3DartMarker edgeMarker(cmap3);

        cmap3.foreach_cell([&](M3Edge m3e){
            Dart e0 = m3e.dart;
            Dart e1 = cmap3.phi2(e0);
            VEC3 pos = (m3pos[e0] + m3pos[e1]) / 2;

            M3Vertex v = cmap3.cut_edge(m3e);
            m3pos[v] = pos;
            vertexMarker.mark_orbit(v);
        });

        cgogn_log_info("integrity: ") << cmap3.check_map_integrity();

        cmap3.foreach_cell([&](M3Face m3f){
            Dart f0 = m3f.dart;
            if(!vertexMarker.is_marked(f0))
                f0 = cmap3.phi1(f0);

            Dart f1 = cmap3.phi<11>(f0);
            Dart f2 = cmap3.phi<11>(f1);
            Dart f3 = cmap3.phi<11>(f2);

            cmap3.cut_face(f0, f2);
            M3Vertex v = cmap3.cut_edge(M3Edge(cmap3.phi_1(f0)));
            m3pos[v] = (m3pos[f0] + m3pos[f1] + m3pos[f2] + m3pos[f3])/4;

            cmap3.cut_face(cmap3.phi_1(f0), f1);
            cmap3.cut_face(cmap3.phi_1(f2), f3);
            cmap3.foreach_incident_edge(v, [&](M3Edge m3e){
                edgeMarker.mark_orbit(m3e);
            });
        });

        cgogn_log_info("integrity: ") << cmap3.check_map_integrity();

        cmap3.foreach_cell([&](M3Volume m3w){
            Dart a = m3w.dart;
            while(!edgeMarker.is_marked(a) || !edgeMarker.is_marked(cmap3.phi1(a)))
                a = cmap3.phi1(a);

            std::vector<Dart> path_a;
            for(uint i = 0; i < 8; i++){
                path_a.push_back(a);
                a = cmap3.phi<121>(a);
            }

            std::vector<Dart> bs;
            bs.push_back(cmap3.phi1(path_a[0]));
            bs.push_back(cmap3.phi<21>(path_a[1]));

            for(Dart d : path_a){
                cgogn_log_info("d: ") << " " << d << " // phi1: " << cmap3.phi1(d) << " // phi2: " << cmap3.phi2(d)<< " // phi3: " << cmap3.phi3(d) << " " << cmap3.is_boundary(cmap3.phi3(d));
                cgogn_log_info("d: ") << " " << d << " " << cmap3.phi<121>(d) << " " << cmap3.phi<121121>(d) << " " << cmap3.phi<121121121>(d) << " " << cmap3.phi<121121121121>(d) << " " << cmap3.phi<121121121121121>(d) << " " << cmap3.phi<121121121121121121>(d);
            }
            cmap3.cut_volume(path_a);

            Dart f0 = cmap3.phi2(path_a[0]);
            Dart f1 = cmap3.phi2(path_a[2]);
            Dart f2 = cmap3.phi2(path_a[4]);
            Dart f3 = cmap3.phi2(path_a[6]);
            cmap3.cut_face(f0, f2);
            M3Vertex v = cmap3.cut_edge(M3Edge(cmap3.phi_1(f0)));
            m3pos[v] = (m3pos[f0] + m3pos[f1] + m3pos[f2] + m3pos[f3])/4;

            cmap3.cut_face(cmap3.phi_1(f0), f3);
            cmap3.cut_face(cmap3.phi_1(f2), f1);

            std::vector<Dart> cs;
            cs.push_back(cmap3.phi<1211>(bs[0]));
            cs.push_back(cmap3.phi<12112121>(bs[0]));
            cs.push_back(cmap3.phi<1211>(bs[1]));
            cs.push_back(cmap3.phi<12112121>(bs[1]));

            for(Dart b0 : bs){
                Dart b = b0;
                std::vector<Dart> path_b;
                for(uint i = 0; i < 6; i++){
                    path_b.push_back(b);
                    b = cmap3.phi<121>(b);
                }
                cmap3.cut_volume(path_b);
                cmap3.cut_face(cmap3.phi2(path_b[1]), cmap3.phi2(path_b[4]));
            }

            for(Dart c0 : cs){
                Dart c = c0;
                std::vector<Dart> path_c;
                for(uint i = 0; i < 4; i++){
                    path_c.push_back(c);
                    c = cmap3.phi<121>(c);
                }
                cmap3.cut_volume(path_c);
            }
        });

        cgogn_log_info("integrity: ") << cmap3.check_map_integrity();
    }

    void add_layer_edge_corner(UGraph& ug, UG_Attributes& ug_attribs, CMap3& cmap3, UGEdge uge){
        CMap3::VertexAttribute<VEC3> m3pos = cmap3.get_attribute<VEC3, M3Vertex>("position");
        Dart d0 = ug_attribs.sections[uge.dart];
//        Dart d1 = cmap3.phi<2112>(d0);

        Dart d = cmap3.phi<12>(d0);
        while(!cmap3.is_boundary(cmap3.phi<3>(d))){
            d = cmap3.phi<32112>(d);
            cgogn_log_info("d: ") << d;
        }
        cgogn_log_info("d on boundary: ") << d;

        M3Builder m3builder(cmap3);

        std::vector<std::pair<Dart, Dart>> cycle;
        Dart e = d;
        Dart h = m3builder.add_prism_topo_fp(4u);
        Dart first = cmap3.phi2(h);

        Dart last = cmap3.phi<112>(h);
        cycle.push_back(std::pair<Dart, Dart>(e, cmap3.phi<12>(h)));
        e = cmap3.phi<12321>(e);
        while(e != d){
            h = m3builder.add_prism_topo_fp(4u);
            m3builder.sew_volumes_fp(cmap3.phi<2>(h), last);
            last = cmap3.phi<112>(h);
            cycle.push_back(std::pair<Dart, Dart>(e, cmap3.phi<12>(h)));
            // corner
            if(cmap3.is_boundary(cmap3.phi<123>(e))){
                e = cmap3.phi<121>(e);
                h = m3builder.add_prism_topo_fp(4u);
                m3builder.sew_volumes_fp(cmap3.phi<2>(h), last);
                last = cmap3.phi<12>(h);
            }
            // not corner
            else
                e = cmap3.phi<12321>(e);
//            cgogn_log_info("cycle: ") << cycle.size();
        }
        m3builder.sew_volumes_fp(first, last);

        if(cmap3.is_boundary(cmap3.phi<3>(h)))
                cgogn_log_info("h on boundary") << 0;
        cgogn_log_info("h not boundary") << h << " " << cmap3.phi3(h);
        M3Volume m3w = M3Volume(m3builder.close_hole_topo(h));
        m3builder.boundary_mark(m3w);

        if(cmap3.is_boundary(cmap3.phi<3>(h)))
                cgogn_log_info("h on boundary") << 0;
        cgogn_log_info("h not boundary") << h << " " << cmap3.phi3(h);

        cmap3.foreach_incident_vertex(CMap3::ConnectedComponent(first), [&](M3Vertex m3v){
            m3builder.new_orbit_embedding(m3v);
            m3pos[m3v] = VEC3(1,1,1);
        });

        for(std::pair<Dart, Dart> vol : cycle){
            cgogn_log_info("vol:")<< 0;
            cmap3.sew_volumes(M3Face(vol.first), M3Face(vol.second));
        }

        if(cmap3.is_boundary(cmap3.phi<3>(d)))
                cgogn_log_info("d on boundary") << 0;

        cgogn_log_info("cycle: ") << cycle.size();
        return;
    }


    void add_layer_edge(UGraph& ug, UG_Attributes& ug_attribs, CMap3& cmap3, UGEdge uge){
        uint nb_bound_face_0 = 0;
        cmap3.foreach_dart([&](Dart d){
            if(cmap3.is_boundary(d))
                nb_bound_face_0++;
        });

        CMap3::VertexAttribute<VEC3> m3pos = cmap3.get_attribute<VEC3, M3Vertex>("position");
        Dart d0 = ug_attribs.sections[uge.dart];

        Dart d = cmap3.phi<12>(d0);
        uint layers = 1;
        while(!cmap3.is_boundary(cmap3.phi<3>(d))){
            d = cmap3.phi<32112>(d);
            layers++;
        }

        M3Builder m3builder(cmap3);

        std::vector<std::pair<Dart, Dart>> cycle;
        std::vector<std::pair<Dart, Dart>> vols;
        Dart e = d;
        Dart h = m3builder.add_prism_topo_fp(4u);
        Dart first = cmap3.phi2(h);

        Dart last = cmap3.phi<112>(h);
        cycle.push_back(std::pair<Dart, Dart>(e, cmap3.phi<12>(h)));
        vols.push_back(std::pair<Dart, Dart>(h, cmap3.phi<211>(h)));
        e = cmap3.phi<12321>(e);
        while(e != d){
            h = m3builder.add_prism_topo_fp(4u);
            m3builder.sew_volumes_fp(cmap3.phi<2>(h), last);
            last = cmap3.phi<112>(h);
            cycle.push_back(std::pair<Dart, Dart>(e, cmap3.phi<12>(h)));
            vols.push_back(std::pair<Dart, Dart>(h, cmap3.phi<211>(h)));
            // corner
            if(cmap3.is_boundary(cmap3.phi<123>(e))){
                e = cmap3.phi<121>(e);
            }
            else
                e = cmap3.phi<12321>(e);
        }
        m3builder.sew_volumes_fp(first, last);

        M3Volume m3w = M3Volume(m3builder.close_hole_topo(h));
        m3builder.boundary_mark(m3w);

        cmap3.foreach_incident_vertex(CMap3::ConnectedComponent(first), [&](M3Vertex m3v){
            m3builder.new_orbit_embedding(m3v);
            m3pos[m3v] = VEC3(1,1,1);
        });

        for(std::pair<Dart, Dart> vol : cycle){
            cmap3.sew_volumes(M3Face(vol.first), M3Face(vol.second));
        }

        Scalar radius0 = ug_attribs.radii[UGVertex(uge.dart)];
        Scalar radius1 = ug_attribs.radii[UGVertex(ug.alpha0(uge.dart))];
        VEC3 pos0 = ug_attribs.position[UGVertex(uge.dart)];
        VEC3 pos1 = ug_attribs.position[UGVertex(ug.alpha0(uge.dart))];

//        cgogn_log_info("new cycle") << 100000000;
        for(std::pair<Dart, Dart> vol : vols){
//            Dart h0, h1;
//            h0 = vol.first;
//            h1 = cmap3.phi<211>(vol.second);

//            bool existing_0 = false;
//            Dart a0 = cmap3.phi<12323>(h0);
//            if(!cmap3.is_boundary(a0)){
////                cgogn_log_info("not extremity") << 0;
//                Dart a1 = cmap3.phi<23>(a0);
//                if(!cmap3.is_boundary(a1)){
////                    cgogn_log_info("extra layer") << 0;
//                    Dart a2 = cmap3.phi<21>(a1);
//                    cmap3.sew_volumes(M3Face(a2), M3Face(h0));
//                    existing_0 = true;
//                }
//            }

//            bool existing_1 = false;
//            a0 = cmap3.phi<12323>(h1);
//            if(!cmap3.is_boundary(a0)){
////                cgogn_log_info("not extremity") << 0;
//                Dart a1 = cmap3.phi<23>(a0);
//                if(!cmap3.is_boundary(a1)){
////                    cgogn_log_info("extra layer") << 0;
//                    Dart a2 = cmap3.phi<21>(a1);
//                    cmap3.sew_volumes(M3Face(a2), M3Face(h1));
//                    existing_1 = true;
//                }
//            }

//            if(!existing_0){
                VEC3 dir0 = m3pos[cmap3.phi1(vol.first)] - pos0;
                m3pos[vol.first] = m3pos[cmap3.phi1(vol.first)] + dir0.normalized() * radius0;
//            }
//            if(!existing_1){
                VEC3 dir1 = m3pos[cmap3.phi1(vol.second)] - pos1;
                m3pos[vol.second] = m3pos[cmap3.phi1(vol.second)] + dir1.normalized() * radius1;
//            }

        }

        uint nb_bound_face_1 = 0;
        cmap3.foreach_dart([&](Dart d){
            if(cmap3.is_boundary(d))
                nb_bound_face_1++;
        });
        cgogn_log_info("nb bound faces: ") << (nb_bound_face_0);
        cgogn_log_info("nb bound faces: ") << (nb_bound_face_1);

        return;
    }

    void quality_build_frames(CMap3& cmap3, CMap3::VertexAttribute<VEC3>& m3pos, CMap3_Quality_Attributes& m3_QA_Attribs){
        m3_QA_Attribs.hexFrame = cmap3.get_attribute<Frame, M3Volume>("hexFrame");
        if(!m3_QA_Attribs.hexFrame.is_valid()){
            cgogn_log_info("added hexframe: ") << 1;
            m3_QA_Attribs.hexFrame = cmap3.add_attribute<Frame, M3Volume>("hexFrame");
        }
        m3_QA_Attribs.cornerFrame = cmap3.get_attribute<Frame, M3Vertex2>("cornerFrame");
        if(!m3_QA_Attribs.cornerFrame.is_valid()){
            m3_QA_Attribs.cornerFrame = cmap3.add_attribute<Frame, M3Vertex2>("cornerFrame");
        }
        m3_QA_Attribs.diagonals = cmap3.get_attribute<Frame, M3Volume>("diagonals");
        if(!m3_QA_Attribs.diagonals.is_valid()){
            m3_QA_Attribs.diagonals = cmap3.add_attribute<Frame, M3Volume>("diagonals");
        }

        cmap3.foreach_cell([&](M3Volume m3w) -> bool{
            quality_build_frames_hexa(cmap3, m3pos, m3_QA_Attribs, m3w);
            return true;
        });

    }

    void quality_build_frames_hexa(CMap3& cmap3, CMap3::VertexAttribute<VEC3>& m3pos, CMap3_Quality_Attributes& m3_QA_Attribs, M3Volume m3w){
        Dart d0 = m3w.dart;
        Dart D[8];
        D[0] = d0;
        D[1] = cmap3.phi<1>(d0);
        D[3] = cmap3.phi<211>(d0);
        D[2] = cmap3.phi<1>(D[3]);
        D[5] = cmap3.phi<1>(D[1]);
        D[4] = cmap3.phi<1>(D[5]);
        D[6] =cmap3.phi<211>(D[5]);
        D[7] = cmap3.phi<1>(D[6]);

        VEC3 P[8];
        P[0] = m3pos[D[0]];
        P[1] = m3pos[D[1]];
        P[2] = m3pos[D[2]];
        P[3] = m3pos[D[3]];
        P[4] = m3pos[D[4]];
        P[5] = m3pos[D[5]];
        P[6] = m3pos[D[6]];
        P[7] = m3pos[D[7]];

        m3_QA_Attribs.hexFrame[m3w] << (P[1] + P[2] + P[5] + P[6] - P[0] - P[3] - P[4] - P[7]),
                (P[3] + P[2] + P[7] + P[6] - P[0] - P[1] - P[4] - P[5]),
                (P[4] + P[5] + P[6] + P[7] - P[0] - P[1] - P[2] - P[3]);

        m3_QA_Attribs.cornerFrame[D[0]] << (P[1] - P[0]), (P[3] - P[0]), (P[4] - P[0]);
        m3_QA_Attribs.cornerFrame[D[1]] << (P[2] - P[1]), (P[0] - P[1]), (P[5] - P[1]);
        m3_QA_Attribs.cornerFrame[D[2]] << (P[3] - P[2]), (P[1] - P[2]), (P[6] - P[2]);
        m3_QA_Attribs.cornerFrame[D[3]] << (P[0] - P[3]), (P[2] - P[3]), (P[7] - P[3]);
        m3_QA_Attribs.cornerFrame[D[4]] << (P[0] - P[4]), (P[7] - P[4]), (P[5] - P[4]);
        m3_QA_Attribs.cornerFrame[D[5]] << (P[7] - P[5]), (P[6] - P[5]), (P[1] - P[5]);
        m3_QA_Attribs.cornerFrame[D[6]] << (P[7] - P[6]), (P[2] - P[6]), (P[5] - P[6]);
        m3_QA_Attribs.cornerFrame[D[7]] << (P[6] - P[7]), (P[4] - P[7]), (P[3] - P[7]);
    }

    void quality_scaled_jacobian(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs){
        m3_QA_Attribs.Scaled_Jacobian = cmap3.get_attribute<Scalar, M3Volume>("Scaled_Jacobian");
        if(!m3_QA_Attribs.Scaled_Jacobian.is_valid()){
            m3_QA_Attribs.Scaled_Jacobian = cmap3.add_attribute<Scalar, M3Volume>("Scaled_Jacobian");
        }

        cmap3.foreach_cell([&](M3Volume m3w) -> bool{
            quality_scaled_jacobian_hexa(cmap3, m3_QA_Attribs, m3w);
            return true;
        });
    }

    void quality_scaled_jacobian_hexa(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs, M3Volume m3w){
        Frame m = m3_QA_Attribs.hexFrame[m3w];
        m.col(0).normalize();
        m.col(1).normalize();
        m.col(2).normalize();

        Scalar jacobian = m.determinant();
        cmap3.foreach_incident_vertex(m3w, [&](M3Vertex m3v){
            m = m3_QA_Attribs.cornerFrame[m3v.dart];
            m.col(0).normalize();
            m.col(1).normalize();
            m.col(2).normalize();
            Scalar temp = m.determinant();
            if(temp < jacobian)
                jacobian = temp;
        });
        m3_QA_Attribs.Scaled_Jacobian[m3w] = jacobian;
    }

    void quality_jacobian(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs){
        m3_QA_Attribs.Jacobian = cmap3.get_attribute<Scalar, M3Volume>("Jacobian");
        if(!m3_QA_Attribs.Jacobian.is_valid()){
            m3_QA_Attribs.Jacobian = cmap3.add_attribute<Scalar, M3Volume>("Jacobian");
        }

        cmap3.foreach_cell([&](M3Volume m3w) -> bool{
            quality_jacobian_hexa(cmap3, m3_QA_Attribs, m3w);
            return true;
        });
    }

    void quality_jacobian_hexa(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs, M3Volume m3w){
        Scalar jacobian = m3_QA_Attribs.hexFrame[m3w].determinant();
        cmap3.foreach_incident_vertex(m3w, [&](M3Vertex m3v){
            Scalar temp = m3_QA_Attribs.cornerFrame[m3v.dart].determinant();
            if(temp < jacobian)
                jacobian = temp;
        });
        m3_QA_Attribs.Jacobian[m3w] = jacobian;
    }

    void quality_mean_frobenius(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs){
        m3_QA_Attribs.Mean_Frobenius = cmap3.get_attribute<Scalar, M3Volume>("Mean_Frobenius");
        if(!m3_QA_Attribs.Mean_Frobenius.is_valid()){
            m3_QA_Attribs.Mean_Frobenius = cmap3.add_attribute<Scalar, M3Volume>("Mean_Frobenius");
        }

        cmap3.foreach_cell([&](M3Volume m3w) -> bool{
            quality_mean_frobenius_hexa(cmap3, m3_QA_Attribs, m3w);
            return true;
        });
    }

    void quality_max_frobenius(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs){
        m3_QA_Attribs.Max_Forbenius = cmap3.get_attribute<Scalar, M3Volume>("Max_Forbenius");
        if(!m3_QA_Attribs.Max_Forbenius.is_valid()){
            m3_QA_Attribs.Max_Forbenius = cmap3.add_attribute<Scalar, M3Volume>("Max_Forbenius");
        }

        cmap3.foreach_cell([&](M3Volume m3w) -> bool{
            quality_max_frobenius_hexa(cmap3, m3_QA_Attribs, m3w);
            return true;
        });
    }

    void quality_mean_frobenius_hexa(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs, M3Volume m3w){
        Scalar jacobian = m3_QA_Attribs.hexFrame[m3w].determinant();
        cmap3.foreach_incident_vertex(m3w, [&](M3Vertex m3v){
            Scalar temp = m3_QA_Attribs.cornerFrame[m3v.dart].determinant();
            if(temp < jacobian)
                jacobian = temp;
        });
        m3_QA_Attribs.Jacobian[m3w] = jacobian;
    }

    void quality_max_frobenius_hexa(CMap3& cmap3, CMap3_Quality_Attributes& m3_QA_Attribs, M3Volume m3w){
        Scalar jacobian = m3_QA_Attribs.hexFrame[m3w].determinant();
        cmap3.foreach_incident_vertex(m3w, [&](M3Vertex m3v){
            Scalar temp = m3_QA_Attribs.cornerFrame[m3v.dart].determinant();
            if(temp < jacobian)
                jacobian = temp;
        });
        m3_QA_Attribs.Jacobian[m3w] = jacobian;
    }
}
}
