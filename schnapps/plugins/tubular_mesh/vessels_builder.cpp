#include <schnapps/plugins/tubular_mesh/vessels_builder.h>

#include <cgogn/geometry/types/eigen.h>
#include <cgogn/geometry/functions/orientation.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/geometry/functions/basics.h>
namespace schnapps
{

Vessels_Builder::Vessels_Builder() {}

void Vessels_Builder::set_skeleton(UGraph* ug){
    ug_ = ug;
    ug_marker_ = new UGDartMarker(*ug_);
    UGNormals_ = ug_->add_attribute<VEC3, UGVertex>("normals");
    UGTangents_ = ug_->add_attribute<VEC3, UGVertex>("tangents");
    UGposition_ = ug_->template get_attribute<VEC3, UGVertex>("position");

}

void Vessels_Builder::set_cmap3(CMap3* cmap3){
    cmap3_ = cmap3;
    m3builder_ = new M3Builder(*cmap3_);
    if (! cmap3_->is_embedded<M3Vertex>())
            m3builder_->create_embedding<M3Vertex::ORBIT>();
    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
    M3position_ = ca.get_chunk_array<VEC3>("position");
    if (M3position_ == nullptr)
        M3position_ = ca.add_chunk_array<VEC3>("position");
}

bool Vessels_Builder::compute_cmap3(){
    subdivide_graph();
    UGConnections_ = ug_->template add_attribute<Dart, UndirectedGraph::CDart>("connections");
    analyse_graph();
    build_cmap3();
//    clean_up();
    return true;
}

void Vessels_Builder::analyse_graph(){
    find_ends();
    find_branches();
    cgogn_log_info("graph analysis:\n") << "nb inter: " << intersections_.size() << "   nb branches: " << branches_.size();
    uint32 nb_branches = branches_.size();
    check_isolated_loops();
    cgogn_log_info("") << "\t\tnb loops: " << branches_.size() - nb_branches;

}

void Vessels_Builder::find_ends(){
    ug_->foreach_cell([&](UGVertex v){
        switch(ug_->nb_darts_of_orbit(v)){
            case 1: // end
                extremities_.push_back(v);
            case 2: // joint
            break;
            case 3: // intersection 3+
            default:
                extremities_.push_back(v);
                intersections_.push_back(v);
            break;
        }
    });
    cgogn_log_info("extremities: ") << extremities_.size();
}

void Vessels_Builder::find_branches(){
    for(UGVertex v : extremities_){
        Dart d0 = v.dart;
        if(!ug_marker_->is_marked(d0)){
            cgogn_log_info("extremities: ") << d0 << " / " << ug_->alpha1(d0);
            ug_marker_->mark(d0);
            Dart d1 = ug_->alpha0(d0);
            ug_marker_->mark(d1);
            while(!ug_->is_boundary(d1) && ug_->nb_darts_of_orbit(UGVertex(d1)) == 2){
                d1 = ug_->alpha1(d1);
                ug_marker_->mark(d1);
                d1 = ug_->alpha0(d1);
                ug_marker_->mark(d1);
            }
            if(!ug_->is_boundary(d1)) // pas un sommet isolé
                branches_.push_back(std::pair<Dart, Dart>(d0, d1));
        }
    }
}

void Vessels_Builder::check_isolated_loops(){
    ug_->foreach_dart([&](Dart d0){
        if(!ug_marker_->is_marked(d0)){
            ug_marker_->mark(d0);
            Dart d1 = ug_->alpha0(d0);
            ug_marker_->mark(d1);
            while(!(ug_->alpha1(d1) == d0) && !ug_->is_boundary(d1) && ug_->nb_darts_of_orbit(UGVertex(d1)) == 2){
                d1 = ug_->alpha1(d1);
                ug_marker_->mark(d1);
                d1 = ug_->alpha0(d1);
                ug_marker_->mark(d1);
            }
            if(!ug_->is_boundary(d1))
                branches_.push_back(std::pair<Dart, Dart>(d0, d1));
        }
    });
}

void Vessels_Builder::subdivide_graph(){
    Scalar radius(0.5);
    ug_marker_->unmark_all();
    bool over = false;
    while(!over){
        over = true;
        ug_->foreach_cell([&](UGEdge e){
            if(!ug_marker_->is_marked(e.dart)){
                Dart v0 = e.dart;
                Dart v1 = ug_->alpha0(e.dart);
                VEC3 edge_vec = UGposition_[v0] - UGposition_[v1];
                cgogn_log_info("edge: ") << edge_vec.norm();
                if(edge_vec.norm() > 2 * radius){
                    over &= false;
                    UGVertex vert = ug_->add_vertex();
                    UGposition_[vert] = (UGposition_[v0] + UGposition_[v1]) * Scalar(0.5);
//                    UGposition_[vert] += VEC3::Random() * 0.1;
                    cgogn_log_info("length longer than radius: ") << 1;
                    ug_->connect_vertices(UGVertex(v0), vert);
                    ug_->connect_vertices(UGVertex(v1), vert);
                    ug_->disconnect_vertices(e);
//                    ug_marker_->mark_orbit(e);
                }
                else{
                    ug_marker_->mark_orbit(e);
                }
            }
        });
    }
    ug_marker_->unmark_all();
}

void Vessels_Builder::build_cmap3(){
    cgogn_log_info("Building cmap3");
    build_intersections();
    compute_tangents();
    build_branches();
}

void Vessels_Builder::build_intersections(){
    intersection_m2builder_.set_cmap2(cmap2_);
    intersection_m2builder_.set_ugraph(ug_);

    for(UGVertex ugv : intersections_){
        build_intersection3(ugv);
//                build_inter3(ugv);
//        intersection_m2builder_.build_intersection(ugv);
    }
}

void Vessels_Builder::build_intersection3(UGVertex ugv){
    Scalar radius = Scalar(0.5f);
    VEC3 Ctr = UGposition_[ugv];
    std::vector<VEC3> P, PV, C, Q, M, E, F;
    std::vector<Dart> PD;

    ug_->foreach_dart([&](Dart d){
        if(ug_->is_boundary(d))
                cgogn_log_info("boundary: ") << ug_->is_boundary(d);
        if(ug_->is_boundary(ug_->alpha0(d)))
                cgogn_log_info("boundary: ") << ug_->is_boundary(d);
        cgogn_log_info("emb:") << d << " " <<ug_->embedding(d, UndirectedGraph::CDart::ORBIT);
    });

    ug_->foreach_adjacent_vertex_through_edge(ugv, [&](UGVertex ugv2){
        VEC3 v = (UGposition_[ugv2] - Ctr).normalized();
        cgogn_log_info("d") << ugv2 ;
        P.push_back(Ctr + v * Scalar(radius));
        PD.push_back(ug_->alpha0(ugv2.dart));
    });

    VEC3 V = (P[1] - P[0]).cross(P[2] - P[0]).normalized();
    Q.push_back(Ctr + V * Scalar(radius));
    Q.push_back(Ctr - V * Scalar(radius));

    PV.push_back((P[1]-P[0]).normalized());
    PV.push_back((P[2]-P[1]).normalized());
    PV.push_back((P[0]-P[2]).normalized());
    M.push_back(Ctr + PV[0].cross(V) * radius);
    M.push_back(Ctr + PV[1].cross(V) * radius);
    M.push_back(Ctr + PV[2].cross(V) * radius);

    C.push_back((M[0]+M[1]+M[2]+Q[0])*Scalar(0.25f));
    C.push_back((M[0]+M[1]+M[2]+Q[1])*Scalar(0.25f));

    E.push_back((M[0] + M[1])*Scalar(0.5f));
    E.push_back((M[1] + M[2])*Scalar(0.5f));
    E.push_back((M[2] + M[0])*Scalar(0.5f));
    E.push_back((M[0] + Q[0])*Scalar(0.5f));
    E.push_back((M[1] + Q[0])*Scalar(0.5f));
    E.push_back((M[2] + Q[0])*Scalar(0.5f));
    E.push_back((M[0] + Q[1])*Scalar(0.5f));
    E.push_back((M[1] + Q[1])*Scalar(0.5f));
    E.push_back((M[2] + Q[1])*Scalar(0.5f));

    F.push_back((M[0] + M[1] + Q[0]) / Scalar(3.0f));
    F.push_back((M[1] + M[2] + Q[0]) / Scalar(3.0f));
    F.push_back((M[2] + M[0] + Q[0]) / Scalar(3.0f));
    F.push_back((M[0] + M[1] + Q[1]) / Scalar(3.0f));
    F.push_back((M[1] + M[2] + Q[1]) / Scalar(3.0f));
    F.push_back((M[2] + M[0] + Q[1]) / Scalar(3.0f));

    /// Projection sur la shère de centre Ctr de l'embranchement
    for(uint32 i = 0; i < E.size(); i++){
        E[i] = (E[i] - Ctr).normalized();
        E[i] = Ctr + E[i] * radius;
    }

    for(uint32 i = 0; i < F.size(); i++){
        F[i] = (F[i] - Ctr).normalized();
        F[i] = Ctr + F[i] * radius;
    }

    /// building CMAP3
    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
    /// Insertion des points dans la 3 carte
    std::vector<uint32> Ci, Qi, Mi, Ei, Fi;
    uint32 Ctri;

    uint32 l = ca.insert_lines<1>();
    (*M3position_)[l] = Ctr;
    Ctri = l;

    for(VEC3 v: C){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Ci.push_back(l);
    }
    C.clear();

    for(VEC3 v: Q){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Qi.push_back(l);
    }
    Q.clear();

    for(VEC3 v: M){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Mi.push_back(l);
    }
    M.clear();

    for(VEC3 v: E){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Ei.push_back(l);
    }
    E.clear();

    for(VEC3 v: F){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Fi.push_back(l);
    }
    F.clear();

    /// Construction des hexaèdres de l'embranchement
    std::vector<Dart> H;
    H.push_back(add_hexa({Mi[0], Ei[2], Fi[2], Ei[3], Ei[0], Ctri, Ci[0], Fi[0]}));
    H.push_back(add_hexa({Ei[0], Ctri, Ci[0], Fi[0], Mi[1], Ei[1], Fi[1], Ei[4]}));
    H.push_back(add_hexa({Ctri, Ei[2], Fi[2], Ci[0], Ei[1], Mi[2], Ei[5], Fi[1]}));
    H.push_back(add_hexa({Ei[3], Fi[2], Ei[5], Qi[0], Fi[0], Ci[0], Fi[1], Ei[4]}));
    H.push_back(add_hexa({Ei[2], Mi[0], Ei[6], Fi[5], Ctri, Ei[0], Fi[3], Ci[1]}));
    H.push_back(add_hexa({Ei[0], Ctri, Ei[1], Mi[1], Fi[3], Ci[1], Fi[4], Ei[7]}));
    H.push_back(add_hexa({Ctri, Ei[2], Mi[2], Ei[1], Ci[1], Fi[5], Ei[8], Fi[4]}));
    H.push_back(add_hexa({Fi[5], Ei[6], Qi[1], Ei[8], Ci[1], Fi[3], Ei[7], Fi[4]}));

    /// couture topologique des hexaèdres de l'embranchement
    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[0]), H[1]);
    m3builder_->sew_volumes_fp(cmap3_->phi<121>(H[0]), H[2]);
    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[0]), cmap3_->phi<2>(H[3]));
    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[0]), cmap3_->phi<2>(H[4]));
    m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[1]), cmap3_->phi<1112>(H[2]));
    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[1]), cmap3_->phi<2112>(H[3]));
    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[1]), H[5]);
    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[2]), cmap3_->phi<121>(H[3]));
    m3builder_->sew_volumes_fp(cmap3_->phi<2>(H[2]), H[6]);
    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[4]), cmap3_->phi<2>(H[5]));
    m3builder_->sew_volumes_fp(cmap3_->phi<212>(H[4]), cmap3_->phi<2>(H[6]));
    m3builder_->sew_volumes_fp(cmap3_->phi<112>(H[4]), cmap3_->phi<2>(H[7]));
    m3builder_->sew_volumes_fp(cmap3_->phi<12>(H[5]), cmap3_->phi<1112>(H[6]));
    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[5]), cmap3_->phi<2112>(H[7]));
    m3builder_->sew_volumes_fp(cmap3_->phi<2112>(H[6]), cmap3_->phi<212>(H[7]));

    /// Mise en place de connections à l'embranchement dans l'attribut
    UGConnections_[PD[1]] = cmap3_->phi<2121>(H[0]);
    UGConnections_[PD[2]] = cmap3_->phi<21121>(H[1]);
    UGConnections_[PD[0]] = cmap3_->phi<123211>(H[0]);

    /// stockage de la normal
    UGNormals_[ugv] = V;
}

void Vessels_Builder::build_branches(){
    ug_marker_->unmark_all();
    for(std::pair<Dart, Dart> branch : branches_){
        if(ug_->nb_darts_of_orbit(UGVertex(branch.first)) >= 3){
            build_branch(branch.first);
        }
        else{
            build_branch(branch.second);
        }
    }
}

void Vessels_Builder::compute_tangents(){
    ug_marker_->unmark_all();
    for(std::pair<Dart, Dart> branch : branches_){
        if(ug_->nb_darts_of_orbit(UGVertex(branch.first)) >= 3){
            compute_tangent_3_1(branch.first);
        }
        else{
            compute_tangent_3_1(branch.second);
        }
    }
}

void Vessels_Builder::compute_tangent_3_1(Dart d){
    cgogn_log_info("branch 3 1");
    VEC3 N = UGNormals_[UGVertex(d)];
    Dart d0 = d;
    Dart d1 = ug_->alpha0(d);
    VEC3 V0 = (UGposition_[d1] - UGposition_[d0]).normalized();
    uint32 valence = ug_->nb_darts_of_orbit(UGVertex(d1));
    while(valence == 2){
        cgogn_log_info("darts: ") << d0 << " " << d1;
        d0 = ug_->alpha1(d1);
        d1 = ug_->alpha0(d0);
        VEC3 V1 = (UGposition_[d1] - UGposition_[d0]).normalized();
        UGNormals_[d0] = N;
        VEC3 T = (V1 + V0).normalized();
        UGTangents_[d0] = T;
        V0 = V1;
        valence = ug_->nb_darts_of_orbit(UGVertex(d1));
    }

    if(valence == 1){
        UGNormals_[d1] = N;
        UGTangents_[d1] = V0;
    }
}

void Vessels_Builder::build_branch(Dart d){
    cgogn_log_info("looping1");
    Scalar radius = Scalar(0.5f);
    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();

    Dart d0 = d;
    Dart d1 = ug_->alpha0(d0);
    uint32 valence = ug_->nb_darts_of_orbit(UGVertex(d1));

    while(true){
        cgogn_log_info("valence: ") << valence << " " << d0 << " " << d1 ;
        VEC3 T = UGTangents_[d1];
        VEC3 N = UGNormals_[d1];
        VEC3 M = (N.cross(T)).normalized();

        Dart D = UGConnections_[d0];
        if(D.is_nil()){
            cgogn_log_info("Branch error: connection is nil") << d0 << " " << D;
            return;
        }

        // entrance connection darts and
        std::vector<Dart> D0 = {D, cmap3_->phi<1232>(D), cmap3_->phi<11232>(D), cmap3_->phi<111232>(D), cmap3_->phi<1112321232>(D), cmap3_->phi<11123211232>(D)};
        std::vector<uint32> E0 = {
            cmap3_->embedding(D0[2], M3Vertex::ORBIT),
            cmap3_->embedding(D0[1], M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<1>(D0[0]), M3Vertex::ORBIT),
            cmap3_->embedding(D0[0], M3Vertex::ORBIT),
            cmap3_->embedding(D0[5], M3Vertex::ORBIT),
            cmap3_->embedding(D0[4], M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<111>(D0[5]), M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<11>(D0[5]), M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<11>(D0[4]), M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<11>(D0[2]), M3Vertex::ORBIT),
            cmap3_->embedding(cmap3_->phi<11>(D0[1]), M3Vertex::ORBIT)
        };

        std::vector<Dart> D1 = {add_hexa(), add_hexa(), add_hexa(), add_hexa(), add_hexa(), add_hexa()};
        m3builder_->sew_volumes_fp(cmap3_->phi<1112>(D1[0]), cmap3_->phi<2>(D1[1]));
        m3builder_->sew_volumes_fp(cmap3_->phi<112>(D1[0]), cmap3_->phi<2>(D1[2]));
        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[0]), cmap3_->phi<2>(D1[3]));
        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[1]), cmap3_->phi<1112>(D1[2]));
        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[2]), cmap3_->phi<1112>(D1[4]));
        m3builder_->sew_volumes_fp(cmap3_->phi<1112>(D1[3]), cmap3_->phi<2>(D1[4]));
        m3builder_->sew_volumes_fp(cmap3_->phi<112>(D1[3]), cmap3_->phi<2>(D1[5]));
        m3builder_->sew_volumes_fp(cmap3_->phi<12>(D1[4]), cmap3_->phi<1112>(D1[5]));

        std::vector<uint32> E1;
        /// gestion des extrémités et articulations
        if(valence <= 2){
            E1 = {ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(),
                  ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>(), ca.insert_lines<1>()
           };
            (*M3position_)[E1[0]] = UGposition_[d1];
            (*M3position_)[E1[1]] = UGposition_[d1] + N * radius * Scalar(0.5f);
            (*M3position_)[E1[2]] = UGposition_[d1] + (N - M).normalized() * radius;
            (*M3position_)[E1[3]] = UGposition_[d1] - M * radius;
            (*M3position_)[E1[4]] = UGposition_[d1] - (N + M).normalized() * radius;
            (*M3position_)[E1[5]] = UGposition_[d1] - N * radius * Scalar(0.5f);
            (*M3position_)[E1[6]] = UGposition_[d1] - N * radius;
            (*M3position_)[E1[7]] = UGposition_[d1] - (N - M).normalized() * radius;
            (*M3position_)[E1[8]] = UGposition_[d1] + M * radius;
            (*M3position_)[E1[9]] = UGposition_[d1] + (N + M).normalized() * radius;
            (*M3position_)[E1[10]] = UGposition_[d1] + N * radius;
        }
        /// gestion des embranchements
        if(valence >= 3){
            Dart E = UGConnections_[d1];
            cgogn_log_info("Branching point") << E ;
            if(E.is_nil()){
                cgogn_log_info("Branch error: connection is nil") << E ;
            }

            /// Inversion des connections si les normales sont trop différentes
            if(UGNormals_[UGVertex(d0)].dot(UGNormals_[UGVertex(d1)]) >= 0)
                E = cmap3_->phi<1123211>(E);
            else
                E = cmap3_->phi<111232111>(E);

            std::vector<Dart> D2 = {E, cmap3_->phi<111232>(E), cmap3_->phi<11232>(E), cmap3_->phi<1232>(E), cmap3_->phi<1232111232>(E), cmap3_->phi<123211232>(E)};
            E1 = {
                cmap3_->embedding(D2[3], M3Vertex::ORBIT),
                cmap3_->embedding(D2[2], M3Vertex::ORBIT),
                cmap3_->embedding(D2[0], M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<1>(D2[0]), M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<2>(D2[5]), M3Vertex::ORBIT),
                cmap3_->embedding(D2[5], M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<11>(D2[5]), M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<12>(D2[4]), M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<11>(D2[2]), M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<11>(D2[1]), M3Vertex::ORBIT),
                cmap3_->embedding(cmap3_->phi<111>(D2[1]), M3Vertex::ORBIT),
            };
        }

        embed_hexa(D1[0], {E0[2], E0[3], E0[0], E0[1], E1[2], E1[3], E1[0], E1[1]});
        embed_hexa(D1[1], {E0[2], E0[1], E0[9], E0[10], E1[2], E1[1], E1[9], E1[10]});
        embed_hexa(D1[2], {E0[1], E0[0], E0[8], E0[9], E1[1], E1[0], E1[8], E1[9]});
        embed_hexa(D1[3], {E0[0], E0[3], E0[4], E0[5], E1[0], E1[3], E1[4], E1[5]});
        embed_hexa(D1[4], {E0[0], E0[5], E0[7], E0[8], E1[0], E1[5], E1[7], E1[8]});
        embed_hexa(D1[5], {E0[5], E0[4], E0[6], E0[7], E1[5], E1[4], E1[6], E1[7]});

        if(valence == 1) break;
        if(valence == 2){
            d0 = ug_->alpha1(d1);
            d1 = ug_->alpha0(d0);
        }
        if(valence == 3) break;


        UGConnections_[d0] = cmap3_->phi<2112>(D1[0]);
        valence = ug_->nb_darts_of_orbit(UGVertex(d1));
    }

}

cgogn::Dart Vessels_Builder::add_hexa(){
    return m3builder_->add_prism_topo_fp(4u);
}

cgogn::Dart Vessels_Builder::add_hexa(std::vector<uint32> verts){
    Dart d = m3builder_->add_prism_topo_fp(4u);
    const std::array<Dart, 8> vertices_of_hexa = {
        d,
        cmap3_->phi1(d),
        cmap3_->phi1(cmap3_->phi1(d)),
        cmap3_->phi_1(d),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi_1(d))))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(d)))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(d))))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(d))))))
    };

    uint32 index = 0;
    for (Dart dv : vertices_of_hexa)
        m3builder_->template set_orbit_embedding<M3Vertex>(M3Vertex2(dv), verts[index++]);
    return d;
}

void Vessels_Builder::embed_hexa(Dart d, std::vector<uint32> verts){
    const std::array<Dart, 8> vertices_of_hexa = {
        d,
        cmap3_->phi1(d),
        cmap3_->phi1(cmap3_->phi1(d)),
        cmap3_->phi_1(d),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi_1(d))))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(d)))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(d))))),
        cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(cmap3_->phi2(cmap3_->phi1(cmap3_->phi1(d))))))
    };

    uint32 index = 0;
    for (Dart dv : vertices_of_hexa)
        m3builder_->template set_orbit_embedding<M3Vertex>(M3Vertex2(dv), verts[index++]);
}

void Vessels_Builder::clean_up(){
    ug_->remove_attribute(UGConnections_);
    ug_->remove_attribute(UGTangents_);
    ug_->remove_attribute(UGNormals_);
}


Intersection_M2Builder::Intersection_M2Builder() : radius_(1.0), cmap2_(nullptr), ug_(nullptr){

}

void Intersection_M2Builder::set_cmap2(CMap2 *cmap2){
    cmap2_ = cmap2;
    cmap2_->clear_and_remove_attributes();

    m2builder_ = new M2Builder(*cmap2_);
}

void Intersection_M2Builder::set_ugraph(UGraph *ug){
    ug_ = ug;
    UGposition_ = ug_->template get_attribute<VEC3, UGVertex>("position");
}

void Intersection_M2Builder::build_intersection(UGVertex ugv){
    cmap2_->clear_and_remove_attributes();

    Ppos_.clear();
    Pdart_.clear();

    ugv_ = ugv;
    center_ = UGposition_[ugv];
    radius_ = Scalar(0.5); /// PLACEHOLDER

    get_intersection_data();
    build_all();
}

bool Intersection_M2Builder::in_quad(Dart face, VEC3 P){
    using cgogn::geometry::intersection_ray_triangle;

    VEC3 A = M2Position_[face];
    VEC3 B = M2Position_[cmap2_->phi<1>(face)];
    VEC3 C = M2Position_[cmap2_->phi<11>(face)];
    VEC3 D = M2Position_[cmap2_->phi<111>(face)];
    VEC3 DIR = P - center_;

    return (intersection_ray_triangle(center_, DIR, A, B, C)
            || intersection_ray_triangle(center_, DIR, A, C, D));
}

void Intersection_M2Builder::get_intersection_data(){
    ug_->foreach_adjacent_vertex_through_edge(ugv_, [&](UGVertex ugv){
        Ppos_.push_back(project_on_sphere(UGposition_[ugv]));
        Pdart_.push_back(ug_->alpha0(ugv.dart));
                UGVertex v = ug_->add_vertex();
                UGposition_[v] = project_on_sphere(UGposition_[ugv]);
    });
}

void Intersection_M2Builder::build_core(){
    std::vector<VEC3> Qp, Mp;
    std::vector<Dart> Fd;

    Fd = {m2builder_->add_face_topo_fp(4), m2builder_->add_face_topo_fp(4), m2builder_->add_face_topo_fp(4)};

    m2builder_->phi2_sew(Fd[0], cmap2_->phi<1>(Fd[2]));
    m2builder_->phi2_sew(cmap2_->phi<1>(Fd[0]), Fd[1]);
    m2builder_->phi2_sew(cmap2_->phi<11>(Fd[0]), cmap2_->phi_1(Fd[1]));
    m2builder_->phi2_sew(cmap2_->phi_1(Fd[0]), cmap2_->phi<11>(Fd[2]));
    m2builder_->phi2_sew(cmap2_->phi<1>(Fd[1]), Fd[2]);
    m2builder_->phi2_sew(cmap2_->phi<11>(Fd[1]), cmap2_->phi_1(Fd[2]));

    M2Position_ = cmap2_->add_attribute<VEC3, M2Vertex>("position");
    M2FacePoints_ = cmap2_->add_attribute<VEC3, M2Face>("face_points");
    M2FaceBranch_ = cmap2_->add_attribute<Dart, M2Face>("face_branches");

    VEC3 V = (Ppos_[1] - Ppos_[0]).cross(Ppos_[2] - Ppos_[0]).normalized();
    Qp = {center_ + V * radius_, center_ - V * radius_};
    Mp = {center_ + (Ppos_[1]-Ppos_[0]).normalized().cross(V) * radius_,
         center_ + (Ppos_[2]-Ppos_[1]).normalized().cross(V) * radius_,
         center_ + (Ppos_[0]-Ppos_[2]).normalized().cross(V) * radius_};

    M2Position_[M2Vertex(Fd[0])] = Mp[2];
    M2Position_[M2Vertex(Fd[1])] = Mp[0];
    M2Position_[M2Vertex(Fd[2])] = Mp[1];
    M2Position_[M2Vertex(cmap2_->phi<1>(Fd[0]))] = Qp[1];
    M2Position_[M2Vertex(cmap2_->phi_1(Fd[0]))] = Qp[0];

    M2FacePoints_[M2Face(Fd[0])] = Ppos_[0];
    M2FacePoints_[M2Face(Fd[1])] = Ppos_[1];
    M2FacePoints_[M2Face(Fd[2])] = Ppos_[2];

    M2FaceBranch_[M2Face(Fd[0])] = Pdart_[0];
    M2FaceBranch_[M2Face(Fd[1])] = Pdart_[1];
    M2FaceBranch_[M2Face(Fd[2])] = Pdart_[2];
}

void Intersection_M2Builder::build_all(){
    build_core();
//    move_points();
    for(uint i = 3; i < Ppos_.size(); i++){
        add_point(i);
    }
    move_points();
}

void Intersection_M2Builder::add_point(uint32 pt_nb){
    VEC3 P0 = Ppos_[pt_nb];
    Dart F0 = Pdart_[pt_nb];
    M2Face face2cut;
    std::vector<VEC3> Quadp; //positions
    std::vector<Dart> Quadd; //darts

    cmap2_->foreach_cell([&](M2Face f) -> bool {
        if(in_quad(f.dart, P0)){
            Quadd = {f.dart, cmap2_->phi<1>(f.dart), cmap2_->phi<11>(f.dart), cmap2_->phi<111>(f.dart)};
            Quadp = {M2Position_[Quadd[0]], M2Position_[Quadd[1]], M2Position_[Quadd[2]], M2Position_[Quadd[3]]};
            face2cut = f; return false;
        }
        return true;
    });

    VEC3 P1 = M2FacePoints_[M2Face(face2cut)];
    Dart F1 = M2FaceBranch_[M2Face(face2cut)];

    Dart cut0, cut1;
    VEC3 AC = (Quadp[2] - Quadp[0]).normalized();
    VEC3 BD = (Quadp[3] - Quadp[1]).normalized();
    VEC3 P0P1 = (P1 - P0).normalized();
    if(abs(AC.dot(P0P1)) < abs(BD.dot(P0P1))){
        cut0 = Quadd[0]; cut1 = Quadd[2];
    } else {
        cut0 = Quadd[1]; cut1 = Quadd[3];
    }
    cmap2_->cut_face(cut0, cut1);
    M2Vertex v = cmap2_->cut_edge(M2Edge(cmap2_->phi_1(cut0)));
    M2Position_[v] = project_on_sphere((P0 + P1) * Scalar(0.5));

    Dart newFace, oldFace;
    VEC3 out0 = M2Position_[v] - center_;
    VEC3 out1 = ((P0 - M2Position_[cut0]).normalized().cross((P1 - M2Position_[cut0]).normalized()));
    if(out1.dot(out0) >= 0){
        newFace = cut0; oldFace = cut1;
    }
    else{
        newFace = cut1; oldFace = cut0;
    }

    M2FaceBranch_[M2Face(newFace)] = F0;
    M2FaceBranch_[M2Face(oldFace)] = F1;
    M2FacePoints_[M2Face(newFace)] = P0;
    M2FacePoints_[M2Face(oldFace)] = P1;

    move_points(newFace);
}

VEC3 Intersection_M2Builder::project_on_sphere(VEC3 P){
    return center_ + (P - center_).normalized() * radius_;
}

VEC3 Intersection_M2Builder::mean_dir(M2Vertex vert){
    using Quat = Eigen::Quaterniond;

    uint32 valence = cmap2_->nb_darts_of_orbit(M2Vertex(vert));

    std::vector<Quat> rotations;
    rotations.reserve(valence);

    VEC3 start = (M2Position_[vert] - center_).normalized();
    cmap2_->foreach_dart_of_orbit(vert, [&](Dart face){
        VEC3 dir = (M2FacePoints_[M2Face(face)] - center_).normalized();
        Quat q = Quat::FromTwoVectors(start, dir);
        q.normalize();
        rotations.push_back(q);
    });

    Eigen::MatrixXd m(4, valence);
    for(uint32 i = 0; i < valence; ++i){
        const Quat& q = rotations[i];
        m.col(i) = VEC4(q.w(), q.x(), q.y(), q.z());
    }

    Eigen::MatrixXd mm = m * m.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mm);
    VEC4 r = es.eigenvectors().col(3);

    Quat mean_rot(r[0], r[1], r[2], r[3]);
    mean_rot.normalize();
    return mean_rot._transformVector(start) * radius_ + center_;

//    VEC3 mean_dir = (M2Position_[vert] - center_).normalized();
//    for(uint32 i = 0; i < 50; ++i){
//    Eigen::Matrix3d rot;
//    rot.setIdentity();
//    cmap2_->foreach_dart_of_orbit(vert, [&](Dart face){
//        VEC3 v = (M2FacePoints_[M2Face(face)] - center_).normalized();
//        double angle = cgogn::geometry::angle(mean_dir, v);
//        angle = -1.0 * (M_PI - angle) / 100.0;
//        Eigen::AngleAxisd aa({ angle, mean_dir.cross(v) });
//        rot *= aa.toRotationMatrix();
//    });
//    mean_dir = rot * mean_dir;
//    mean_dir.normalize();
//    }
//    return -mean_dir * radius_ + center_;
}

//VEC3 Intersection_M2Builder::mean_dir2(M2Vertex vert){

//}

void Intersection_M2Builder::move_point(Dart vert){
//    uint32 valence = cmap2_->nb_darts_of_orbit(M2Vertex(vert));
//    if(valence == 2){
//        VEC3 V0 = M2FacePoints_[vert];
//        VEC3 V1 = M2FacePoints_[cmap2_->phi<2>(vert)];
//        M2Position_[vert] = project_on_sphere((V0 + V1) / Scalar(2));
//        return;
//    }
//    if(valence == 3){
//        VEC3 V0 = M2FacePoints_[vert];
//        VEC3 V1 = M2FacePoints_[cmap2_->phi<2>(vert)];
//        VEC3 V2 = M2FacePoints_[cmap2_->phi<212>(vert)];
//        VEC3 N = ((V2 - V0).normalized()).cross((V1 - V0).normalized()).normalized();
//        M2Position_[vert] = center_ + N * radius_;
//        return;
//    }
//    if(valence > 3){
////        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic > coord(3, valence);
////        uint32 i = 0;
////        cmap2_->foreach_dart_of_orbit(M2Vertex(vert), [&](Dart face){
////            coord.col(i++) = M2FacePoints_[face];
////        });

////        VEC3 centroid (coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());
////        coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);
////        auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
////        VEC3 N = svd.matrixU().rightCols<1>();
////        cgogn_log_info("result ") << N[0] << " " << N[1] << " " << N[2];
////        if(N.dot(centroid - center_) < 0) N = -N;

////        M2Position_[vert] = center_ + N * radius_;
//        for(uint32 i = 0; i < 10; ++i){
//        M2Position_[vert] = mean_dir(M2Vertex(vert));
//        }
//        return;
//    }
//    for(uint32 i = 0; i < 30; ++i){
    M2Position_[vert] = mean_dir(M2Vertex(vert));
//    }
}

cgogn::Dart Intersection_M2Builder::convex_quad(Dart f){
    Dart res;
    VEC3 A = M2Position_[M2Vertex(f)];
    VEC3 B = M2Position_[M2Vertex(cmap2_->phi<1>(f))];
    VEC3 C = M2Position_[M2Vertex(cmap2_->phi<11>(f))];
    VEC3 D = M2Position_[M2Vertex(cmap2_->phi<111>(f))];
    VEC3 AB = B - A;
    VEC3 AC = C - A;
    VEC3 AD = D - A;
    VEC3 N0 = AC.cross(AD);
    VEC3 N1 = AB.cross(AC);
    VEC3 N = N0.cross(N1).normalized();

    if(N.dot(AC) < 0){
        cgogn_log_info("concave: ") << f;
        res = cmap2_->phi<1>(f);
    }
    else{
        res = f;
        cgogn_log_info("convexe: ") << f;
    }
    return res;
}

void Intersection_M2Builder::move_points(Dart face){
    cmap2_->foreach_dart_of_orbit(M2Face(face), [&](Dart vert){
        move_point(vert);
    });
}

void Intersection_M2Builder::move_points(){
    cmap2_->foreach_cell([&](M2Face face){
//        convex_quad(face.dart);
    });
}

void Intersection_M2Builder::clear(){

}

} // namespace schnapps
