#include <schnapps/plugins/tubular_mesh/vessels_builder.h>

#include <cgogn/geometry/types/eigen.h>

namespace schnapps
{

Vessels_Builder::Vessels_Builder() {}

void Vessels_Builder::set_skeleton(UGraph* ug){
    ug_ = ug;
    ug_marker_ = new UGDartMarker(*ug_);
    UGNormals_ = ug_->add_attribute<VEC3, UGVertex>("normals");
    UGTangents_ = ug_->add_attribute<VEC3, UGVertex>("tangents");
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
    analyse_graph();
    build_cmap3();
    clean_up();
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
            case 3: // intersection 3
                extremities_.push_back(v);
                intersections_.push_back(v);
            break;
            default:
            break;
        }
    });
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

void Vessels_Builder::build_cmap3(){
    cgogn_log_info("Building cmap3");
    UGposition_ = ug_->template get_attribute<VEC3, UGVertex>("position");
    UGConnections_ = ug_->template add_attribute<Dart, UndirectedGraph::CDart>("connections");
    build_intersections();
    compute_tangents();
    build_branches();
}

void Vessels_Builder::build_intersections(){
    for(UGVertex ugv : intersections_){
        build_intersection3(ugv);
    }
}

void Vessels_Builder::build_intersection3(UGVertex ugv){
    Scalar radius = Scalar(0.75f);
    VEC3 Ctr = UGposition_[ugv];
    std::vector<VEC3> P, PV, C, Q, M, E, F;
    std::vector<Dart> PD;

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

//    /// Projection sur la shère de centre Ctr de l'embranchement
//    for(uint32 i = 0; i < E.size(); i++){
//        E[i] = (E[i] - Ctr).normalized();
//        E[i] = Ctr + E[i] * radius;
//    }

//    for(uint32 i = 0; i < F.size(); i++){
//        F[i] = (F[i] - Ctr).normalized();
//        F[i] = Ctr + F[i] * radius;
//    }

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
    Scalar radius = Scalar(0.75f);
    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
//    // entrance connection darts and
//    std::vector<Dart> D0 = {D, cmap3_->phi<1232>(D), cmap3_->phi<11232>(D), cmap3_->phi<111232>(D), cmap3_->phi<1112321232>(D), cmap3_->phi<11123211232>(D)};
//    std::vector<uint32> E0;
//    E0.push_back(cmap3_->embedding(D0[2], M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(D0[1], M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<1>(D0[0]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(D0[0], M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<111>(D0[3]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(D0[4], M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<111>(D0[5]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<11>(D0[5]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<11>(D0[4]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<11>(D0[2]), M3Vertex::ORBIT));
//    E0.push_back(cmap3_->embedding(cmap3_->phi<11>(D0[1]), M3Vertex::ORBIT));

    Dart d0 = d;
    Dart d1 = ug_->alpha0(d0);
    uint32 valence = ug_->nb_darts_of_orbit(UGVertex(d1));

    while(valence <= 2){
        cgogn_log_info("valence: ") << valence << " " << d0 << " " << d1 ;
        VEC3 T = UGTangents_[d1];
        VEC3 N = UGNormals_[d1];
        VEC3 M = (N.cross(T)).normalized();
        cgogn_log_info("N: ") << N[0] << " " << N[1] << " " << N[2];
        cgogn_log_info("M: ") << M[0] << " " << M[1] << " " << M[2];
        cgogn_log_info("T: ") << T[0] << " " << T[1] << " " << T[2];
        UGVertex v0 = ug_->add_vertex();
        UGVertex v1 = ug_->add_vertex();
        UGVertex v2 = ug_->add_vertex();
        UGVertex v3 = ug_->add_vertex();
        UGVertex v4 = ug_->add_vertex();
        UGVertex v5 = ug_->add_vertex();
        UGVertex v6 = ug_->add_vertex();
        UGVertex v7 = ug_->add_vertex();
        UGVertex v8 = ug_->add_vertex();
        UGVertex v9 = ug_->add_vertex();
        UGposition_[v0] = UGposition_[d1] + M * radius;
        UGposition_[v1] = UGposition_[d1] - M * radius;
        UGposition_[v2] = UGposition_[d1] + N * radius;
        UGposition_[v3] = UGposition_[d1] - N * radius;
        UGposition_[v4] = UGposition_[d1] + (M + N).normalized() * radius;
        UGposition_[v5] = UGposition_[d1] - (M + N).normalized() * radius;
        UGposition_[v6] = UGposition_[d1] + (M - N).normalized() * radius;
        UGposition_[v7] = UGposition_[d1] - (M - N).normalized() * radius;
        UGposition_[v8] = UGposition_[d1] + N * radius * Scalar(0.5f);
        UGposition_[v9] = UGposition_[d1] - N * radius * Scalar(0.5f);


        Dart D = UGConnections_[d0];
//        cgogn_log_info("d: ") << a
        if(D.is_nil()){
            cgogn_log_info("Branch error: connection is nil") << d0 << " " << D;
            return;
        }
        else{
            cgogn_log_info("Branch good") << d0 << " " << D;
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


        UGConnections_[d0] = cmap3_->phi<2112>(D1[0]);
        cgogn_log_info("next: ") << d0 << " " << D1[0] << " " << UGConnections_[d0];
        valence = ug_->nb_darts_of_orbit(UGVertex(d1));
    }

//    if(valence == 1){
//        UGVertex v0 = ug_->add_vertex();
//        UGVertex v1 = ug_->add_vertex();
//        UGVertex v2 = ug_->add_vertex();
//        UGVertex v3 = ug_->add_vertex();
//        UGVertex v4 = ug_->add_vertex();
//        UGVertex v5 = ug_->add_vertex();
//        UGVertex v6 = ug_->add_vertex();
//        UGVertex v7 = ug_->add_vertex();
//        UGposition_[v0] = UGposition_[d1] + M * radius;
//        UGposition_[v1] = UGposition_[d1] - M * radius;
//        UGposition_[v2] = UGposition_[d1] + N * radius;
//        UGposition_[v3] = UGposition_[d1] - N * radius;
//        UGposition_[v4] = UGposition_[d1] + (M + N).normalized() * radius;
//        UGposition_[v5] = UGposition_[d1] - (M + N).normalized() * radius;
//        UGposition_[v6] = UGposition_[d1] + (M - N).normalized() * radius;
//        UGposition_[v7] = UGposition_[d1] - (M - N).normalized() * radius;
//    }
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

} // namespace schnapps
