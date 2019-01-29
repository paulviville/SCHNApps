#include <schnapps/plugins/tubular_mesh/vessels_builder.h>

#include <cgogn/geometry/types/eigen.h>

namespace schnapps
{

Vessels_Builder::Vessels_Builder() {}

void Vessels_Builder::set_skeleton(UGraph* ug){
    ug_ = ug;
    ug_marker_ = new UGDartMarker(*ug_);
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

    ug_->foreach_adjacent_vertex_through_edge(ugv, [&](UGVertex ugv2){
        VEC3 v = UGposition_[ugv2] - Ctr;
        v.normalize();
        P.push_back(Ctr + v * Scalar(radius));
    });

    VEC3 V = (P[1] - P[0]).cross(P[2] - P[0]);
    V.normalize();
    Q.push_back(Ctr + V * Scalar(radius));
    Q.push_back(Ctr - V * Scalar(radius));

    PV.push_back((P[1]-P[0]));
    PV.push_back((P[2]-P[1]));
    PV.push_back((P[0]-P[2]));
    PV[0].normalize();
    PV[1].normalize();
    PV[2].normalize();
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


    /// building CMAP3
    auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
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

    for(VEC3 v: Q){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Qi.push_back(l);
    }

    for(VEC3 v: M){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Mi.push_back(l);
    }

    for(VEC3 v: E){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Ei.push_back(l);
    }

    for(VEC3 v: F){
        uint32 l = ca.insert_lines<1>();
        (*M3position_)[l] = v;
        Fi.push_back(l);
    }

    std::vector<Dart> H;

    H.push_back(add_hexa({Mi[0], Ei[2], Fi[2], Ei[3], Ei[0], Ctri, Ci[0], Fi[0]}));
    H.push_back(add_hexa({Ei[0], Ctri, Ci[0], Fi[0], Mi[1], Ei[1], Fi[1], Ei[4]}));
    H.push_back(add_hexa({Ctri, Ei[2], Fi[2], Ci[0], Ei[1], Mi[2], Ei[5], Fi[1]}));
    H.push_back(add_hexa({Ei[3], Fi[2], Ei[5], Qi[0], Fi[0], Ci[0], Fi[1], Ei[4]}));
    H.push_back(add_hexa({Ei[2], Mi[0], Ei[6], Fi[5], Ctri, Ei[0], Fi[3], Ci[1]}));
    H.push_back(add_hexa({Ei[0], Ctri, Ei[1], Mi[1], Fi[3], Ci[1], Fi[4], Ei[7]}));
    H.push_back(add_hexa({Ctri, Ei[2], Mi[2], Ei[1], Ci[1], Fi[5], Ei[8], Fi[4]}));
    H.push_back(add_hexa({Fi[5], Ei[6], Qi[1], Ei[8], Ci[1], Fi[3], Ei[7], Fi[4]}));

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
}

void Vessels_Builder::build_branches(){

}

void Vessels_Builder::build_branch(uint32 branch_nb){

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
    {
        m3builder_->template set_orbit_embedding<M3Vertex>(M3Vertex2(dv), verts[index++]);
    }

    return d;
}

void Vessels_Builder::clean_up(){
    ug_->remove_attribute("connections");
}

} // namespace schnapps
