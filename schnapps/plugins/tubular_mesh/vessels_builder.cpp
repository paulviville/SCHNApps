#include <schnapps/plugins/tubular_mesh/vessels_builder.h>

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
}

bool Vessels_Builder::compute_cmap3(){ // filler for tests
        if (! cmap3_->is_embedded<M3Vertex>())
                m3builder_->create_embedding<M3Vertex::ORBIT>();
        // recup du container d'attributs de sommet
        auto& ca = m3builder_->attribute_container<M3Vertex::ORBIT>();
        // recup de l'attribut
        CMap3::ChunkArray<VEC3>* position = ca.get_chunk_array<VEC3>("position");
        if (position == nullptr)
            position = ca.add_chunk_array<VEC3>("position");

        VEC3 Ps[] = {{ 0.5, -0.5, -0.5},
                    { 0.5, 0.5, -0.5},
                    { 0.5, 0.5, 0.5},
                    { 0.5, -0.5, 0.5},
                    { -0.5, -0.5, -0.5},
                    { -0.5, 0.5, -0.5},
                    { -0.5, 0.5, 0.5},
                    { -0.5, -0.5, 0.5}};

        std::vector<uint32> verts;
        for (const auto& p:Ps)
        {
            auto l = ca.insert_lines<1>();
            (*position)[l] = p;
            verts.push_back(l);
        }

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
    return true;
}

void Vessels_Builder::analyse_graph(){
    find_ends();
    find_branches();
    cgogn_log_info("graph analysis:\n") << "nb inter: " << intersections_.size() << "   nb branches: " << branches_.size();
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
            if(!ug_->is_boundary(d1))
                branches_.push_back(std::pair<Dart, Dart>(d0, d1));
        }
    }
}

void Vessels_Builder::check_isolated_loops(){

}

void Vessels_Builder::build_intersections(){

}

void Vessels_Builder::build_intersection3(uint32 inter_nb){

}

void Vessels_Builder::build_branches(){

}

void Vessels_Builder::build_branch(uint32 branch_nb){

}

} // namespace schnapps
