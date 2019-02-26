#ifndef VESSELS_BUILDER_H_
#define VESSELS_BUILDER_H_

#include <cgogn/core/cmap/cmap3.h>
#include <cgogn/core/cmap/cmap3_builder.h>
#include <cgogn/core/graph/undirected_graph.h>
#include <cgogn/core/graph/undirected_graph_builder.h>
#include <cgogn/core/cmap/cmap2.h>
#include <cgogn/core/cmap/cmap2_builder.h>

#include <schnapps/core/schnapps.h>

namespace schnapps
{
class Intersection_M2Builder{
    using Scalar = cgogn::geometry::ScalarOf<VEC3>;
    using UGraph = cgogn::UndirectedGraph;
    using UGVertex = UGraph::Vertex;
    using UGEdge = UGraph::Edge;
    using Dart = cgogn::Dart;
    using UGCellCache = cgogn::CellCache<UGraph>;
    using CMap2 = cgogn::CMap2;
    using M2Builder = CMap2::Builder;
    using M2Vertex = CMap2::Vertex;
    using M2Face = CMap2::Face;
    using M2Edge = CMap2::Edge;

public:
    Intersection_M2Builder();
    void set_cmap2(CMap2* cmap2);
    void set_ugraph(UGraph* ug);
    void build_intersection(UGVertex ugv);
    Dart convex_quad(Dart f);
    VEC3 project_on_sphere(VEC3 p);

private:
    CMap2* cmap2_;
    UGraph* ug_;
    M2Builder* m2builder_;
    UndirectedGraph::VertexAttribute<VEC3> UGposition_;
    UndirectedGraph::CDartAttribute<Dart> UGConnections_;
    CMap2::VertexAttribute<VEC3> M2Position_;
    CMap2::FaceAttribute<VEC3> M2FacePoints_;
    CMap2::FaceAttribute<Dart> M2FaceBranch_;

    UGVertex ugv_;
    VEC3 center_;
    Scalar radius_;
    std::vector<VEC3> Ppos_;
    std::vector<Dart> Pdart_;

    bool in_quad(Dart face, VEC3 P);
    VEC3 mean_dir(M2Vertex vert);
    VEC3 mean_dir2(M2Vertex vert);
    void get_intersection_data();
    void build_core();
    void build_all();
    void add_point(uint32 pt_nb);
    void move_points(Dart face);
    void move_point(Dart vert);
    void move_points();
    void clear();
};

class Vessels_Builder{
    using Scalar = cgogn::geometry::ScalarOf<VEC3>;
    using CMap3 = cgogn::CMap3;
    using M3Builder = CMap3::Builder;
    using UGraph = cgogn::UndirectedGraph;
    using UGVertex = UGraph::Vertex;
    using UGEdge = UGraph::Edge;
    using Dart = cgogn::Dart;
    using UGDartMarker = cgogn::DartMarker<UGraph>;
    using UGCellCache = cgogn::CellCache<UGraph>;
    using M3Vertex = CMap3::Vertex;
    using M3Vertex2 = CMap3::Vertex2;
    using CMap2 = cgogn::CMap2;
    using M2Builder = CMap2::Builder;
    using M2Vertex = CMap2::Vertex;
    using M2Face = CMap2::Face;
    using M2Edge = CMap2::Edge;
public:
    Vessels_Builder();
    void set_skeleton(UGraph* ug);
    void set_cmap3(CMap3* cmap3);
    bool compute_cmap3();

    CMap2* cmap2_;
    M2Builder* m2builder_;

private:
    CMap3* cmap3_;
    M3Builder* m3builder_;
    UGraph* ug_;
    UGDartMarker* ug_marker_;
    Intersection_M2Builder intersection_m2builder_;

    std::vector<UGVertex> extremities_;
    std::vector<UGVertex> intersections_;
    std::vector<std::pair<Dart, Dart>> branches_;
    UndirectedGraph::VertexAttribute<VEC3> UGposition_;
//    CMap3::ChunkArrayContainer* M3Cac_;
    CMap3::ChunkArray<VEC3>* M3position_;
    UndirectedGraph::CDartAttribute<Dart> UGConnections_;
//    UndirectedGraph::VertexAttribute<VEC3> UGNormals_;
    UndirectedGraph::CDartAttribute<VEC3> UGNormals_;
    UndirectedGraph::VertexAttribute<VEC3> UGTangents_;
//    UndirectedGraph::EdgeAttribute
    void analyse_graph();
    void find_ends();
    void find_branches();
    void check_isolated_loops();

    void subdivide_graph();
    void build_cmap3();
    void build_intersections();
    void build_intersection3(UGVertex ugv);
    void build_intersectionN(UGVertex ugv);
    void build_branches();
    void build_branch(Dart d);
    void compute_tangents();
    void compute_tangent_3_1(Dart d);
    void propagate_normals();
//    void compute_joint_normals();

    Dart create_block(M2Face f);
    Dart add_hexa();
    Dart add_hexa(std::vector<uint32> verts);
    void embed_hexa(Dart d, std::vector<uint32> verts);
    bool in_quad(VEC3 A, VEC3 B, VEC3 C, VEC3 D, VEC3 P, VEC3 CTR);

    void clean_up();
};



} // namespace schnapps

#endif // VESSELS_BUILDER_H_
