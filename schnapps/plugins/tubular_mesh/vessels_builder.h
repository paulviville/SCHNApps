#ifndef VESSELS_BUILDER_H_
#define VESSELS_BUILDER_H_

#include <cgogn/core/cmap/cmap3.h>
#include <cgogn/core/cmap/cmap3_builder.h>
#include <cgogn/core/graph/undirected_graph.h>
#include <cgogn/core/graph/undirected_graph_builder.h>

#include <cgogn/geometry/functions/orientation.h>
#include <schnapps/core/schnapps.h>

namespace schnapps
{
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

public:
    Vessels_Builder();
    void set_skeleton(UGraph* ug);
    void set_cmap3(CMap3* cmap3);
    bool compute_cmap3();

private:
    CMap3* cmap3_;
    M3Builder* m3builder_;
    UGraph* ug_;
    UGDartMarker* ug_marker_;

    std::vector<UGVertex> extremities_;
    std::vector<UGVertex> intersections_;
    std::vector<std::pair<Dart, Dart>> branches_;
    UndirectedGraph::VertexAttribute<VEC3> UGposition_;
//    CMap3::ChunkArrayContainer* M3Cac_;
    CMap3::ChunkArray<VEC3>* M3position_;
    UndirectedGraph::CDartAttribute<Dart> UGConnections_;
    UndirectedGraph::VertexAttribute<VEC3> UGNormals_;
    UndirectedGraph::VertexAttribute<VEC3> UGTangents_;
//    UndirectedGraph::EdgeAttribute
    void analyse_graph();
    void find_ends();
    void find_branches();
    void check_isolated_loops();

    void build_cmap3();
    void build_intersections();
    void build_intersection3(UGVertex ugv);
    void build_branches();
    void build_branch(Dart d);
    void compute_tangents();
    void compute_tangent_3_1(Dart d);
    void propagate_normals();
//    void compute_joint_normals();

    Dart add_hexa();
    Dart add_hexa(std::vector<uint32> verts);
    void embed_hexa(Dart d, std::vector<uint32> verts);
    void clean_up();
};

} // namespace schnapps

#endif // VESSELS_BUILDER_H_
