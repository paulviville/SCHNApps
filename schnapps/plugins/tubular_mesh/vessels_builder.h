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

namespace vessels_building{
    using Scalar = cgogn::geometry::ScalarOf<VEC3>;
    using Frame = Eigen::Matrix3d;
    using Eigen::AngleAxisd;
    using Dart = cgogn::Dart;

    using UGraph = cgogn::UndirectedGraph;
    using UGDartMarker = cgogn::DartMarker<UGraph>;
    using UGDart = UGraph::CDart;
    using UGVertex = UGraph::Vertex;
    using UGEdge = UGraph::Edge;
    using UGCache = cgogn::CellCache<UGraph>;

    using CMap2 = cgogn::CMap2;
    using M2Builder = CMap2::Builder;
    using M2Dart = CMap2::CDart;
    using M2Vertex = CMap2::Vertex;
    using M2Edge = CMap2::Edge;
    using M2Face = CMap2::Face;
    using M2Volume = CMap2::Volume;

    using CMap3 = cgogn::CMap3;
    using M3DartMarker = cgogn::DartMarker<CMap3>;
    using M3Edge = CMap3::Edge;
    using M3Face = CMap3::Face;
    using M3Volume = CMap3::Volume;
    using M3Vertex = CMap3::Vertex;
    using M3Vertex2 = CMap3::Vertex2;
    using M3Builder = CMap3::Builder;

    using Branch = std::pair<Dart, Dart>;

    typedef struct Graph {
        std::vector<Branch> branches;
        std::vector<UGVertex> intersections;
        bool valid;
    } Graph;

    typedef struct UG_Attributes {
        UGraph::VertexAttribute<VEC3> position;
        UGraph::VertexAttribute<Scalar> radii;
        UGraph::CDartAttribute<Dart> m2_interface;
        UGraph::VertexAttribute<Dart> m2_CC;
        UGraph::CDartAttribute<Frame> frames;
        UGraph::CDartAttribute<Dart> sections;
//        UGraph::EdgeAttribute<uint> layers;
    } UG_Attributes;

    typedef struct CMap2_Attributes{
        CMap2::VertexAttribute<VEC3> position;
        CMap2::EdgeAttribute<VEC3> edge_pos;
        CMap2::CDartAttribute<Dart> connections;
        CMap2::VolumeAttribute<VEC3> center;

        CMap2::FaceAttribute<VEC3> f_point;
        CMap2::FaceAttribute<Dart> f_branch;
    } CMap2_Attributes;



    /// Build hexmesh from input undirected graph
    /// Cmap3 empty cmap3 to build in
    /// Cmap2 empty cmap2 to build in for debugging purposes
    /// returns false if build fails
    bool build_hexmesh(UGraph& ug, CMap3& cmap3, CMap2& cmap2);
    /// interpolation of frames test
    bool build_hexmesh2(UGraph& ug, CMap3& cmap3, CMap2& cmap2);

    /// PLACEHOLDER
    /// remove when graphs come with radius attribute
    void add_radius(UGraph& ug, Scalar radius);

    /// prints all pertinent information of the graph
    /// returns false if invalid graph
    bool ug_stats(const UGraph& ug);

    /// Gather undirected graph attributes
    bool ug_get_attributes(UGraph& ug, UG_Attributes& ug_attribs);
    /// verification, subdivision, and processing, of the graph pre construction
    /// returns false if invalid graph
    bool ug_preprocessing(UGraph& ug, UG_Attributes& ug_attribs);
    /// eliminate undesirable elements
    bool ug_clean(UGraph& ug, UG_Attributes& ug_attribs);
    /// subdivides the graph
    bool ug_sudivide(UGraph& ug, UG_Attributes& ug_attribs);

    /// Collection of relevant information: branches and intersections
    /// returns graph with simplified information
    /// if error occurs: graph.valid == false
    Graph ug_analysis(const UGraph& ug);

    /// create all required attributes for the cmap2
    bool m2_add_attributes(CMap2& cmap2, CMap2_Attributes& m2_attribs);

    /// remove all temporary attributes
//    void attribs_clean_up(UG_Attributes& ug_attribs, CMap2_Attributes& m2_attribs);
    void m2_attribs_clean_up(CMap2& cmap2, CMap2_Attributes& m2_attribs);
    void ug_attribs_clean_up(UGraph& ug,UG_Attributes& ug_attribs);

    /// Build connection interfaces topology in cmap2 (+geometry for branching point vertices)
    /// returns false if an interface fails to be created
    bool build_connection_interfaces(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2);
    bool build_interface_1(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv);
    bool build_interface_2(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv);
    bool build_interface_n(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, UGVertex ugv);
    bool build_intersection_interfaces(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2);
    bool build_1_2_interfaces(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2);

    bool complete_intersections(UGraph& ug, UG_Attributes& ug_attribs, const Graph& graph, CMap2& cmap2, CMap2_Attributes& m2_attribs);
    bool complete_intersection(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv);
    bool complete_intersection_3(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv);
    bool complete_intersection_n(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv);
    bool create_intersection_frames(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, UGVertex ugv);

    /// Propagate frames from branching points using double reflection
    bool propagate_frames(UGraph& ug, UG_Attributes& ug_attribs, const Graph& graph, CMap2& cmap2);
    bool propagate_frame_n_1(UGraph& ug, UG_Attributes& ug_attribs, Dart d);
    bool propagate_frame_n_n(UGraph& ug, UG_Attributes& ug_attribs, Dart d, CMap2& cmap2);

    /// use frames to set geometry of connection interfaces
    bool set_interfaces_geometry(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs);
    bool set_edge_geometry(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs);

    /// builds a section for each branch in the graph and stores connecting darts in cmap2
    bool build_branch_sections(UGraph& ug, UG_Attributes& ug_attribs, CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3);
    /// creates a section of a branch and returns d0
    Dart new_section(CMap3& cmap3);
    /// sews all sections of the cmap3
    bool sew_sections(CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3);
    /// set geometry of cmap3 from cmap2
    bool set_m3_geometry(CMap2& cmap2, CMap2_Attributes& m2_attribs, CMap3& cmap3);

    /// returns first dart for which the diagonal starting at it's base would produce a convex oriented set of triangles from a quad
    Dart convex_quad(CMap2& cmap2, CMap2_Attributes& m2_attribs, Dart f);
    /// average of points on a sphere
    VEC3 mean_dir(VEC3 center, Scalar radius, VEC3 point, std::vector<VEC3> points);
    /// returns coordinates of the projection of P on the sphere of center C ands radius R
    VEC3 project_on_sphere(VEC3 P, Scalar R, VEC3 C);

    Dart shift_interface(CMap2& cmap2, Dart m2f, uint nb_shifts);
    Frame shift_frame(Frame frame, uint nb_shifts);


    /// POST treatment
    void subdivide(CMap3& cmap3);
    void add_layer_edge(UGraph& ug, UG_Attributes& ug_attribs, CMap3& cmap3, UGEdge uge);
    void add_layer_edge_corner(UGraph& ug, UG_Attributes& ug_attribs, CMap3& cmap3, UGEdge uge);
}

} // namespace schnapps

#endif // VESSELS_BUILDER_H_
