/*******************************************************************************
* SCHNApps                                                                     *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#include <schnapps/plugins/tubular_mesh/tubular_mesh.h>
#include <cgogn/core/cmap/cmap3.h>
#include <cgogn/core/cmap/cmap3_builder.h>
#include <schnapps/plugins/cmap_provider/cmap_provider.h>
#include <schnapps/plugins/import/import.h>
#include <schnapps/plugins/polyline_render/polyline_render.h>

#include <schnapps/core/schnapps.h>
#include <cgogn/geometry/types/geometry_traits.h>
#include <cgogn/geometry/types/eigen.h>
#include <cgogn/geometry/functions/orientation.h>

#include <schnapps/plugins/tubular_mesh/vessels_builder.h>
#include <cgogn/core/graph/undirected_graph_builder.h>

namespace schnapps
{

namespace plugin_tubular_mesh
{
using UGVertex = UndirectedGraph::Vertex;
using UGEdge = UndirectedGraph::Edge;
using Dart = cgogn::Dart;
using UGDartMarker = cgogn::DartMarker<UndirectedGraph>;
using UGCellCache = cgogn::CellCache<UndirectedGraph>;
using M3Vertex = CMap3::Vertex;
using M3Vertex2 = CMap3::Vertex2;

Plugin_TubularMesh::Plugin_TubularMesh() :
	plugin_import_(nullptr),
	plugin_cmap_provider_(nullptr),
    ugh_(nullptr),
    map3h_(nullptr)
{
	this->name_ = SCHNAPPS_PLUGIN_NAME;
}

QString Plugin_TubularMesh::plugin_name()
{
	return SCHNAPPS_PLUGIN_NAME;
}


bool Plugin_TubularMesh::enable()
{
	plugin_import_ = static_cast<plugin_import::Plugin_Import*>(schnapps_->enable_plugin(plugin_import::Plugin_Import::plugin_name()));
	plugin_cmap_provider_ = static_cast<plugin_cmap_provider::Plugin_CMapProvider*>(schnapps_->enable_plugin(plugin_cmap_provider::Plugin_CMapProvider::plugin_name()));
    plugin_polyline_render_ = static_cast<plugin_polyline_render::Plugin_PolylineRender*>(schnapps_->enable_plugin(plugin_polyline_render::Plugin_PolylineRender::plugin_name()));

    ugh_ = plugin_import_->import_graph_from_file("/home/viville/Data/two_intersections.cg");
    plugin_polyline_render_->set_edge_color(schnapps_->selected_view(), ugh_, QColor(255,255,255), true);
    plugin_polyline_render_->set_vertex_scale_factor(schnapps_->selected_view(), ugh_, 0.1f, true);
//    schnapps_->selected_view().link_object(ugh_test_, true);
    UndirectedGraph* ug = ugh_->map();
    UndirectedGraph::VertexAttribute<VEC3> UGposition = ug->template get_attribute<VEC3, UGVertex>("position");
//    UndirectedGraph::CDartAttribute
    map3h_ = plugin_cmap_provider_->add_cmap3("vessels");
    CMap3* map3 = map3h_->map();

    Vessels_Builder VBuilder;
    VBuilder.set_skeleton(ug);
    VBuilder.set_cmap3(map3);


    UGDartMarker dart_marker(*ug);
    cgogn_log_info("NB Darts:") << ug->nb_darts();

    UGCellCache initial_cache(*ug);
    initial_cache.template build<UGEdge>();
    initial_cache.template build<UGVertex>();

    std::vector<UGVertex> extremities;
    std::vector<UGVertex> intersections;
    std::vector<std::pair<Dart, Dart>> branches;

    VBuilder.compute_cmap3();
    map3h_->create_vbo("position");
    map3h_->set_bb_vertex_attribute(add_setting("Bounding box attribute", "position").toString());
    map3h_->notify_attribute_added(M3Vertex::ORBIT, "position");
    map3h_->notify_connectivity_change();

    return true;
}

void Plugin_TubularMesh::disable()
{

}

} // namespace plugin_tubular_mesh

} // namespace schnapps
