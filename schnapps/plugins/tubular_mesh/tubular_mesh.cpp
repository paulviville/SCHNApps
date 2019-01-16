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

#include <schnapps/plugins/cmap_provider/cmap_provider.h>
#include <schnapps/plugins/import/import.h>

#include <schnapps/core/schnapps.h>

namespace schnapps
{

namespace plugin_tubular_mesh
{

Plugin_TubularMesh::Plugin_TubularMesh() :
	plugin_import_(nullptr),
	plugin_cmap_provider_(nullptr),
	ugh_(nullptr)
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

	ugh_ = plugin_import_->import_graph_from_file("/home/kraemer/Media/Data/graph/graph.cg");

	return true;
}

void Plugin_TubularMesh::disable()
{

}

} // namespace plugin_tubular_mesh

} // namespace schnapps
