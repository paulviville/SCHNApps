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

#ifndef SCHNAPPS_PLUGIN_TUBULAR_MESH_H_
#define SCHNAPPS_PLUGIN_TUBULAR_MESH_H_

#include <schnapps/plugins/tubular_mesh/dll.h>

#include <schnapps/core/plugin_processing.h>

#include <QAction>

namespace schnapps
{

namespace plugin_cmap_provider
{
class Plugin_CMapProvider;
class CMap2Handler;
class CMap3Handler;
class UndirectedGraphHandler;
}

namespace plugin_import { class Plugin_Import; }

namespace plugin_tubular_mesh
{

using CMap2Handler = plugin_cmap_provider::CMap2Handler;
using CMap3Handler = plugin_cmap_provider::CMap3Handler;
using UndirectedGraphHandler = plugin_cmap_provider::UndirectedGraphHandler;

/**
* @brief Plugin for tubular mesh generation
*/
class SCHNAPPS_PLUGIN_TUBULAR_MESH_API Plugin_TubularMesh : public PluginProcessing
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "SCHNApps.Plugin")
	Q_INTERFACES(schnapps::Plugin)

public:

	Plugin_TubularMesh();
	inline ~Plugin_TubularMesh() override {}
	static QString plugin_name();

private:

	bool enable() override;
	void disable() override;

public:


private:

	plugin_import::Plugin_Import* plugin_import_;
	plugin_cmap_provider::Plugin_CMapProvider* plugin_cmap_provider_;

	UndirectedGraphHandler* ugh_;
};

} // namespace plugin_tubular_mesh

} // namespace schnapps

#endif // SCHNAPPS_PLUGIN_TUBULAR_MESH_H_
