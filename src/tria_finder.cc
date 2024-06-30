#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools.h>

#include <map>
#include <queue>

#include "tria_finder.h"

/*
 * This constructor initializes the TriaFinder with a given triangulation.
*/
TriaFinder::TriaFinder(const dealii::Triangulation<3> &t) : mesh(&t) {
  tria_signal = t.signals.any_change.connect([&]() { need_update = true; });

  need_update = true;
}

/*
 * This destructor disconnects the signal from the triangulation.
*/
TriaFinder::~TriaFinder() {
  if (tria_signal.connected()) {
    tria_signal.disconnect();
  }
}

/*
 * This function generates an R-tree index of the used vertices.
*/
void TriaFinder::generate_used_vertices_rtree() const {
  // Extract the vertices used in the active cells
  std::map<unsigned int, dealii::Point<3>> used_vertices;
  for (const auto &cell : mesh->active_cell_iterators()) {
    for (unsigned int v = 0; v < cell->n_vertices(); ++v)
      used_vertices[cell->vertex_index(v)] = cell->vertex(v);
  }

  // Create an R-tree index of the used vertices
  unsigned int i = 0;
  std::vector<std::pair<dealii::Point<3>, unsigned int>> vertices(used_vertices.size());
  for (const auto &it : used_vertices)
    vertices[i++] = std::make_pair(it.second, it.first);
  used_vertices_rtree = dealii::pack_rtree(vertices);
}

/*
 * This function generates a map that associates each vertex with the cells it belongs to.
*/
void TriaFinder::generate_vertex_to_cell_map() const {
  vertex_to_cells.resize(mesh->n_vertices());
  for (const auto &cell : mesh->active_cell_iterators()) {
    for (unsigned int v = 0; v < cell->n_vertices(); v++) {
      vertex_to_cells[cell->vertex_index(v)] = cell;
    }
  }
}

/*
 * This function updates the data structures used for efficient vertex and cell lookup.
*/
void TriaFinder::update() const {
  generate_used_vertices_rtree();
  generate_vertex_to_cell_map();

  need_update = false;
}

/*
 * This function finds the closest vertex to the given point.
 * It uses an R-tree index to efficiently find the closest vertex.
*/
int TriaFinder::find_closest_vertex(const dealii::Point<3> &p) const {
  // Update the data structures in case the triangulation has changed.
  if (need_update) {
    update();
  }

  // Find the closest vertex to the given point using the R-tree index.
  std::vector<std::pair<dealii::Point<3>, unsigned int>> closest_vertex;
  used_vertices_rtree.query(boost::geometry::index::nearest(p, 1), std::back_inserter(closest_vertex));

  return closest_vertex[0].second;
}

/*
 * This function finds the active cell in a triangulation around a given point.
 * It uses a breadth-first search to find the cell containing the point.
 * It begins by finding the closest vertex to the point and then searches the
 * neighboring cells until the point is found.
*/
std::pair<dealii::Triangulation<3>::active_cell_iterator, dealii::Point<3>>
TriaFinder::find_active_cell_around_point(const dealii::Point<3> &p) const {
  // Update the data structures in case the triangulation has changed.
  if (need_update) {
    update();
  }

  // Initialize the cell_point pair
  auto cell_point = std::make_pair(dealii::Triangulation<3>::active_cell_iterator(), dealii::Point<3>());

  // Initialize a vector to keep track of the cells that have been visited
  std::vector<bool> touched_cells(mesh->n_active_cells(), false);

  bool inside_cell = false;

  // Initialize a queue to perform the breadth-first search,
  // the first cell to be visited is the one containing the closest vertex to the point.
  std::queue<dealii::Triangulation<3>::active_cell_iterator> q;
  q.push(vertex_to_cells[find_closest_vertex(p)]);

  // Perform the breadth-first search
  while (!q.empty()) {
    // Get the next cell from the queue
    auto cell = q.front();
    q.pop();

    // If the cell has already been visited, we skip it.
    if (touched_cells[cell->active_cell_index()]) {
      continue;
    }
    touched_cells[cell->active_cell_index()] = true;

    // Check if the point is inside the current cell for a tolerance of 1E-6.
    try {
      inside_cell = dealii::GeometryInfo<3>::is_inside_unit_cell(
          dealii::StaticMappingQ1<3>::mapping.transform_real_to_unit_cell(cell, p), 1E-6);
    } catch (const dealii::Mapping<3>::ExcTransformationFailed &) {
      inside_cell = false;
    }

    // If the point is inside the current cell, we have found the cell containing the point.
    // We store the cell and the reference point in the cell and break the loop.
    if (inside_cell) {
      cell_point = std::make_pair(cell, dealii::GeometryInfo<3>::project_to_unit_cell(
                                            dealii::StaticMappingQ1<3>::mapping.transform_real_to_unit_cell(cell, p)));
      break;
    }

    // Add the neighbors of the current cell to the queue if they have not been visited yet.
    std::vector<dealii::Triangulation<3>::active_cell_iterator> neighbors;
    dealii::GridTools::get_active_neighbors<dealii::Triangulation<3>>(cell, neighbors);
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
      if (!touched_cells[(*it)->active_cell_index()]) {
        q.push(*it);
      }
    }
  }

  return cell_point;
}
