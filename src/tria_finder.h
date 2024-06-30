#ifndef _TRIA_FINDER_H_
#define _TRIA_FINDER_H_ 1

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/rtree.h>

/**
 * Class for finding information related to a triangulation.
 *
 * This class provides methods for finding the closest vertex to a given point,
 * as well as finding the active cell around a given point in a triangulation.
 * It also includes methods for updating the triangulation and generating data structures
 * for efficient vertex and cell lookup.
 */
class TriaFinder {
public:
  /**
   * Constructs a TriaFinder object.
   */
  TriaFinder(const dealii::Triangulation<3> &);

  /**
   * Destructor for the TriaFinder class.
   */
  ~TriaFinder();

  /**
   * Finds the closest vertex to the given point.
   *
   * @param point The point for which to find the closest vertex.
   * @return The index of the closest vertex.
   */
  int find_closest_vertex(const dealii::Point<3> &) const;

  /**
   * Finds the active cell in a triangulation around a given point.
   *
   * @param point The point around which to find the active cell.
   * @return A pair containing the active cell iterator and the reference point in the cell.
   */
  std::pair<dealii::Triangulation<3>::active_cell_iterator, dealii::Point<3>>
  find_active_cell_around_point(const dealii::Point<3> &) const;

private:
  void update() const;

  /**
   * Generates an R-tree index of the used vertices.
   */
  void generate_used_vertices_rtree() const;

  /**
   * Generates a map that associates each vertex with the cells it belongs to.
   */
  void generate_vertex_to_cell_map() const;

private:
  const dealii::Triangulation<3> *mesh;    // Pointer to the triangulation.
  boost::signals2::connection tria_signal; // Connection to the signal emitted by the triangulation.
  mutable bool need_update;                // Flag to indicate if the data structures need to be updated.

  mutable std::vector<dealii::Triangulation<3>::active_cell_iterator> vertex_to_cells;  // Map from vertex to cells.
  mutable dealii::RTree<std::pair<dealii::Point<3>, unsigned int>> used_vertices_rtree; // R-tree index of used vertices.
};

#endif
