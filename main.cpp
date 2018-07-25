#include "lagrange_quad.hpp"
#include "lagrange_hex.hpp"

int main()
{
  // Lagrange quads
  for(unsigned order = 1; order <= 10; ++order)
  {
    std::vector<unsigned> connectivity;
    std::vector<std::valarray<double>> points;
    std::vector<double> point_data;

    generate_lagrange_quad(order, connectivity, points);
    dst_to_center_quad(points, point_data);
    write_quad_vtk("quad_p=" + std::to_string(order) + ".vtk",
                   connectivity, points, point_data);
  }

  // Lagrange hexes
  for(unsigned order = 1; order <= 10; ++order)
  {
    std::vector<unsigned> connectivity;
    std::vector<std::valarray<double>> points;
    std::vector<double> point_data;

    generate_lagrange_hex(order, connectivity, points);
    dst_to_center_hex(points, point_data);
    write_hex_vtk("hex_p=" + std::to_string(order) + ".vtk",
                   connectivity, points, point_data);
  }

  return 0;
}
