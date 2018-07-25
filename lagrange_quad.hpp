// Alexander Grayver, 2018
#ifndef LAGRANGE_QUAD_H
#define LAGRANGE_QUAD_H

#include <vector>
#include <cmath>
#include <valarray>
#include <fstream>

/*
 * Modified from
 * https://github.com/Kitware/VTK/blob/265ca48a79a36538c95622c237da11133608bbe5/Common/DataModel/vtkLagrangeQuadrilateral.cxx#L558
 */
int PointIndexFromIJK_Quad(int i, int j, const int order[2])
{
  bool ibdy = (i == 0 || i == order[0]);
  bool jbdy = (j == 0 || j == order[1]);
  // How many boundaries do we lie on at once?
  int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

  if (nbdy == 2) // Vertex DOF
  { // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0));
  }

  int offset = 4;
  if (nbdy == 1) // Edge DOF
  {
    if (!ibdy)
    { // On i axis
      return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;
      }

      if (!jbdy)
      { // On j axis
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;
    }
  }

  offset += 2 * (order[0] - 1 + order[1] - 1);
  // nbdy == 0: Face DOF
  return offset + (i - 1) + (order[0] - 1) * ((j - 1));
}

/*
 * Modified from
 * https://github.com/Kitware/VTK/blob/675adbc0feeb3f62730ecacb2af87917af124543/Filters/Sources/vtkCellTypeSource.cxx#L1307
 */
void generate_lagrange_quad(const int cell_order,
                            std::vector<unsigned>& conn,
                            std::vector<std::valarray<double>>& points)
{
  const int numPtsPerCell = (cell_order + 1) * (cell_order + 1);

  conn.resize(numPtsPerCell);
  points.resize(numPtsPerCell, std::valarray<double>(2));
  const int order[2] = { cell_order, cell_order };
  const int numPtsPerDim[2] = { cell_order + 1, cell_order + 1};
  const std::vector<std::valarray<double>>& nodes = {{0., 0.},
                                                     {1., 0.},
                                                     {1., 1.},
                                                     {0., 1.}};

  for (int n = 0; n <= order[1]; ++n)
  {
    for (int m = 0; m <= order[0]; ++m)
    {
      int connectivity_idx = PointIndexFromIJK_Quad(m, n, order);
      int lexicographic_idx = n * numPtsPerDim[0] + m;

      const double r = static_cast<double>(m) / order[0];
      const double s = static_cast<double>(n) / order[1];
      const std::valarray<double>& pm =
          (1.0 - r) * (nodes[3] * s + nodes[0] * (1.0 - s)) +
          r *         (nodes[2] * s + nodes[1] * (1.0 - s));

      conn[connectivity_idx] = lexicographic_idx;
      points[lexicographic_idx] = pm;
    }
  }
}

void dst_to_center_quad(const std::vector<std::valarray<double>>& points,
                        std::vector<double> &distance_to_center)
{
  distance_to_center.reserve(points.size());

  const std::valarray<double> center = {0.5, 0.5};
  for(auto &point: points)
  {
    auto v = point - center;
    distance_to_center.emplace_back(sqrt(v[0]*v[0] + v[1]*v[1]));
  }
}

void write_quad_vtk(const std::string &file_name,
                    const std::vector<unsigned>& connectivity,
                    const std::vector<std::valarray<double>>& points,
                    const std::vector<double>& point_data)
{
  std::ofstream ofs(file_name);

  ofs << "# vtk DataFile Version 4.2\n"
      << "vtk output\n"
      << "ASCII\n"
      << "DATASET UNSTRUCTURED_GRID\n"
      << "POINTS " << points.size() << " double\n";

  for(auto &p: points)
    ofs << p[0] << " " <<  p[1] << " 0\n";

  ofs << "CELLS 1 " << connectivity.size() + 1 << "\n";
  ofs << connectivity.size();

  for(auto &c: connectivity)
    ofs << " " << c;

  const int VTK_LAGRANGE_QUAD = 70;
  ofs << "\nCELL_TYPES 1\n" << VTK_LAGRANGE_QUAD << "\n";

  ofs << "POINT_DATA " << point_data.size() << "\n"
      << "FIELD FieldData 1\n"
      << "DistanceToCenter 1 " << point_data.size() << " float\n";
  for(auto &d: point_data)
    ofs << d << " ";

  ofs << "\n";

  ofs.close();
}

#endif
