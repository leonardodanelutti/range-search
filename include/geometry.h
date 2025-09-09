#pragma once
#define CGAL_EIGEN3_ENABLED

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
// #include <CGAL/Kernel_d/Iso_box_d.h>

// Use CGAL's built-in dimension-parameterized kernels

// Inexact kernel
// template<int D>
// using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<D>>;

// Exact kernel
template<int D>
using Kernel = CGAL::Epeck_d<CGAL::Dimension_tag<D>>;

// Point types with compile-time dimension
template<int D>
using Point = typename Kernel<D>::Point_d;

using Point2D = Point<2>;
using Point3D = Point<3>;

// Iso_box (axis-aligned bounding box) types
template<int D>
using Iso_box = typename Kernel<D>::Iso_box_d;

using Iso_box2D = Iso_box<2>;
using Iso_box3D = Iso_box<3>;

// Helper functions for Iso_box coordinate access
template<int D>
auto min_coord(const Iso_box<D>& box, int dimension) {
    Kernel<D> kernel;
    auto min_pt = kernel.construct_min_vertex_d_object()(box);
    return min_pt.cartesian(dimension);
}

template<int D>
auto max_coord(const Iso_box<D>& box, int dimension) {
    Kernel<D> kernel;
    auto max_pt = kernel.construct_max_vertex_d_object()(box);
    return max_pt.cartesian(dimension);
}

// Get the range (min, max) for a specific dimension of a D-dimensional box
template<int D>
std::pair<typename Kernel<D>::FT, typename Kernel<D>::FT> get_dimension_range(const Iso_box<D>& box, int dimension) {
    auto min_val = ::min_coord<D>(box, dimension);
    auto max_val = ::max_coord<D>(box, dimension);
    return std::make_pair(min_val, max_val);
}

// Helper function to check if a point is within a range
template<int D>
bool is_point_in_range(const Point<D>& point, const Iso_box<D>& range) {
    for (int i = 0; i < D; ++i) {
        auto coord = point.cartesian(i);
        if (coord < ::min_coord<D>(range, i) || coord > ::max_coord<D>(range, i)) {
            return false;
        }
    }
    return true;
}

// Helper function to create Iso_box from coordinate ranges more concisely
template<int D>
Iso_box<D> make_iso_box(std::initializer_list<double> min_coords, std::initializer_list<double> max_coords) {
    Kernel<D> kernel;
    return kernel.construct_iso_box_d_object()(
        Point<D>(min_coords.begin(), min_coords.end()),
        Point<D>(max_coords.begin(), max_coords.end())
    );
}
