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
