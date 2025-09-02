#include <iostream>
#include "geometry.h"
#include "algorithms.h"

int main() {
    std::cout << "Range Search Algorithms Comparison" << std::endl;

    // Create a new kd-tree instance for 2D points
    std::vector<Point2D> points = { Point2D(1, 2), Point2D(3, 4), Point2D(5, 6) };
    KDTree<2> tree(points);

    // Create a 2D bounding box using proper CGAL construction
    Point2D min_point(0, 0);
    Point2D max_point(6, 7);
    
    // Create box using kernel construct function
    Kernel<2> kernel;
    Iso_box2D box = kernel.construct_iso_box_d_object()(min_point, max_point);
    
    std::vector<Point<2>> result = tree.range_search(box);

    std::cout << "Found " << result.size() << " points in range" << std::endl;

    // TODO: Implement algorithm selection and benchmarking
    return 0;
}
