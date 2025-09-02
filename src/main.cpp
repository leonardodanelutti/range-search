#include <iostream>
#include "geometry.h"
#include "algorithms.h"

int main() {
    std::cout << "Range Search Algorithms Comparison" << std::endl;

    // Create a new kd-tree instance for 2D points
    std::vector<Point2D> points = { 
        Point2D(1, 2), Point2D(3, 4), Point2D(5, 6), 
        Point2D(2, 1), Point2D(4, 3), Point2D(6, 5),
        Point2D(0, 0), Point2D(7, 7)
    };
    KDTree<2> tree(points);

    // Create a range query using CGAL
    Point2D min_point(1, 1);
    Point2D max_point(5, 5);
    
    // Create box using kernel construct function
    Kernel<2> kernel;
    Iso_box2D box = kernel.construct_iso_box_d_object()(min_point, max_point);
    
    std::vector<Point<2>> result = tree.range_search(box);

    std::cout << "Found " << result.size() << " points in range" << std::endl;
    
    std::cout << "\nKD-Tree structure (values only in leaves):" << std::endl;
    tree.print_tree_structure();

    // TODO: Implement algorithm selection and benchmarking
    return 0;
}
