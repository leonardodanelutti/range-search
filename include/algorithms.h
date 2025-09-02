// algorithms.h
// Declarations for range search algorithms

#pragma once

#include "geometry.h"
#include <memory>
#include <vector>
#include <algorithm>

// KD-Tree
template<int D>
class KDTree {
public:
    struct Node {
        Point<D> point;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
        // Keep track of the current dimension
        int split_dimension;
        
        Node(Point<D> p, int dim) : point(p), split_dimension(dim) {}
    };

private:
    std::unique_ptr<Node> root;
    size_t size_;
    
    // Build the KD-Tree from a set of points
    std::unique_ptr<Node> build_tree(std::vector<Point<D>>& points, int depth, int start, int end) {
        if (start >= end) return nullptr;
        
        int dimension = depth % D;

        int median = start + (end - start) / 2;
        
        // Use nth_element to find the median point in the current dimension
        std::nth_element(points.begin() + start, points.begin() + median, points.begin() + end,
                         [dimension](const Point<D>& a, const Point<D>& b) {
                             return a.cartesian(dimension) < b.cartesian(dimension);
                         });

        // Create node and recursively build subtrees
        std::unique_ptr<Node> node = std::make_unique<Node>(points[median], dimension);
        
        node->left = build_tree(points, depth + 1, start, median);
        node->right = build_tree(points, depth + 1, median + 1, end);
        
        return node;
    }
    
    // Helper function for range search
    void range_search(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;
        
        // Check if current point is in range
        bool in_range = true;
        for (int i = 0; i < D; ++i) {
            auto coord = node->point.cartesian(i);
            auto min_val = ::min_coord<D>(range, i);
            auto max_val = ::max_coord<D>(range, i);
            if (coord < min_val || coord > max_val) {
                in_range = false;
                break;
            }
        }
        
        if (in_range) {
            result.push_back(node->point);
        }
        
        // Recursively search children
        int dim = node->split_dimension;
        auto split_value = node->point.cartesian(dim);
        
        // Search left subtree if range overlaps
        if (::min_coord<D>(range, dim) <= split_value) {
            range_search(node->left.get(), range, result);
        }
        
        // Search right subtree if range overlaps
        if (::max_coord<D>(range, dim) >= split_value) {
            range_search(node->right.get(), range, result);
        }
    }

public:
    // Default constructor - creates empty tree
    KDTree(): size_(0) {}
    
    // Constructor that builds tree from vector of points
    KDTree(std::vector<Point<D>> points) : size_(points.size()) {
        if (size_ == 0) {
            root = nullptr;
            return;
        }
        root = build_tree(points, 0, 0, static_cast<int>(size_));
    }
    
    // Range search: find all points within the range query
    std::vector<Point<D>> range_search(const Iso_box<D>& range) const {
        std::vector<Point<D>> result;
        range_search(root.get(), range, result);
        return result;
    }
    
    // Get size of the tree
    size_t size() const { return size_; }
    
    // Check if tree is empty
    bool empty() const { return size_ == 0; }
};
