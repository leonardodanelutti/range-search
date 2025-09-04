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
        virtual ~Node() = default;
        virtual bool is_leaf() const = 0;
    };
    
    struct InternalNode : public Node {
        int split_dimension;
        double split_coordinate;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
        
        InternalNode(int dim, double coord) : split_dimension(dim), split_coordinate(coord) {}
        
        bool is_leaf() const override { return false; }
    };
    
    struct LeafNode : public Node {
        Point<D> point;
        
        LeafNode(Point<D> pt) : point(pt) {}
        
        bool is_leaf() const override { return true; }
    };

private:
    std::unique_ptr<Node> root;
    size_t size_;
    
    // Helper function to check if a point is within a range
    bool is_point_in_range(const Point<D>& point, const Iso_box<D>& range) const {
        for (int i = 0; i < D; ++i) {
            auto coord = point.cartesian(i);
            if (coord < ::min_coord<D>(range, i) || coord > ::max_coord<D>(range, i)) {
                return false;
            }
        }
        return true;
    }
    
    // Build the KD-Tree from a set of points
    std::unique_ptr<Node> build_tree(std::vector<Point<D>>& points, int depth, int start, int end) {
        if (start >= end) return nullptr;
        
        // If we have only one point, create a leaf
        if (end - start == 1) {
            return std::make_unique<LeafNode>(points[start]);
        }
        
        int dimension = depth % D;
        int median = start + (end - start) / 2;
        
        // Use nth_element to find the median point in the current dimension
        std::nth_element(points.begin() + start, points.begin() + median, points.begin() + end,
                         [dimension](const Point<D>& a, const Point<D>& b) {
                             return a.cartesian(dimension) < b.cartesian(dimension);
                         });

        // Get the split coordinate from the median point
        double split_coord = CGAL::to_double(points[median].cartesian(dimension));
        
        // Create internal node with split information
        auto internal_node = std::make_unique<InternalNode>(dimension, split_coord);
        
        // Build subtrees recursively
        internal_node->left = build_tree(points, depth + 1, start, median);
        internal_node->right = build_tree(points, depth + 1, median, end);
        
        return internal_node;
    }
    
    // Helper function for range search
    void range_search(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;
        
        if (node->is_leaf()) {
            // This is a leaf node, check the single point in it
            LeafNode* leaf = static_cast<LeafNode*>(node);
            const auto& point = leaf->point;
            
            // Check if current point is in range
            if (is_point_in_range(point, range)) {
                result.push_back(point);
            }
        } else {
            // This is an internal node, use split coordinate to decide which subtrees to search
            InternalNode* internal = static_cast<InternalNode*>(node);
            
            int dim = internal->split_dimension;
            double split_value = internal->split_coordinate;
            
            // Search left subtree if range overlaps
            auto min_val = CGAL::to_double(::min_coord<D>(range, dim));
            if (min_val <= split_value) {
                range_search(internal->left.get(), range, result);
            }
            
            // Search right subtree if range overlaps
            auto max_val = CGAL::to_double(::max_coord<D>(range, dim));
            if (max_val >= split_value) {
                range_search(internal->right.get(), range, result);
            }
        }
    }

public:
    // Empty tree
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
    
    // Get root node for visualization purposes
    Node* getRoot() const { return root.get(); }
    
    // Debug function to print tree structure
    void print_tree_structure(Node* node = nullptr, int depth = 0) const {
        if (node == nullptr) node = root.get();
        if (!node) return;
        
        std::string indent(depth * 2, ' ');
        
        if (node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            std::cout << indent << "Leaf: (";
            for (int i = 0; i < D; ++i) {
                if (i > 0) std::cout << ", ";
                std::cout << CGAL::to_double(leaf->point.cartesian(i));
            }
            std::cout << ")" << std::endl;
        } else {
            InternalNode* internal = static_cast<InternalNode*>(node);
            std::cout << indent << "Internal: dim=" << internal->split_dimension 
                      << " split=" << internal->split_coordinate << std::endl;
            
            if (internal->left) {
                std::cout << indent << "Left:" << std::endl;
                print_tree_structure(internal->left.get(), depth + 1);
            }
            if (internal->right) {
                std::cout << indent << "Right:" << std::endl;
                print_tree_structure(internal->right.get(), depth + 1);
            }
        }
    }
};
