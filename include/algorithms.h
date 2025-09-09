// algorithms.h
// Declarations for range search algorithms

#pragma once

#include "geometry.h"
#include <memory>
#include <vector>
#include <algorithm>

template<int D>
class OrthogonalRangeSearchAlgorithm {
public:
    virtual ~OrthogonalRangeSearchAlgorithm() = default;
    virtual std::vector<Point<D>> range_search(const Iso_box<D>& range) const = 0;
    virtual size_t size() const = 0;
    virtual bool empty() const = 0;
};

// KD-Tree
template<int D>
class KDTree : public OrthogonalRangeSearchAlgorithm<D> {
public:
    struct Node {
        virtual ~Node() = default;
        virtual bool is_leaf() const = 0;
    };
    
    struct InternalNode : public Node {
        int split_dimension;
        typename Kernel<D>::FT split_coordinate;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
        
        InternalNode(int dim, typename Kernel<D>::FT coord) : split_dimension(dim), split_coordinate(coord) {}
        
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
        auto split_coord = points[median].cartesian(dimension);
        
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
            LeafNode* leaf = static_cast<LeafNode*>(node);
            const auto& point = leaf->point;
            
            // Check if current point is in range
            if (::is_point_in_range<D>(point, range)) {
                result.push_back(point);
            }
        } else {
            // This is an internal node, use split coordinate to decide which subtrees to search
            InternalNode* internal = static_cast<InternalNode*>(node);
            
            int dim = internal->split_dimension;
            auto split_value = internal->split_coordinate;
            
            // Search left subtree if range overlaps
            auto min_val = ::min_coord<D>(range, dim);
            if (min_val <= split_value) {
                range_search(internal->left.get(), range, result);
            }
            
            // Search right subtree if range overlaps
            auto max_val = ::max_coord<D>(range, dim);
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
                      << " split=" << CGAL::to_double(internal->split_coordinate) << std::endl;
            
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

template<int D>
class RangeTree : public OrthogonalRangeSearchAlgorithm<D> {
public:
    struct Node {
        virtual ~Node() = default;
        virtual bool is_leaf() const = 0;

        // Dimension this node splits on
        int dimension;
    };

    struct InternalNode : public Node {
        typename Kernel<D>::FT split_value;
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
        // Associated structure
        std::unique_ptr<RangeTree<D>> assoc;

        InternalNode(typename Kernel<D>::FT val, int dim) : split_value(val) {
            this->dimension = dim;
        }
        bool is_leaf() const override { return false; }
    };

    struct LeafNode : public Node {
        Point<D> point;

        LeafNode(Point<D> pt, int dim) : point(pt) {
            this->dimension = dim;
        }
        bool is_leaf() const override { return true; }
    };

private:
    std::unique_ptr<Node> root;
    size_t size_;

    // Build the Range Tree from a set of points starting at given dimension
    std::unique_ptr<Node> build_tree(std::vector<Point<D>>& points, int dim = 0) {
        if (points.empty()) return nullptr;

        if (dim >= D) return nullptr;

        // If we have only one point, create a leaf
        if (points.size() == 1) {
            auto leaf = std::make_unique<LeafNode>(points[0], dim);
            return leaf;
        }

        // Use nth_element to find median in current dimension and partition points
        int mid = points.size() / 2;

        std::nth_element(points.begin(), points.begin() + mid, points.end(),
                         [dim](const Point<D>& a, const Point<D>& b) {
                             return a.cartesian(dim) < b.cartesian(dim);
                         });

        typename Kernel<D>::FT split_val = points[mid].cartesian(dim);

        // Create internal node
        auto node = std::make_unique<InternalNode>(split_val, dim);

        std::vector<Point<D>> left_points(points.begin(), points.begin() + mid);
        std::vector<Point<D>> right_points(points.begin() + mid, points.end());

        // Build subtrees recursively
        node->left = build_tree(left_points, dim);
        node->right = build_tree(right_points, dim);

        // Build associated structure for next dimension if not at last dimension
        if (dim + 1 < D) {
            node->assoc = std::make_unique<RangeTree<D>>(points, dim + 1);
        }

        return node;
    }

    Node* find_split_node(Node* node, typename Kernel<D>::FT min_val, typename Kernel<D>::FT max_val) const {
        while (node && !node->is_leaf()) {
            InternalNode* internal = static_cast<InternalNode*>(node);
            auto x_v = internal->split_value;
            if (max_val <= x_v || min_val >= x_v) {
                if (max_val <= x_v) {
                    node = internal->left.get();
                } else {
                    node = internal->right.get();
                }
            } else {
                break;
            }
        }
        return node;
    }

    // Helper function for range search
    void range_search_helper(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;

        // Get the range for the current dimension
        int dim = node->dimension;
        auto [min_val, max_val] = ::get_dimension_range<D>(range, dim);
        
        // Find split node
        Node* v_split = find_split_node(node, min_val, max_val);

        if (v_split->is_leaf()) {
            // Report the point if in range
            LeafNode* leaf = static_cast<LeafNode*>(v_split);
            if (::is_point_in_range<D>(leaf->point, range)) {
                result.push_back(leaf->point);
            }
        } else {
            InternalNode* internal = static_cast<InternalNode*>(v_split);
            auto split_value = internal->split_value;
            
            // Follow the path to min_val
            Node* v = internal->left.get();
            while (v && !v->is_leaf()) {
                InternalNode* v_internal = static_cast<InternalNode*>(v);
                if (min_val <= v_internal->split_value) {
                    // Call the query on associated structure of right child
                    if (v_internal->right) {
                        if (v_internal->dimension < D - 1) {
                            InternalNode* right_internal = static_cast<InternalNode*>(v_internal->right.get());
                            if (right_internal->assoc) {
                                auto assoc_results = right_internal->assoc->range_search(range);
                                result.insert(result.end(), assoc_results.begin(), assoc_results.end());
                            }
                        } else {
                            // Last dimension, search directly in right subtree
                            range_search_helper(v_internal->right.get(), range, result);
                        }
                    }
                    v = v_internal->left.get();
                } else {
                    v = v_internal->right.get();
                }
            }
            // Check if the point stored at v must be reported
            if (v && v->is_leaf()) {
                LeafNode* leaf = static_cast<LeafNode*>(v);
                if (::is_point_in_range<D>(leaf->point, range)) {
                    result.push_back(leaf->point);
                }
            }

            // Follow the path to max_val
            v = internal->right.get();
            while (v && !v->is_leaf()) {
                InternalNode* v_internal = static_cast<InternalNode*>(v);
                if (max_val >= v_internal->split_value) {
                    // Call the query on associated structure of left child
                    if (v_internal->left) {
                        if (v_internal->dimension < D - 1) {
                            InternalNode* left_internal = static_cast<InternalNode*>(v_internal->left.get());
                            if (left_internal->assoc) {
                                auto assoc_results = left_internal->assoc->range_search(range);
                                result.insert(result.end(), assoc_results.begin(), assoc_results.end());
                            }
                        } else {
                            // Last dimension, search directly in left subtree
                            range_search_helper(v_internal->left.get(), range, result);
                        }
                    }
                    v = v_internal->right.get();
                } else {
                    v = v_internal->left.get();
                }
            }
            // Check if the point stored at v must be reported
            if (v && v->is_leaf()) {
                LeafNode* leaf = static_cast<LeafNode*>(v);
                if (::is_point_in_range<D>(leaf->point, range)) {
                    result.push_back(leaf->point);
                }
            }
        }
    }

public:
    // Empty tree
    RangeTree() : size_(0) {}
    
    // Constructor that builds tree from vector of points starting at given dimension
    RangeTree(std::vector<Point<D>> points, int start_dim = 0) : size_(points.size()) {
        if (size_ == 0) {
            root = nullptr;
            return;
        }
        root = build_tree(points, start_dim);
    }
    
    // Range search: find all points within the range query
    std::vector<Point<D>> range_search(const Iso_box<D>& range) const override {
        std::vector<Point<D>> result;
        range_search_helper(root.get(), range, result);
        return result;
    }
    
    // Get size of the tree
    size_t size() const override { return size_; }
    
    // Check if tree is empty
    bool empty() const override { return size_ == 0; }
    
    // Get root node for visualization purposes
    Node* getRoot() const { return root.get(); }
};
