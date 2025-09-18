// algorithms.h
// Declarations for range search algorithms

#pragma once

#include "geometry.h"
#include <memory>
#include <vector>
#include <algorithm>
#include <optional>

template<int D>
class OrthogonalRangeSearch {
public:
    virtual ~OrthogonalRangeSearch() = default;
    virtual std::vector<Point<D>> range_search(const Iso_box<D>& range) const = 0;
    virtual size_t size() const = 0;
    virtual bool empty() const = 0;
    virtual size_t memory_usage() const = 0;  // Memory usage in bytes
};

// KD-Tree
template<int D>
class KDTree : public OrthogonalRangeSearch<D> {
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
    
    // Calculate memory usage of the tree
    size_t memory_usage() const override {
        return calculate_node_memory(root.get()) + sizeof(*this);
    }
    
private:
    // Helper function to calculate memory usage of nodes
    size_t calculate_node_memory(Node* node) const {
        if (!node) return 0;
        
        if (node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            return sizeof(LeafNode);
        } else {
            InternalNode* internal = static_cast<InternalNode*>(node);
            return sizeof(InternalNode) + 
                   calculate_node_memory(internal->left.get()) +
                   calculate_node_memory(internal->right.get());
        }
    }
    
public:
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
class RangeTree : public OrthogonalRangeSearch<D> {
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

    // Find the node where the search paths to min_val and max_val diverge
    Node* find_split_node(Node* node, typename Kernel<D>::FT min_val, typename Kernel<D>::FT max_val) const {
        while (node && !node->is_leaf()) {
            InternalNode* internal = static_cast<InternalNode*>(node);
            auto x_v = internal->split_value;
            if (max_val <= x_v) {
                node = internal->left.get();
            } else if (min_val >= x_v) {
                node = internal->right.get();
            } else {
                break;
            }
        }
        return node;
    }

    // Helper to check and report leaf node if in range
    bool check_and_report_leaf(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (node && node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            if (::is_point_in_range<D>(leaf->point, range)) {
                result.push_back(leaf->point);
            }
            return true;
        }
        return false;
    }

    // Helper to report all the values in a subtree
    void report_tree(Node* node, std::vector<Point<D>>& result) const {
        if (!node) return;
        
        if (node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            result.push_back(leaf->point);
        } else {
            InternalNode* internal = static_cast<InternalNode*>(node);
            report_tree(internal->left.get(), result);
            report_tree(internal->right.get(), result);
        }
    }

    // Helper to search in the associated structure 
    void search_assoc(Node* subtree, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!subtree) return;
        
        if (!check_and_report_leaf(subtree, range, result)) {
            InternalNode* subtree_internal = static_cast<InternalNode*>(subtree);
            if (subtree_internal->dimension < D - 1 && subtree_internal->assoc) {
                auto assoc_results = subtree_internal->assoc->range_search(range);
                result.insert(result.end(), assoc_results.begin(), assoc_results.end());
            } else {
                // Last dimension, report all nodes in this subtree
                report_tree(subtree, result);
            }
        }
    }

    // Helper function for range search
    void range_search_helper(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;

        // Get the range for the current dimension
        int dim = node->dimension;
        auto [min_val, max_val] = ::get_dimension_range<D>(range, dim);
        
        // Find split node
        Node* v_split = find_split_node(node, min_val, max_val);

        if (check_and_report_leaf(v_split, range, result)) {
            return;
        }

        InternalNode* internal = static_cast<InternalNode*>(v_split);
        auto split_value = internal->split_value;
        
        // Follow the path to min_val
        Node* v = internal->left.get();
        while (v && !v->is_leaf()) {
            InternalNode* v_internal = static_cast<InternalNode*>(v);
            if (min_val <= v_internal->split_value) {
                // Call the query on associated structure of right child
                search_assoc(v_internal->right.get(), range, result);
                v = v_internal->left.get();
            } else {
                v = v_internal->right.get();
            }
        }
        check_and_report_leaf(v, range, result);

        // Follow the path to max_val
        v = internal->right.get();
        while (v && !v->is_leaf()) {
            InternalNode* v_internal = static_cast<InternalNode*>(v);
            if (max_val >= v_internal->split_value) {
                // Call the query on associated structure of left child
                search_assoc(v_internal->left.get(), range, result);
                v = v_internal->right.get();
            } else {
                v = v_internal->left.get();
            }
        }
        check_and_report_leaf(v, range, result);
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
    
    // Calculate memory usage of the tree
    size_t memory_usage() const override {
        return calculate_node_memory_rangetree(root.get()) + sizeof(*this);
    }

    // Get root node for visualization purposes
    Node* getRoot() const { return root.get(); }
    
private:
    // Helper function to calculate memory usage of nodes for RangeTree
    size_t calculate_node_memory_rangetree(Node* node) const {
        if (!node) return 0;
        
        if (node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            return sizeof(LeafNode);
        } else {
            InternalNode* internal = static_cast<InternalNode*>(node);
            size_t memory = sizeof(InternalNode);
            memory += calculate_node_memory_rangetree(internal->left.get());
            memory += calculate_node_memory_rangetree(internal->right.get());
            // Add memory for associated structure
            if (internal->assoc) {
                memory += internal->assoc->memory_usage();
            }
            return memory;
        }
    }
};

// Base class for associated structures
template<int D>
class AssociatedStructure {
public:
    virtual ~AssociatedStructure() = default;
    virtual size_t memory_usage() const = 0;
};

// CascadeArray - A cascade array is an array of points sorted by one dimension
// Each element of the array has two pointer to another elements in another cascade arrays
template<int D>
class CascadeArray : public AssociatedStructure<D> {
public:
    struct Element {
        Point<D> point;
        // Pointer is just the index of the element in the array, or std::nullopt if none
        std::optional<size_t> left;
        std::optional<size_t> right;
        Element(Point<D> pt) : point(pt), left(std::nullopt), right(std::nullopt) {}
    };

    std::vector<Element> elements;

    // Builder for single element
    CascadeArray(Point<D> pt) {
        elements.emplace_back(pt);
    }

    // Given two CascadeArrays, merge them into one
    CascadeArray(CascadeArray<D>* left_array, CascadeArray<D>* right_array) {
        
        elements.reserve(left_array->elements.size() + right_array->elements.size());

        size_t i = 0, j = 0;
        size_t p_i = 0, p_j = 0;
        while (i < left_array->elements.size() && j < right_array->elements.size()) {
            if (left_array->elements[i].point.cartesian(D - 1) < right_array->elements[j].point.cartesian(D - 1)) {
                Element new_elem(left_array->elements[i].point);
                
                // Update pointers
                new_elem.left = p_i;
                new_elem.right = p_j;

                elements.push_back(std::move(new_elem));
                ++i;
                // If consecutive points have same coordinate in this dimension, keep the pointer to the first one
                if (left_array->elements[i].point.cartesian(D - 1) != left_array->elements[i - 1].point.cartesian(D - 1)) {
                    p_i = i;
                }
            } else {
                Element new_elem(right_array->elements[j].point);

                // Update pointers
                new_elem.left = p_i;
                new_elem.right = p_j;

                elements.push_back(std::move(new_elem));
                ++j;
                // If consecutive points have same coordinate in this dimension, keep the pointer to the first one
                if (right_array->elements[j].point.cartesian(D - 1) != right_array->elements[j - 1].point.cartesian(D - 1)) {
                    p_j = j;
                }
            }
        }
        
        // Handle remaining elements from left child array
        while (i < left_array->elements.size()) {
            Element new_elem(left_array->elements[i].point);

            // Update pointers
            new_elem.left = p_i;

            elements.push_back(std::move(new_elem));
            ++i;
            // If consecutive points have same coordinate in this dimension, keep the pointer to the first one
            if (i < left_array->elements.size() && left_array->elements[i].point.cartesian(D - 1) != left_array->elements[i - 1].point.cartesian(D - 1)) {
                p_i = i;
            }
        }
        
        // Handle remaining elements from right child array
        while (j < right_array->elements.size()) {
            Element new_elem(right_array->elements[j].point);

            // Update pointers
            new_elem.right = p_j;
            
            elements.push_back(std::move(new_elem));
            ++j;
            // If consecutive points have same coordinate in this dimension, keep the pointer to the first one
            if (j < right_array->elements.size() && right_array->elements[j].point.cartesian(D - 1) != right_array->elements[j - 1].point.cartesian(D - 1)) {
                p_j = j;
            }
        }
    }

    // Returns pointer to the smallest element greater or equal that value
    // If none found, returns nullptr
    std::optional<size_t> find_smallest_greater_or_equal(typename Kernel<D>::FT value, int dim = D-1) {
        int left = 0, right = static_cast<int>(elements.size());
        std::optional<size_t> result = std::nullopt;
        while (left < right) {
            int mid = left + (right - left) / 2;
            if (elements[mid].point.cartesian(dim) >= value) {
                result = mid;
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        return result;
    }
    
    // Calculate memory usage of the cascade array
    size_t memory_usage() const override {
        return sizeof(*this) + elements.size() * sizeof(Element);
    }
};

// LayeredRangeTree
template<int D>
class LayeredRangeTree : public OrthogonalRangeSearch<D>, public AssociatedStructure<D> {
    static_assert(D > 1, "LayeredRangeTree requires D > 1 (use KDTree or RangeTree for D=1)");
    
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
        // Associated structure:
        // - LayeredRangeTree when dimension < D-2
        // - CascadeArray when dimension == D-2
        std::unique_ptr<AssociatedStructure<D>> assoc;

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

        // Build associated structure for next dimension
        if (dim + 2 == D) {
            // NEW
            // Create CascadeArray for the last dimension
            CascadeArray<D>* left_cascade = nullptr;
            CascadeArray<D>* right_cascade = nullptr;
            
            // Cascade array for the left subtree
            if (node->left && node->left->is_leaf()) {
                LeafNode* left_leaf = static_cast<LeafNode*>(node->left.get());
                // Create a temporary cascade array for the single point
                left_cascade = new CascadeArray<D>(left_leaf->point);
            } else if (node->left) {
                InternalNode* left_internal = static_cast<InternalNode*>(node->left.get());
                if (left_internal->assoc) {
                    left_cascade = static_cast<CascadeArray<D>*>(left_internal->assoc.get());
                }
            }
            
            // Cascade array for the right subtree
            if (node->right && node->right->is_leaf()) {
                LeafNode* right_leaf = static_cast<LeafNode*>(node->right.get());
                // Create a temporary cascade array for the single point
                right_cascade = new CascadeArray<D>(right_leaf->point);
            } else if (node->right) {
                InternalNode* right_internal = static_cast<InternalNode*>(node->right.get());
                if (right_internal->assoc) {
                    right_cascade = static_cast<CascadeArray<D>*>(right_internal->assoc.get());
                }
            }
            
            // Merge the cascade arrays
            if (left_cascade && right_cascade) {
                node->assoc = std::make_unique<CascadeArray<D>>(left_cascade, right_cascade);
                // Clean up temporary arrays if they were created for leaf nodes
                if (node->left && node->left->is_leaf()) {
                    delete left_cascade;
                }
                if (node->right && node->right->is_leaf()) {
                    delete right_cascade;
                }
            } else if (left_cascade) {
                if (node->left && node->left->is_leaf()) {
                    node->assoc = std::unique_ptr<CascadeArray<D>>(left_cascade);
                } else {
                    // This shouldn't happen in a properly constructed tree
                    node->assoc = nullptr;
                }
            } else if (right_cascade) {
                if (node->right && node->right->is_leaf()) {
                    node->assoc = std::unique_ptr<CascadeArray<D>>(right_cascade);
                } else {
                    // This shouldn't happen in a properly constructed tree
                    node->assoc = nullptr;
                }
            }
        } else if (dim + 2 < D) {
            // Create LayeredRangeTree for other dimensions
            node->assoc = std::make_unique<LayeredRangeTree<D>>(points, dim + 1);
        }

        return node;
    }

    // Find the node where the search paths to min_val and max_val diverge
    Node* find_split_node(Node* node, typename Kernel<D>::FT min_val, typename Kernel<D>::FT max_val) const {
        while (node && !node->is_leaf()) {
            InternalNode* internal = static_cast<InternalNode*>(node);
            auto x_v = internal->split_value;
            if (max_val <= x_v) {
                node = internal->left.get();
            } else if (min_val >= x_v) {
                node = internal->right.get();
            } else {
                break;
            }
        }
        return node;
    }

    // Helper to check and report leaf node if in range
    bool check_and_report_leaf(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (node && node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            if (::is_point_in_range<D>(leaf->point, range)) {
                result.push_back(leaf->point);
            }
            return true;
        }
        return false;
    }

    // Helper to search in the associated structure 
    void search_assoc(Node* node, std::optional<size_t> cascade_pointer, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;
        
        if (check_and_report_leaf(node, range, result)) return;

        InternalNode* node_internal = static_cast<InternalNode*>(node);
        AssociatedStructure<D>* assoc = node_internal->assoc.get();

        if (node_internal->dimension < D - 2 && assoc) {
            // For LayeredRangeTree associated structures
            LayeredRangeTree<D>* layered_tree = static_cast<LayeredRangeTree<D>*>(assoc);
            auto assoc_results = layered_tree->range_search(range);
            result.insert(result.end(), assoc_results.begin(), assoc_results.end());
        } else if (node_internal->dimension == D - 2 && assoc) {
            // NEW
            // dim D-2, use cascade array
            // Report all points in the cascade array from the cascade pointer to max_val
            CascadeArray<D>* cascade_array = static_cast<CascadeArray<D>*>(assoc);
            if (cascade_pointer.has_value()) {
                size_t i = cascade_pointer.value();
                auto last_dim_max = ::max_coord<D>(range, D - 1);
                while (i < cascade_array->elements.size() && cascade_array->elements[i].point.cartesian(D - 1) <= last_dim_max) {
                    result.push_back(cascade_array->elements[i].point);
                    ++i;
                }
            }
        }
    }

    // NEW
    // Follow the cascade pointer to the left
    std::optional<size_t> get_pointer_left(InternalNode* internal, std::optional<size_t> cascade_pointer) const {
        if (internal->dimension == D - 2 && cascade_pointer.has_value()) {
            CascadeArray<D>* cascade_array = static_cast<CascadeArray<D>*>(internal->assoc.get());
            if (cascade_array && cascade_pointer.value() < cascade_array->elements.size()) {
                return cascade_array->elements[cascade_pointer.value()].left;
            }
        }
        return std::nullopt;
    }

    // Follow the cascade pointer to the right
    std::optional<size_t> get_pointer_right(InternalNode* internal, std::optional<size_t> cascade_pointer) const {
        if (internal->dimension == D - 2 && cascade_pointer.has_value()) {
            CascadeArray<D>* cascade_array = static_cast<CascadeArray<D>*>(internal->assoc.get());
            if (cascade_array && cascade_pointer.value() < cascade_array->elements.size()) {
                return cascade_array->elements[cascade_pointer.value()].right;
            }
        }
        return std::nullopt;
    }

    // Helper function for range search
    void range_search_helper(Node* node, const Iso_box<D>& range, std::vector<Point<D>>& result) const {
        if (!node) return;

        // Get the range for the current dimension
        int dim = node->dimension;
        auto [min_val, max_val] = ::get_dimension_range<D>(range, dim);
        
        // Find split node
        Node* v_split = find_split_node(node, min_val, max_val);

        if (check_and_report_leaf(v_split, range, result)) {
            return;
        }

        InternalNode* internal = static_cast<InternalNode*>(v_split);
        auto split_value = internal->split_value;

        // NEW
        // When dimension is D-2 find the point with the smallest value that is greater
        // than the last dimension value of the range
        std::optional<size_t> cascade_pointer = std::nullopt;
        if (dim == D - 2) {
            CascadeArray<D>* cascade_array = static_cast<CascadeArray<D>*>(internal->assoc.get());
            auto last_dim_min = ::min_coord<D>(range, D - 1);
            cascade_pointer = cascade_array->find_smallest_greater_or_equal(last_dim_min);
        }
        
        // Follow the path to min_val
        Node* v = internal->left.get();
        std::optional<size_t> cascade_pointer_v = get_pointer_left(internal, cascade_pointer); // NEW
        while (v && !v->is_leaf()) {
            InternalNode* v_internal = static_cast<InternalNode*>(v);
            if (min_val <= v_internal->split_value) {
                search_assoc(v_internal->right.get(), get_pointer_right(v_internal, cascade_pointer_v), range, result);
                v = v_internal->left.get();
                cascade_pointer_v = get_pointer_left(v_internal, cascade_pointer_v); // NEW
            } else {
                v = v_internal->right.get();
                cascade_pointer_v = get_pointer_right(v_internal, cascade_pointer_v); // NEW
            }
        }
        check_and_report_leaf(v, range, result);

        // Follow the path to max_val
        v = internal->right.get();
        cascade_pointer_v = get_pointer_right(internal, cascade_pointer); // NEW
        while (v && !v->is_leaf()) {
            InternalNode* v_internal = static_cast<InternalNode*>(v);
            if (max_val >= v_internal->split_value) {
                // Call the query on associated structure of left child
                search_assoc(v_internal->left.get(), get_pointer_left(v_internal, cascade_pointer_v), range, result);
                v = v_internal->right.get();
                cascade_pointer_v = get_pointer_right(v_internal, cascade_pointer_v); // NEW
            } else {
                v = v_internal->left.get();
                cascade_pointer_v = get_pointer_left(v_internal, cascade_pointer_v); // NEW
            }
        }
        check_and_report_leaf(v, range, result);
    }

public:
    // Empty tree
    LayeredRangeTree() : size_(0) {}
    
    // Constructor that builds tree from vector of points starting at given dimension
    LayeredRangeTree(std::vector<Point<D>> points, int start_dim = 0) : size_(points.size()) {
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
    
    // Calculate memory usage of the tree
    size_t memory_usage() const override {
        return calculate_node_memory_layered(root.get()) + sizeof(*this);
    }

    // Get root node for visualization purposes
    Node* getRoot() const { return root.get(); }
    
private:
    // Helper function to calculate memory usage of nodes for LayeredRangeTree
    size_t calculate_node_memory_layered(Node* node) const {
        if (!node) return 0;
        
        if (node->is_leaf()) {
            LeafNode* leaf = static_cast<LeafNode*>(node);
            return sizeof(LeafNode);
        } else {
            InternalNode* internal = static_cast<InternalNode*>(node);
            size_t memory = sizeof(InternalNode);
            memory += calculate_node_memory_layered(internal->left.get());
            memory += calculate_node_memory_layered(internal->right.get());
            // Add memory for associated structure
            if (internal->assoc) {
                memory += internal->assoc->memory_usage();
            }
            return memory;
        }
    }
};
