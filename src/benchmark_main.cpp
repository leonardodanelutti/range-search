#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <map>
#include <cstdlib>
#include <algorithm>
#include "geometry.h"
#include "algorithms.h"

// Benchmark results
struct BenchmarkResult {
    std::string algorithm;
    int dimension;
    size_t point_count;
    std::string distribution;
    int iterations;
    double construction_time_us = 0.0;
    size_t memory_usage_bytes = 0;
    double average_query_time_us = 0.0;
};

// Timer utility
class BenchmarkTimer {
public:
    template<typename Func>
    static double time_operation(Func&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        return static_cast<double>(duration.count());
    }
};

// Distribution types
enum class DistributionType {
    UNIFORM,
    GAUSSIAN,
    CLUSTERED
};

// Algorithm types
enum class AlgorithmType {
    KDTREE,
    RANGETREE,
    LAYEREDRANGETREE
};

// Point generation utilities
template<int D>
class PointGenerator {
public:
    std::vector<Point<D>> generate_points(DistributionType dist, size_t count) {
        std::vector<Point<D>> points;
        points.reserve(count);
        
        std::random_device rd;
        std::mt19937 gen(rd());
        
        switch (dist) {
            case DistributionType::UNIFORM:
                return generate_uniform_points(count, gen);
            case DistributionType::GAUSSIAN:
                return generate_gaussian_points(count, gen);
            case DistributionType::CLUSTERED:
                return generate_clustered_points(count, gen);
        }
        return points;
    }

private:
    std::vector<Point<D>> generate_uniform_points(size_t count, std::mt19937& gen) {
        std::vector<Point<D>> points;
        std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
        
        for (size_t i = 0; i < count; ++i) {
            std::vector<double> coords(D);
            for (int d = 0; d < D; ++d) {
                coords[d] = dist(gen);
            }
            points.emplace_back(coords.begin(), coords.end());
        }
        return points;
    }
    
    std::vector<Point<D>> generate_gaussian_points(size_t count, std::mt19937& gen) {
        std::vector<Point<D>> points;
        std::normal_distribution<double> dist(0.0, 200.0);
        
        for (size_t i = 0; i < count; ++i) {
            std::vector<double> coords(D);
            for (int d = 0; d < D; ++d) {
                coords[d] = dist(gen);
            }
            points.emplace_back(coords.begin(), coords.end());
        }
        return points;
    }
    
    std::vector<Point<D>> generate_clustered_points(size_t count, std::mt19937& gen) {
        std::vector<Point<D>> points;
        std::uniform_real_distribution<double> cluster_center_dist(0.0, 800.0);
        std::normal_distribution<double> point_offset_dist(0.0, 30.0);
        std::uniform_int_distribution<int> cluster_selector(0, 4); // 5 clusters
        
        // Generate cluster centers
        std::vector<Point<D>> cluster_centers;
        for (int i = 0; i < 5; ++i) {
            std::vector<double> coords(D);
            for (int d = 0; d < D; ++d) {
                coords[d] = cluster_center_dist(gen);
            }
            cluster_centers.emplace_back(coords.begin(), coords.end());
        }
        
        for (size_t i = 0; i < count; ++i) {
            int cluster_id = cluster_selector(gen);
            Point<D> center_point = cluster_centers[cluster_id];
            
            std::vector<double> coords(D);
            for (int d = 0; d < D; ++d) {
                double offset = point_offset_dist(gen);
                coords[d] = CGAL::to_double(center_point[d]) + offset;
            }
            points.emplace_back(coords.begin(), coords.end());
        }
        return points;
    }
};

// Query generation utilities
template<int D>
class QueryGenerator {
public:
    Iso_box<D> generate_query_box(const std::vector<Point<D>>& points, size_t target_results = 1000) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        // Find bounding box of all points
        using FT = typename Kernel<D>::FT;
        std::vector<FT> min_coords(D), max_coords(D);
        
        // Initialize with first point coordinates
        for (int d = 0; d < D; ++d) {
            min_coords[d] = points[0][d];
            max_coords[d] = points[0][d];
        }
        
        // Find actual min/max coordinates
        for (const auto& point : points) {
            for (int d = 0; d < D; ++d) {
                min_coords[d] = std::min(min_coords[d], point[d]);
                max_coords[d] = std::max(max_coords[d], point[d]);
            }
        }
        
        // Generate random query box
        std::vector<double> query_min_coords(D), query_max_coords(D);

        // Create a vector of D-1 random number uniformly distributed numbers between 0.1 and 0.1*range
        std::vector<double> random_offsets;
        for (int d = 0; d < D - 1; ++d) {
            double range = CGAL::to_double(max_coords[d] - min_coords[d]);
            std::uniform_real_distribution<double> offset_dist(0.1, 0.1 * range);
            random_offsets.push_back(offset_dist(gen));
        }

        // c = target_results * range1 * range2 * ... * rangeD / number_of_points
        double c = target_results;
        for (int d = 0; d < D; ++d)
            c *= CGAL::to_double(max_coords[d] - min_coords[d]);
        c /= points.size();

    
        // Find x: 
        //  - calc x_0 = c^(1/(D))
        //  - one step of Newton's method for the function f(x) = x(x+random_offsets[0])(x+random_offsets[1])...(x+random_offsets[D-2]) - c
        double x = std::pow(c, 1.0 / D);
        double f_x = x;
        for (int d = 0; d < D - 1; ++d)
            f_x *= (x + random_offsets[d]);
        f_x -= c;
        double f_prime_x = 1.0;
        for (int d = 0; d < D - 1; ++d)
            f_prime_x *= (x + random_offsets[d]);
        double f_prime_1_x = 0.0;
        for (int d = 0; d < D - 1; ++d)
            f_prime_1_x += 1 / (x + random_offsets[d]);
        f_prime_x = f_prime_x * (1 + x*f_prime_1_x);
        x = x - f_x / f_prime_x;
        double box_size_x = x;

        double result = 1.0;
        for (int d = 0; d < D; ++d) {
            result *= box_size_x + (d == 0 ? 0.0 : random_offsets[d - 1]);
        }

        
        for (int d = 0; d < D; ++d) {
            double box_size = box_size_x;
            if (d != 0) {
                box_size += random_offsets[d - 1];
            }
            
            double min_coord_double = CGAL::to_double(min_coords[d]);
            double max_coord_double = CGAL::to_double(max_coords[d]);
            std::uniform_real_distribution<double> start_dist(min_coord_double, max_coord_double - box_size);
            query_min_coords[d] = start_dist(gen);
            query_max_coords[d] = query_min_coords[d] + box_size;
        }
        
        Point<D> query_min(query_min_coords.begin(), query_min_coords.end());
        Point<D> query_max(query_max_coords.begin(), query_max_coords.end());
        
        // Use kernel to construct Iso_box properly
        Kernel<D> kernel;
        return kernel.construct_iso_box_d_object()(query_min, query_max);
    }
};

std::string distribution_to_string(DistributionType dist) {
    switch (dist) {
        case DistributionType::UNIFORM: return "uniform";
        case DistributionType::GAUSSIAN: return "gaussian";
        case DistributionType::CLUSTERED: return "clustered";
        default: return "unknown";
    }
}

std::string algorithm_to_string(AlgorithmType alg) {
    switch (alg) {
        case AlgorithmType::KDTREE: return "KDTree";
        case AlgorithmType::RANGETREE: return "RangeTree";
        case AlgorithmType::LAYEREDRANGETREE: return "LayeredRangeTree";
        default: return "Unknown";
    }
}

DistributionType string_to_distribution(const std::string& str) {
    if (str == "uniform") return DistributionType::UNIFORM;
    if (str == "gaussian") return DistributionType::GAUSSIAN;
    if (str == "clustered") return DistributionType::CLUSTERED;
    return DistributionType::UNIFORM;
}

AlgorithmType string_to_algorithm(const std::string& str);

AlgorithmType string_to_algorithm(const std::string& str) {
    if (str == "KDTree" || str == "kdtree") return AlgorithmType::KDTREE;
    if (str == "RangeTree" || str == "rangetree") return AlgorithmType::RANGETREE;
    if (str == "LayeredRangeTree" || str == "layeredrangetree") return AlgorithmType::LAYEREDRANGETREE;
    return AlgorithmType::KDTREE;
}

// Function declarations for benchmark analysis
void generate_benchmark_analysis(const std::vector<BenchmarkResult>& results);
void create_memory_analysis_csv(const std::vector<BenchmarkResult>& results);
void create_construction_time_csv(const std::vector<BenchmarkResult>& results);
void create_query_time_csv(const std::vector<BenchmarkResult>& results);

// Main benchmark function
template<int D>
BenchmarkResult run_benchmark(DistributionType distribution, size_t point_count, 
                             AlgorithmType algorithm) {
    
    BenchmarkResult result;
    result.distribution = distribution_to_string(distribution);
    result.point_count = point_count;
    result.dimension = D;
    result.algorithm = algorithm_to_string(algorithm);
    result.iterations = 0; // Will be updated with actual count
    
    // Generate points
    PointGenerator<D> generator;
    auto points = generator.generate_points(distribution, point_count);
    
    // Generate query generator for runtime queries
    QueryGenerator<D> query_gen;
    
    // Benchmark based on algorithm type
    switch (algorithm) {
        case AlgorithmType::KDTREE: {
            // Measure construction time
            KDTree<D>* constructed_tree = nullptr;
            result.construction_time_us = BenchmarkTimer::time_operation([&]() {
                constructed_tree = new KDTree<D>(points);
            });
            
            // Measure memory usage
            result.memory_usage_bytes = constructed_tree->memory_usage();
            
            // Measure average query time (run for 10 seconds)
            KDTree<D> tree(points);
            double total_query_time = 0.0;
            int query_count = 0;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            auto end_time = start_time + std::chrono::seconds(5);
            
            while (std::chrono::high_resolution_clock::now() < end_time) {
                auto query = query_gen.generate_query_box(points, std::min(0.05 * point_count, 500.0));
                total_query_time += BenchmarkTimer::time_operation([&]() {
                    auto result_points = tree.range_search(query);
                    size_t size = result_points.size();
                    (void)size;
                });
                query_count++;
            }
            
            result.average_query_time_us = total_query_time / query_count;
            result.iterations = query_count; // Update with actual count
            
            // Clean up
            delete constructed_tree;
            break;
        }
        
        case AlgorithmType::RANGETREE: {
            // Measure construction time
            RangeTree<D>* constructed_tree = nullptr;
            result.construction_time_us = BenchmarkTimer::time_operation([&]() {
                constructed_tree = new RangeTree<D>(points);
            });
            
            // Measure memory usage
            result.memory_usage_bytes = constructed_tree->memory_usage();
            
            // Measure average query time (run for 10 seconds)
            RangeTree<D> tree(points);
            double total_query_time = 0.0;
            int query_count = 0;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            auto end_time = start_time + std::chrono::seconds(10);
            
            while (std::chrono::high_resolution_clock::now() < end_time) {
                auto query = query_gen.generate_query_box(points, std::min(0.05 * point_count, 200.0));
                total_query_time += BenchmarkTimer::time_operation([&]() {
                    auto result_points = tree.range_search(query);
                    size_t size = result_points.size();
                    (void)size;
                });
                query_count++;
            }
            
            result.average_query_time_us = total_query_time / query_count;
            result.iterations = query_count; // Update with actual count
            
            // Clean up
            delete constructed_tree;
            break;
        }
        
        case AlgorithmType::LAYEREDRANGETREE: {
            // Measure construction time
            LayeredRangeTree<D>* constructed_tree = nullptr;
            result.construction_time_us = BenchmarkTimer::time_operation([&]() {
                constructed_tree = new LayeredRangeTree<D>(points);
            });
            
            // Measure memory usage
            result.memory_usage_bytes = constructed_tree->memory_usage();
            
            // Measure average query time (run for 10 seconds)
            LayeredRangeTree<D> tree(points);
            double total_query_time = 0.0;
            int query_count = 0;
            
            auto start_time = std::chrono::high_resolution_clock::now();
            auto end_time = start_time + std::chrono::seconds(15);
            
            while (std::chrono::high_resolution_clock::now() < end_time) {
                auto query = query_gen.generate_query_box(points, std::min(0.05 * point_count, 200.0));
                total_query_time += BenchmarkTimer::time_operation([&]() {
                    auto result_points = tree.range_search(query);
                    size_t size = result_points.size();
                    (void)size;
                });
                query_count++;
            }
            
            result.average_query_time_us = total_query_time / query_count;
            result.iterations = query_count; // Update with actual count
            
            // Clean up
            delete constructed_tree;
            break;
        }
    }
    
    return result;
}

// Template dispatcher for different dimensions
BenchmarkResult run_benchmark_dispatch(int dimension, DistributionType distribution, 
                                      size_t point_count, AlgorithmType algorithm, int) {
    switch (dimension) {
        case 2:
            return run_benchmark<2>(distribution, point_count, algorithm);
        case 3:
            return run_benchmark<3>(distribution, point_count, algorithm);
        case 4:
            return run_benchmark<4>(distribution, point_count, algorithm);
        case 5:
            return run_benchmark<5>(distribution, point_count, algorithm);
        case 6:
            return run_benchmark<6>(distribution, point_count, algorithm);
        default:
            throw std::invalid_argument("Unsupported dimension: " + std::to_string(dimension));
    }
}

// Print usage information
void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n"
              << "Options:\n"
              << "  --distribution DIST    Distribution type: uniform, gaussian, clustered (default: uniform)\n"
              << "  --points COUNT         Number of points (default: 1000)\n"
              << "  --dimension DIM        Dimension: 2 or 3 (default: 2)\n"
              << "  --algorithm ALG        Algorithm: KDTree, RangeTree, LayeredRangeTree (default: KDTree)\n"
              << "  --help                 Show this help message\n"
              << "\nNote: Query benchmarking runs for 10 seconds and calculates average time.\n"
              << "\nExample:\n"
              << "  " << program_name << " --algorithm LayeredRangeTree --dimension 2 --points 5000\n";
}

// Print benchmark result
void print_result(const BenchmarkResult& result) {
    std::cout << "\n=== Benchmark Results ===\n";
    std::cout << "Algorithm: " << result.algorithm << "\n";
    std::cout << "Dimension: " << result.dimension << "D\n";
    std::cout << "Distribution: " << result.distribution << "\n";
    std::cout << "Points: " << result.point_count << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Memory usage: " << (result.memory_usage_bytes / 1024.0) << " KB\n";
    std::cout << "Construction time: " << (result.construction_time_us / 1000.0) << " ms\n";
    std::cout << "Average query time: " << (result.average_query_time_us) << " μs\n";
    std::cout << "========================\n";
}

// Enhanced default benchmark test with comprehensive analysis
void run_default_test() {
    std::cout << "Running comprehensive benchmark analysis (2D only)...\n";
    
    // Test configurations
    // std::vector<size_t> point_counts = {1125, 1688, 2531, 3797, 5695, 8543, 12814, 19222, 28833, 43249, 64873, 97310, 145965, 218948, 328421, 656842, 828421, 1000000};
    std::vector<size_t> point_counts = {10000, 20000, 40000, 60000, 80000, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000};
    std::vector<AlgorithmType> algorithms = {
        AlgorithmType::KDTREE,
        AlgorithmType::RANGETREE,
        AlgorithmType::LAYEREDRANGETREE
    };
    
    // Store all results for analysis
    std::vector<BenchmarkResult> all_results;
    
    // Run benchmarks for 2D only
    std::cout << "\n=== Testing 2D algorithms ===\n";
    
    for (size_t points : point_counts) {
        std::cout << "Points: " << points << " - ";
        
        for (auto alg : algorithms) {
            std::cout << algorithm_to_string(alg);
            auto result = run_benchmark_dispatch(2, DistributionType::UNIFORM, points, alg, 0);
            std::cout << "(" << result.iterations << " queries) ";
            all_results.push_back(result);
        }
        std::cout << "✓\n";
    }
    
    // Generate analysis and graphs
    generate_benchmark_analysis(all_results);
}

// Generate comprehensive analysis and data files for graphing
void generate_benchmark_analysis(const std::vector<BenchmarkResult>& results) {
    std::cout << "\n=== Benchmark Analysis ===\n";
    
    // Create CSV files for each metric
    create_memory_analysis_csv(results);
    create_construction_time_csv(results);
    create_query_time_csv(results);
    
    std::cout << "\nAnalysis complete! Generated files:\n";
    std::cout << "- memory_analysis.csv (Memory usage data)\n";
    std::cout << "- construction_time.csv (Construction time data)\n";
    std::cout << "- query_time.csv (Query time data)\n";
}

// Create memory usage analysis CSV
void create_memory_analysis_csv(const std::vector<BenchmarkResult>& results) {
    std::ofstream file("memory_analysis.csv");
    file << "Points,Dimension,KDTree_KB,RangeTree_KB,LayeredRangeTree_KB\n";
    
    // Group by points and dimension
    std::map<std::pair<size_t, int>, std::map<std::string, double>> memory_data;
    
    for (const auto& result : results) {
        auto key = std::make_pair(result.point_count, result.dimension);
        memory_data[key][result.algorithm] = result.memory_usage_bytes / 1024.0;
    }
    
    for (const auto& entry : memory_data) {
        size_t points = entry.first.first;
        int dim = entry.first.second;
        const auto& algorithms = entry.second;
        
        file << points << "," << dim << ",";
        file << algorithms.at("KDTree") << ",";
        file << algorithms.at("RangeTree") << ",";
        file << algorithms.at("LayeredRangeTree") << "\n";
    }
}

// Create construction time analysis CSV
void create_construction_time_csv(const std::vector<BenchmarkResult>& results) {
    std::ofstream file("construction_time.csv");
    file << "Points,Dimension,KDTree_ms,RangeTree_ms,LayeredRangeTree_ms\n";
    
    // Group by points and dimension
    std::map<std::pair<size_t, int>, std::map<std::string, double>> time_data;
    
    for (const auto& result : results) {
        auto key = std::make_pair(result.point_count, result.dimension);
        time_data[key][result.algorithm] = result.construction_time_us / 1000.0;
    }
    
    for (const auto& entry : time_data) {
        size_t points = entry.first.first;
        int dim = entry.first.second;
        const auto& algorithms = entry.second;
        
        file << points << "," << dim << ",";
        file << algorithms.at("KDTree") << ",";
        file << algorithms.at("RangeTree") << ",";
        file << algorithms.at("LayeredRangeTree") << "\n";
    }
}

// Create query time analysis CSV
void create_query_time_csv(const std::vector<BenchmarkResult>& results) {
    std::ofstream file("query_time.csv");
    file << "Points,Dimension,KDTree_us,RangeTree_us,LayeredRangeTree_us\n";
    
    // Group by points and dimension
    std::map<std::pair<size_t, int>, std::map<std::string, double>> query_data;
    
    for (const auto& result : results) {
        auto key = std::make_pair(result.point_count, result.dimension);
        query_data[key][result.algorithm] = result.average_query_time_us;
    }
    
    for (const auto& entry : query_data) {
        size_t points = entry.first.first;
        int dim = entry.first.second;
        const auto& algorithms = entry.second;
        
        file << points << "," << dim << ",";
        file << algorithms.at("KDTree") << ",";
        file << algorithms.at("RangeTree") << ",";
        file << algorithms.at("LayeredRangeTree") << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Default parameters
    DistributionType distribution = DistributionType::UNIFORM;
    size_t point_count = 1000;
    int dimension = 2;
    AlgorithmType algorithm = AlgorithmType::KDTREE;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
        else if (arg == "--distribution" && i + 1 < argc) {
            distribution = string_to_distribution(argv[++i]);
        }
        else if (arg == "--points" && i + 1 < argc) {
            point_count = std::stoul(argv[++i]);
        }
        else if (arg == "--dimension" && i + 1 < argc) {
            dimension = std::stoi(argv[++i]);
            if (dimension < 2 || dimension > 3) {
                std::cerr << "Error: Dimension must be 2 or 3\n";
                return 1;
            }
        }
        else if (arg == "--algorithm" && i + 1 < argc) {
            algorithm = string_to_algorithm(argv[++i]);
        }
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    try {
        if (argc == 1) {
            // No arguments provided, run default test
            run_default_test();
        } else {
            // Run specific benchmark (queries run for 10 seconds)
            auto result = run_benchmark_dispatch(dimension, distribution, point_count, algorithm, 0);
            print_result(result);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}