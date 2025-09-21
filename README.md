# Range Search Algorithms

This project implements and compares different range searching algorithms in C++ with both interactive visualization and benchmarking capabilities.

## Features

- **Interactive Visualization**: Qt-based GUI for visualizing range search operations
- **Performance Benchmarking**: Command-line tool for measuring algorithm performance
- **KD-Tree Implementation**: Efficient range searching with CGAL integration

## Structure
- `src/`: Source files for algorithm implementations
  - `visualization_main.cpp`: Qt GUI application entry point
  - `benchmark_main.cpp`: Benchmark tool entry point
  - `visualization.cpp/h`: Qt visualization components
- `include/`: Header files for algorithms and geometry
- `docs/`: Documentation
- `benchmark_results/`: Generated benchmark data and plots

## Build

```bash
mkdir build
cd build
cmake ..
make
```

This creates two executables:
- `range_search_viz`: Interactive visualization (requires Qt5)
- `range_search_benchmark`: Performance benchmarking tool

## Usage

### Visualization Tool
```bash
./range_search_viz
```
- Interactive 2D point visualization
- Draggable query rectangles
- Real-time range search results
- Optional KD-tree partition display

### Benchmark Tool
```bash
# Run benchmarks for different algorithms and dimensions
./range_search_benchmark
```

The benchmark tool tests:
- Construction time performance
- Query time performance  
- Memory usage analysis
- Multiple dimensions (2D to 6D)
- Different data distributions (uniform, gaussian, clustered)
- KD-Tree and Range Tree algorithms

Some results are saved to CSV files in the `benchmark_results/` directory.

## Dependencies

- **Required**: CGAL, Eigen3, CMake 3.10+
- **Optional**: Qt5 (for visualization)

## Algorithms
- KD-Tree
- Range Tree
- Layered Range Tree
