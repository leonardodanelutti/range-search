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
- `tests/`: Unit and performance tests
- `docs/`: Documentation

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
# Run all benchmarks
./range_search_benchmark

# Run specific benchmarks
./range_search_benchmark --construction  # Tree construction only
./range_search_benchmark --queries       # Query performance only
./range_search_benchmark --help          # Show options
```

## Dependencies

- **Required**: CGAL, Eigen3, CMake 3.10+
- **Optional**: Qt5 (for visualization)

## Algorithms
- KD-Tree
- Range Tree

## Contributing
Feel free to open issues or submit pull requests.
