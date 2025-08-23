# CGAL Installation Guide

CGAL (Computational Geometry Algorithms Library) is required for this project.

## Linux Installation

1. **Install CGAL and dependencies:**
   ```bash
   sudo apt-get update
   sudo apt-get install libcgal-dev libcgal-qt5-dev
   ```
2. **Verify installation:**
   ```bash
   pkg-config --modversion cgal
   ```

## CMake Integration

The project CMakeLists.txt is configured to find and link CGAL automatically if installed system-wide.

For more details, see: https://www.cgal.org/download.html
