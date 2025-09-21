#!/bin/bash

# Generate 3 essential benchmark graphs directly from main CSV files

echo "Generating 3 essential graphs..."

# Memory Usage Graph
gnuplot <<EOF
set terminal png size 1200,800
set output 'memory_usage.png'
set title 'Memory Usage vs Number of Points (2D)'
set xlabel 'Number of Points'
set ylabel 'Memory Usage (KB)'
set logscale xy
set grid
set datafile separator ','
set key top left
plot 'memory_analysis.csv' using 1:3 with linespoints linewidth 2 pointtype 7 title 'KDTree', \
     'memory_analysis.csv' using 1:4 with linespoints linewidth 2 pointtype 5 title 'RangeTree', \
     'memory_analysis.csv' using 1:5 with linespoints linewidth 2 pointtype 9 title 'LayeredRangeTree'
EOF

# Construction Time Graph  
gnuplot <<EOF
set terminal png size 1200,800
set output 'construction_time.png'
set title 'Construction Time vs Number of Points (2D)'
set xlabel 'Number of Points'
set ylabel 'Construction Time (ms)'
set logscale xy
set grid
set datafile separator ','
set key top left
plot 'construction_time.csv' using 1:3 with linespoints linewidth 2 pointtype 7 title 'KDTree', \
     'construction_time.csv' using 1:4 with linespoints linewidth 2 pointtype 5 title 'RangeTree', \
     'construction_time.csv' using 1:5 with linespoints linewidth 2 pointtype 9 title 'LayeredRangeTree'
EOF

# Query Time Graph
gnuplot <<EOF
set terminal png size 1200,800
set output 'query_time.png'
set title 'Query Time vs Number of Points (2D)'
set xlabel 'Number of Points'
set ylabel 'Average Query Time (μs)'
#set logscale xy
set grid
set datafile separator ','
set key top left
plot 'query_time.csv' using 1:3 with linespoints linewidth 2 pointtype 7 title 'KDTree', \
     'query_time.csv' using 1:4 with linespoints linewidth 2 pointtype 5 title 'RangeTree', \
     'query_time.csv' using 1:5 with linespoints linewidth 2 pointtype 9 title 'LayeredRangeTree'
EOF

echo "Generated 3 essential graphs:"
echo "  - memory_usage.png"
echo "  - construction_time.png" 
echo "  - query_time.png"
ls -la memory_usage.png construction_time.png query_time.png

echo "Graphs generated successfully!"
ls -la *.png
