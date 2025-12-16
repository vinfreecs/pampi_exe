# Set the output format to PNG
set terminal pngcairo size 1600,1200 enhanced font 'Verdana,12'
set output 'performance_scaling_internode.png'

# Set titles and labels
set title "Strong Scaling Performance 4 Node "
set xlabel "MPI Ranks"
set ylabel "Performance (MLUPS/s)"

# Enable grid lines for readability
set grid

# IMPORTANT: Tell gnuplot the file is comma-separated
set datafile separator ","
set yrange [0:*]   # Y-axis starts at 0, top is auto-scaled
set xrange [0:350]   # X-axis starts at 0, right is auto-scaled
# Improve the style of the lines and points
set style line 1 lc rgb '#06216bff' lt 1 lw 2 pt 7 ps 1.5   # Blue, circle points

# Fix the X-axis to show integers only (since ranks are whole numbers)

# Plot the data
plot 'result_bench_internode.csv' using 1:2 every ::1 with linespoints ls 1 title "Measured Performance"