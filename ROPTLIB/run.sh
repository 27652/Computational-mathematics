rm -rf build
rm -rf CMakeFiles
rm -f CMakeCache.txt
cmake -S . -B build
cmake --build build --parallel "$(nproc)"
