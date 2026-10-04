rm -rf build
rm -rf CMakeFiles
rm -f CMakeCache.txt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel "$(nproc)"
