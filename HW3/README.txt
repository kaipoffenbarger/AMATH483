g++ -std=c++17 -c -fPIC ref_daxpy.cpp
g++ -std=c++17 -c -fPIC ref_axpyt.cpp
g++ -std=c++17 -c -fPIC ref_gemmt.cpp
g++ -std=c++17 -c -fPIC ref_dgemm.cpp
g++ -std=c++17 -c -fPIC ref_dgemv.cpp
g++ -std=c++17 -c -fPIC ref_gemvt.cpp
g++ -shared -o librefBLAS.so ref_daxpy.o ref_dgemv.o ref_dgemm.o ref_axpyt.o ref_gemvt.o ref_gemmt.o