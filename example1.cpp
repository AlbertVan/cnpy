#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "cnpy.h"

const int Nx = 128;
const int Ny = 64;
const int Nz = 32;

static double seconds() {
#if defined(_WIN32) || defined(_WIN64)
    LARGE_INTEGER frequency, now;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&now);
    return now.QuadPart / double(frequency.QuadPart);
#else
    timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return now.tv_sec + now.tv_nsec / 1000000000.0;
#endif
}

void test_type_id() {
    std::cout << "bool   : " << typeid(bool).name() << std::endl;
    std::cout << "char   : " << typeid(char).name() << std::endl;
    std::cout << "short  : " << typeid(short).name() << std::endl;
    std::cout << "int    : " << typeid(int).name() << std::endl;
    std::cout << "long   : " << typeid(long).name() << std::endl;
    std::cout << "float  : " << typeid(float).name() << std::endl;
    std::cout << "double : " << typeid(double).name() << std::endl;
    std::cout << "complex: " << typeid(std::complex<float>).name() << std::endl;
    bool flag = true;
    int int_flag = 2;
    long long_flag = 2L;
    std::cout << std::hex;
    std::cout << "bool   : " << typeid(flag).hash_code() << std::endl;
    std::cout << "bool   : " << typeid(bool).hash_code() << std::endl;
    std::cout << "char   : " << typeid(char).hash_code() << std::endl;
    std::cout << "short  : " << typeid(short).hash_code() << std::endl;
    std::cout << "int    : " << typeid(int).hash_code() << std::endl;
    std::cout << "int flag   : " << typeid(int_flag).hash_code() << std::endl;
    std::cout << "long  : " << typeid(long).hash_code() << std::endl;
    std::cout << "long flag : " << typeid(long_flag).hash_code() << std::endl;
    std::cout << "float  : " << typeid(float).hash_code() << std::endl;
    std::cout << "double : " << typeid(double).hash_code() << std::endl;
    std::cout << "complex: " << typeid(std::complex<float>).hash_code() << std::endl;
    std::cout << "complex: " << typeid(std::complex<double>).hash_code() << std::endl;
    std::cout << std::dec;
}

int main() {
    test_type_id();
    double startTime, duration;

    startTime = seconds();
    for (int i = 0; i < 200000; ++i) {
        std::string name = "var1";
        auto out = cnpy::create_local_header(name, 0xffffaaaa, 1024);
        // std::cout << "size = " << out.size();
        // std::cout << std::string(out.begin(), out.end()) << std::endl;
    }
    duration = seconds() - startTime;
    printf("new      :  %.3fs, \n", duration);

    startTime = seconds();
    for (int i = 0; i < 200000; ++i) {
        std::string name = "var1";
        auto out = cnpy::create_local_header(name, 0xffffaaaa, 1024);
        // std::cout << "size = " << out.size();
        // std::cout << std::string(out.begin(), out.end()) << std::endl;
    }
    duration = seconds() - startTime;
    printf("old      :  %.3fs, \n", duration);

    // set random seed so that result is reproducible (for testing)
    srand(0);
    // create random data
    std::vector<std::complex<double>> data(Nx * Ny * Nz);
    for (int i = 0; i < Nx * Ny * Nz; i++) data[i] = std::complex<double>(rand(), rand());

    // save it to file
    cnpy::npy_save("arr1.npy", &data[0], {Nz, Ny, Nx}, "w");

    // load it into a new array
    cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    std::complex<double>* loaded_data = arr.data<std::complex<double>>();

    // make sure the loaded data matches the saved data
    assert(arr.word_size == sizeof(std::complex<double>));
    assert(arr.shape.size() == 3 && arr.shape[0] == Nz && arr.shape[1] == Ny && arr.shape[2] == Nx);
    for (int i = 0; i < Nx * Ny * Nz; i++) assert(data[i] == loaded_data[i]);

    std::cout << "npy load success " << std::endl;

    // append the same data to file
    // npy array on file now has shape (Nz+Nz,Ny,Nx)
    cnpy::npy_save("arr2.npy", &data[0], {Nz, Ny, Nx}, "a");
    // crc32()
    // uint32_t self_crc =
    //     uzlib_crc32((const void *)&data[0], arr.num_bytes(), 0);

    // //   self_crc = self_crc ^ 0xffffffff;
    // uint32_t zlib_crc =
    //     crc32(0, (const unsigned char *)&data[0], arr.num_bytes());
    // std::cout << "mycrc " << std::hex << self_crc << std::endl;
    // std::cout << "zlibcrc " << zlib_crc << std::dec << std::endl;

    // now write to an npz file
    // non-array variables are treated as 1D arrays with 1 element
    double myVar1 = 1.2;
    char myVar2 = 'a';
    cnpy::npz_save("out.npz", "myVar1", &myVar1, {1}, "w");  //"w" overwrites any existing file
    cnpy::npz_save("out.npz", "myVar2", &myVar2, {1},
                   "a");  //"a" appends to the file we created above
    cnpy::npz_save("out.npz", "arr1", &data[0], {Nz, Ny, Nx},
                   "a");  //"a" appends to the file we created above

    // load a single var from the npz file
    cnpy::NpyArray arr2 = cnpy::npz_load("out.npz", "arr1");

    // load the entire npz file
    cnpy::npz_t my_npz = cnpy::npz_load("out.npz");

    // check that the loaded myVar1 matches myVar1
    cnpy::NpyArray arr_mv1 = my_npz["myVar1"];
    double* mv1 = arr_mv1.data<double>();
    assert(arr_mv1.shape.size() == 1 && arr_mv1.shape[0] == 1);
    assert(mv1[0] == myVar1);
}
