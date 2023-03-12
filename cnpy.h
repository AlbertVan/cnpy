// Copyright (C) 2011  Carl Rogers
// Released under MIT License
// license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#ifndef LIBCNPY_H_
#define LIBCNPY_H_

#include <stdint.h>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

namespace cnpy {

using tp_hash_key = size_t;
using np_tp_tb = std::map<tp_hash_key, std::pair<char, size_t>>;

extern np_tp_tb g_np_tp_tb;

struct NpyArray {
    NpyArray(const std::vector<size_t>& _shape, size_t _word_size, bool _fortran_order,
             char _type_class)
        : shape(_shape),
          word_size(_word_size),
          fortran_order(_fortran_order),
          type_class(_type_class) {
        num_vals = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
        data_holder.reset(new std::vector<char>(num_vals * word_size));
    }

    NpyArray() : shape(0), word_size(0), fortran_order(0), type_class('?'), num_vals(0) {}

    template <typename T>
    T* data() {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template <typename T>
    const T* data() const {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template <typename T>
    std::vector<T> as_vec() const {
        const T* p = data<T>();
        return std::vector<T>(p, p + num_vals);
    }

    size_t num_bytes() const { return data_holder->size(); }

    std::shared_ptr<std::vector<char>> data_holder;
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    char type_class;
    size_t num_vals;
};

using npz_t = std::map<std::string, NpyArray>;

uint32_t crc32(uint32_t crc, const void* data, size_t length);
char map_type(const std::type_info& t);
void parse_npy_header(std::istream& is, size_t& word_size, std::vector<size_t>& shape,
                      char& type_class, bool& fortran_order);
void parse_zip_footer(std::istream& is, uint16_t& nrecs, size_t& global_header_size,
                      size_t& global_header_offset);

std::string create_npy_header(const std::vector<size_t>& shape, char type_class, size_t word_size);
std::string create_local_header(std::string& varname, uint32_t crc, uint32_t nbytes);
std::string create_global_header(std::string& varname, std::string& local, uint32_t offset);
std::string create_footer(uint16_t nrecs, uint32_t gh_size, uint32_t gh_offset);

NpyArray npy_load(const std::string& fname);
npz_t npz_load(const std::string& fname);
NpyArray npz_load(const std::string& fname, const std::string& varname);
npz_t npz_load_buffer(std::string& serilize_data);
std::string npz_save_buffer(const npz_t& arrays);

template <typename T>
void npy_save(std::string fname, const T* data, const std::vector<size_t> shape,
              std::string mode = "w") {
    std::vector<size_t> true_data_shape;  // if appending, the shape of existing + new data
    std::fstream fs;
    if (mode == "a") {
        fs.open(fname, std::ios::in | std::ios::out | std::ios::binary);
    }
    if (fs.is_open()) {
        // file exists. we need to append to it. read the header, modify the array size
        size_t word_size;
        bool fortran_order;
        char type_class;
        parse_npy_header(fs, word_size, true_data_shape, type_class, fortran_order);
        assert(!fortran_order);

        if (word_size != sizeof(T)) {
            std::cout << "libnpy error: " << fname << " has word size " << word_size
                      << " but npy_save appending data sized " << sizeof(T) << "\n";
            assert(word_size == sizeof(T));
        }
        if (true_data_shape.size() != shape.size()) {
            std::cout << "libnpy error: npy_save attempting to append misdimensioned data to "
                      << fname << "\n";
            assert(true_data_shape.size() != shape.size());
        }

        for (size_t i = 1; i < shape.size(); i++) {
            if (shape[i] != true_data_shape[i]) {
                std::cout << "libnpy error: npy_save attempting to append misshaped data to "
                          << fname << "\n";
                assert(shape[i] == true_data_shape[i]);
            }
        }
        true_data_shape[0] += shape[0];
    } else {
        true_data_shape = shape;
        fs.open(fname, std::ios::out | std::ios::binary);
        if (!fs.is_open()) {
            throw std::runtime_error("Can't open file: " + fname + "for write.");
        }
    }
    std::string header = create_npy_header(true_data_shape, map_type(typeid(T)), sizeof(T));
    size_t nels = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    fs.write(&header[0], header.size());
    fs.write((const char*)data, sizeof(T) * nels);
    fs.close();
}

template <typename T>
void npz_save(std::string zipname, std::string fname, const T* data,
              const std::vector<size_t>& shape, std::string mode = "w") {
    // first, append a .npy to the fname
    fname += ".npy";
    // now, on with the show
    uint16_t nrecs = 0;
    size_t global_header_offset = 0;
    std::string global_header;

    std::fstream fs;
    if (mode == "a") {
        fs.open(zipname, std::ios::in | std::ios::out | std::ios::binary);
    }
    if (fs.is_open()) {
        // zip file exists. we need to add a new npy file to it.
        // first read the footer. this gives us the offset and size of the global header
        // then read and store the global header.
        // below, we will write the the new data at the start of the global header then append the
        // global header and footer below it
        size_t global_header_size;
        parse_zip_footer(fs, nrecs, global_header_size, global_header_offset);
        global_header.resize(global_header_size);
        std::cout << "npz_save: global header size = " << global_header_size << std::endl;
        // move get to global_header
        fs.seekg(global_header_offset, std::ios::beg);
        fs.read(&global_header[0], global_header_size);
        // move put to global_header
        fs.seekp(global_header_offset, std::ios::beg);
    } else {
        fs.open(zipname, std::ios::out | std::ios::binary);
        if (!fs.is_open()) {
            throw std::runtime_error("Cannot open " + zipname + " for writing");
        }
    }

    std::string npy_header = create_npy_header(shape, map_type(typeid(T)), sizeof(T));

    size_t nels = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    size_t nbytes = nels * sizeof(T) + npy_header.size();

    // get the CRC of the data to be added
    uint32_t crc = crc32(0L, (void*)&npy_header[0], npy_header.size());
    crc = crc32(crc, (uint8_t*)data, nels * sizeof(T));
    // build the local header
    auto local_header = create_local_header(fname, crc, nbytes);
    // build global header
    auto cur_global_header = create_global_header(fname, local_header, global_header_offset);
    global_header += cur_global_header;
    global_header_offset += (nbytes + local_header.size());
    // build footer
    auto footer = create_footer(nrecs + 1, global_header.size(), global_header_offset);
    // write everything
    fs.write(&local_header[0], local_header.size());
    fs.write(&npy_header[0], npy_header.size());
    fs.write((const char*)data, nels * sizeof(T));
    fs.write(&global_header[0], global_header.size());
    fs.write(&footer[0], footer.size());
    fs.close();
}

template <typename T>
void npy_save(std::string fname, const std::vector<T> data, std::string mode = "w") {
    std::vector<size_t> shape;
    shape.push_back(data.size());
    npy_save(fname, &data[0], shape, mode);
}

template <typename T>
void npz_save(std::string zipname, std::string fname, const std::vector<T> data,
              std::string mode = "w") {
    std::vector<size_t> shape;
    shape.push_back(data.size());
    npz_save(zipname, fname, &data[0], shape, mode);
}
}  // namespace cnpy

#endif