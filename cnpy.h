// Copyright (C) 2011  Carl Rogers
// Released under MIT License
// license available in LICENSE file, or at
// http://www.opensource.org/licenses/mit-license.php

#ifndef LIBCNPY_H_
#define LIBCNPY_H_

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <typeinfo>
#include <vector>

#include <fstream>
#include <sstream>

namespace cnpy {

struct NpyArray {
  NpyArray(const std::vector<size_t> &_shape, size_t _word_size,
           bool _fortran_order, char _type_class)
      : shape(_shape), word_size(_word_size), fortran_order(_fortran_order),
        type_class(_type_class) {
    num_vals = 1;
    for (size_t i = 0; i < shape.size(); i++)
      num_vals *= shape[i];
    data_holder = std::shared_ptr<std::vector<char>>(
        new std::vector<char>(num_vals * word_size));
  }

  NpyArray()
      : shape(0), word_size(0), fortran_order(0), type_class('?'), num_vals(0) {
  }

  template <typename T> T *data() {
    return reinterpret_cast<T *>(&(*data_holder)[0]);
  }

  template <typename T> const T *data() const {
    return reinterpret_cast<T *>(&(*data_holder)[0]);
  }

  template <typename T> std::vector<T> as_vec() const {
    const T *p = data<T>();
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

char BigEndianTest();
char map_type(const std::type_info &t);

uint32_t crc32(uint32_t crc, const void *data, unsigned int length) {
  static const uint32_t tinf_crc32tab[16] = {
      0x00000000, 0x1db71064, 0x3b6e20c8, 0x26d930ac, 0x76dc4190, 0x6b6b51f4,
      0x4db26158, 0x5005713c, 0xedb88320, 0xf00f9344, 0xd6d6a3e8, 0xcb61b38c,
      0x9b64c2b0, 0x86d3d2d4, 0xa00ae278, 0xbdbdf21c};
  const unsigned char *buf = (const unsigned char *)data;
  unsigned int i;
  crc = ~crc;
  for (i = 0; i < length; ++i) {
    crc ^= buf[i];
    crc = tinf_crc32tab[crc & 0x0f] ^ (crc >> 4);
    crc = tinf_crc32tab[crc & 0x0f] ^ (crc >> 4);
  }
  // return value suitable for passing in next time, for final value invert it
  return crc ^ 0xffffffff;
}
std::vector<char> create_npy_header(const std::vector<size_t> &shape,
                                    char type_class, size_t word_size);
// void parse_npy_header(FILE* fp,size_t& word_size, std::vector<size_t>& shape,
// char& type_class, bool& fortran_order); void parse_npy_header(unsigned char*
// buffer,size_t& word_size, std::vector<size_t>& shape, char& type_class, bool&
// fortran_order); void parse_zip_footer(FILE* fp, uint16_t& nrecs, size_t&
// global_header_size, size_t& global_header_offset);
npz_t npz_load(std::string fname);
NpyArray npz_load(std::string fname, std::string varname);
NpyArray npy_load(std::string fname);
// refine
void parse_npy_header(std::istream &is, size_t &word_size,
                      std::vector<size_t> &shape, char &type_class,
                      bool &fortran_order);
void parse_zip_footer(std::istream &is, uint16_t &nrecs,
                      size_t &global_header_size, size_t &global_header_offset);

npz_t npz_load_buffer(char *buffer, size_t length);
// std::string npz_save_buffer(const npz_t &da_dict) {
//     // std::string buffer;
//     std::stringstream ss_buffer;
//     uint16_t nrecs = 0;
//     size_t global_header_offset = 0;
//     std::vector<char> global_header;
//     for (auto iter = da_dict.begin(); iter != da_dict.end(); ++iter) {
//         std::string var_name = iter->first;
//         NpyArray var_data = iter->second;

//         std::vector<size_t> shape = var_data.shape;
//         char type_class = var_data.type_class;
//         size_t word_size = var_data.word_size;
//         char* data = var_data.data<char>();

//         std::vector<char> npy_header = create_npy_header(shape, type_class,
//         word_size); size_t nels =
//         std::accumulate(shape.begin(),shape.end(),1,std::multiplies<size_t>());
//         size_t nbytes = nels * word_size + npy_header.size();

//         //get the CRC of the data to be added
//         uint32_t crc = crc32(0L,(uint8_t*)&npy_header[0],npy_header.size());
//         crc = crc32(crc, (uint8_t*)data, nels * word_size);

//         //build the local header
//         std::vector<char> local_header;
//         local_header += "PK"; //first part of sig
//         local_header += (uint16_t) 0x0403; //second part of sig
//         local_header += (uint16_t) 20; //min version to extract
//         local_header += (uint16_t) 0; //general purpose bit flag
//         local_header += (uint16_t) 0; //compression method
//         local_header += (uint16_t) 0; //file last mod time
//         local_header += (uint16_t) 0;     //file last mod date
//         local_header += (uint32_t) crc; //crc
//         local_header += (uint32_t) nbytes; //compressed size
//         local_header += (uint32_t) nbytes; //uncompressed size
//         local_header += (uint16_t) var_name.size(); //fname length
//         local_header += (uint16_t) 0; //extra field length
//         local_header += var_name;

//         //build global header
//         global_header += "PK"; //first part of sig
//         global_header += (uint16_t) 0x0201; //second part of sig
//         global_header += (uint16_t) 20; //version made by
//         global_header.insert(global_header.end(), local_header.begin() + 4,
//         local_header.begin() + 30); global_header += (uint16_t) 0; //file
//         comment length global_header += (uint16_t) 0; //disk number where
//         file starts global_header += (uint16_t) 0; //internal file attributes
//         global_header += (uint32_t) 0; //external file attributes
//         global_header += (uint32_t) global_header_offset; //relative offset
//         of local file header, since it begins where the global header used to
//         begin global_header += var_name;

//         //build footer
//         std::vector<char> footer;
//         footer += "PK"; //first part of sig
//         footer += (uint16_t) 0x0605; //second part of sig
//         footer += (uint16_t) 0; //number of this disk
//         footer += (uint16_t) 0; //disk where footer starts
//         footer += (uint16_t) (nrecs+1); //number of records on this disk
//         footer += (uint16_t) (nrecs+1); //total number of records
//         footer += (uint32_t) global_header.size(); //nbytes of global headers
//         footer += (uint32_t) (global_header_offset + nbytes +
//         local_header.size()); //offset of start of global headers, since
//         global header now starts after newly written array footer +=
//         (uint16_t) 0; //zip file comment length

//         //write everything
//         ss_buffer.write(&local_header[0], sizeof(char) *
//         local_header.size()); ss_buffer.write(&npy_header[0], sizeof(char) *
//         npy_header.size()); ss_buffer.write(data, word_size * nels);
//         ss_buffer.write(&global_header[0], sizeof(char) *
//         global_header.size()); ss_buffer.write(&footer[0], sizeof(char) *
//         footer.size());
//         // fwrite(&local_header[0],sizeof(char),local_header.size(),fp);
//         // fwrite(&npy_header[0],sizeof(char),npy_header.size(),fp);
//         // fwrite(data,sizeof(T),nels,fp);
//         // fwrite(&global_header[0],sizeof(char),global_header.size(),fp);
//         // fwrite(&footer[0],sizeof(char),footer.size(),fp);
//     }
//     return ss_buffer.str();
// }

template <typename T>
std::vector<char> &operator+=(std::vector<char> &lhs, const T rhs) {
  // write in little endian
  for (size_t byte = 0; byte < sizeof(T); byte++) {
    char val = *((char *)&rhs + byte);
    lhs.push_back(val);
  }
  return lhs;
}

template <>
std::vector<char> &operator+=(std::vector<char> &lhs, std::string &rhs);
template <>
std::vector<char> &operator+=(std::vector<char> &lhs, const std::string rhs);
template <>
std::vector<char> &operator+=(std::vector<char> &lhs, const char *rhs);

template <typename T>
std::vector<char> create_npy_header(const std::vector<size_t> &shape) {
  char type_class = map_type(typeid(T));
  size_t word_size = sizeof(T);
  return create_npy_header(shape, type_class, word_size);
}

template <typename T>
void npy_save(std::string fname, const T *data, const std::vector<size_t> shape,
              std::string mode = "w") {
  // FILE* fp = NULL;
  std::vector<size_t>
      true_data_shape; // if appending, the shape of existing + new data

  // if(mode == "a") fp = fopen(fname.c_str(),"r+b");
  std::fstream fs;
  if (mode == "a") {
    fs.open(fname, std::ios::in | std::ios::binary);
  }
  if (fs.is_open()) {
    // file exists. we need to append to it. read the header, modify the array
    // size
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
      std::cout << "libnpy error: npy_save attempting to append misdimensioned "
                   "data to "
                << fname << "\n";
      assert(true_data_shape.size() != shape.size());
    }

    for (size_t i = 1; i < shape.size(); i++) {
      if (shape[i] != true_data_shape[i]) {
        std::cout
            << "libnpy error: npy_save attempting to append misshaped data to "
            << fname << "\n";
        assert(shape[i] == true_data_shape[i]);
      }
    }
    true_data_shape[0] += shape[0];
  } else {
    // fp = fopen(fname.c_str(),"wb");
    std::cout << "open npy file as mode append failed." << std::endl;
    fs.open(fname, std::ios::out | std::ios::binary);
    true_data_shape = shape;
  }

  std::vector<char> header = create_npy_header<T>(true_data_shape);
  size_t nels =
      std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());

  // fseek(fp,0,SEEK_SET);
  // fwrite(&header[0],sizeof(char),header.size(),fp);
  // fseek(fp,0,SEEK_END);
  // fwrite(data,sizeof(T),nels,fp);
  // fclose(fp);

  // fs.seekp(0);
  fs.write(&header[0], header.size());
  // fs.seekp(0, std::ios::end);
  fs.write((const char *)data, sizeof(T) * nels);
  fs.close();
}

std::vector<char> make_local_header(std::string varname, uint32_t crc,
                                    uint32_t nbytes) {
  const uint16_t sig = 0x0403;
  const uint16_t ver = 20;
  const uint16_t null_val = 0;
  uint16_t varname_size = varname.size();
  std::stringstream ss;
  ss << "PK";
  ss.write((char *)&sig, 2);
  ss.write((char *)&ver, 2);
  ss.write((char *)&null_val, 2);
  ss.write((char *)&null_val, 2);
  ss.write((char *)&null_val, 2);
  ss.write((char *)&null_val, 2);
  ss.write((char *)&crc, 4);
  ss.write((char *)&nbytes, 4);
  ss.write((char *)&nbytes, 4);
  ss.write((char *)&varname_size, 2);
  ss.write((char *)&null_val, 2);
  std::string temp = ss.str();
  return std::vector<char>(temp.begin(), temp.end());
}

std::vector<char> make_local_header_old(std::string varname, uint32_t crc,
                                        uint32_t nbytes) {
  std::vector<char> local_header;
  local_header += "PK";                     // first part of sig
  local_header += (uint16_t)0x0403;         // second part of sig
  local_header += (uint16_t)20;             // min version to extract
  local_header += (uint16_t)0;              // general purpose bit flag
  local_header += (uint16_t)0;              // compression method
  local_header += (uint16_t)0;              // file last mod time
  local_header += (uint16_t)0;              // file last mod date
  local_header += (uint32_t)crc;            // crc
  local_header += (uint32_t)nbytes;         // compressed size
  local_header += (uint32_t)nbytes;         // uncompressed size
  local_header += (uint16_t)varname.size(); // fname length
  local_header += (uint16_t)0;              // extra field length
  local_header += varname;
  return local_header;
}

template <typename T>
void npz_save(std::string zipname, std::string fname, const T *data,
              const std::vector<size_t> &shape, std::string mode = "w") {
  std::cout << "=== begine to save var: " << fname << std::endl;
  // first, append a .npy to the fname
  fname += ".npy";

  // now, on with the show
  //  FILE* fp = NULL;
  uint16_t nrecs = 0;
  size_t global_header_offset = 0;
  std::vector<char> global_header;

  // if(mode == "a") fp = fopen(zipname.c_str(),"r+b");
  std::fstream fs;
  if (mode == "a") {
    fs.open(zipname, std::ios::in | std::ios::out | std::ios::binary);
  }
  if (fs.is_open()) {
    // if(fp) {
    // zip file exists. we need to add a new npy file to it.
    // first read the footer. this gives us the offset and size of the global
    // header then read and store the global header. below, we will write the
    // the new data at the start of the global header then append the global
    // header and footer below it
    size_t global_header_size;
    parse_zip_footer(fs, nrecs, global_header_size, global_header_offset);
    // fseek(fp,global_header_offset,SEEK_SET);
    fs.seekg(global_header_offset);
    global_header.resize(global_header_size);

    std::cout << "global header size = " << global_header_size << std::endl;
    // size_t res = fread(&global_header[0],sizeof(char),global_header_size,fp);
    // if(res != global_header_size){
    //     throw std::runtime_error("npz_save: header read error while adding to
    //     existing zip");
    // }
    fs.read(&global_header[0], global_header_size);
    fs.seekp(global_header_offset);

    std::cout << "global header size = " << global_header_size << std::endl;
    // fseek(fp,global_header_offset,SEEK_SET);
    std::cout << "cur seek " << fs.tellg() << std::endl;
    std::cout << "cur seek p " << fs.tellp() << std::endl;
  } else {
    // fp = fopen(zipname.c_str(),"wb");
    std::cout << "npz save load file name as mode 'a' failed " << std::endl;
    fs.open(zipname, std::ios::out | std::ios::binary);
  }

  std::vector<char> npy_header = create_npy_header<T>(shape);

  size_t nels =
      std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
  size_t nbytes = nels * sizeof(T) + npy_header.size();

  // get the CRC of the data to be added
  uint32_t crc = crc32(0L, (uint8_t *)&npy_header[0], npy_header.size());
  crc = crc32(crc, (uint8_t *)data, nels * sizeof(T));

  std::cout << "crc32 : " << crc << std::endl;

  // build the local header
  std::vector<char> local_header;
  local_header += "PK";                   // first part of sig
  local_header += (uint16_t)0x0403;       // second part of sig
  local_header += (uint16_t)20;           // min version to extract
  local_header += (uint16_t)0;            // general purpose bit flag
  local_header += (uint16_t)0;            // compression method
  local_header += (uint16_t)0;            // file last mod time
  local_header += (uint16_t)0;            // file last mod date
  local_header += (uint32_t)crc;          // crc
  local_header += (uint32_t)nbytes;       // compressed size
  local_header += (uint32_t)nbytes;       // uncompressed size
  local_header += (uint16_t)fname.size(); // fname length
  local_header += (uint16_t)0;            // extra field length
  local_header += fname;

  std::cout << std::hex << "local_header: nbytes:" << nbytes
            << " , fname size: " << fname.size() << " fname: " << fname
            << std::endl;

  // build global header
  global_header += "PK";             // first part of sig
  global_header += (uint16_t)0x0201; // second part of sig
  global_header += (uint16_t)20;     // version made by
  global_header.insert(global_header.end(), local_header.begin() + 4,
                       local_header.begin() + 30);
  global_header += (uint16_t)0; // file comment length
  global_header += (uint16_t)0; // disk number where file starts
  global_header += (uint16_t)0; // internal file attributes
  global_header += (uint32_t)0; // external file attributes
  global_header += (uint32_t)
      global_header_offset; // relative offset of local file header, since it
                            // begins where the global header used to begin
  global_header += fname;
  std::cout << "global_header: " << global_header_offset << std::endl;

  // build footer
  std::vector<char> footer;
  footer += "PK";                           // first part of sig
  footer += (uint16_t)0x0605;               // second part of sig
  footer += (uint16_t)0;                    // number of this disk
  footer += (uint16_t)0;                    // disk where footer starts
  footer += (uint16_t)(nrecs + 1);          // number of records on this disk
  footer += (uint16_t)(nrecs + 1);          // total number of records
  footer += (uint32_t)global_header.size(); // nbytes of global headers
  footer += (uint32_t)(global_header_offset + nbytes +
                       local_header.size()); // offset of start of global
                                             // headers, since global header now
                                             // starts after newly written array
  footer += (uint16_t)0; // zip file comment length
  std::cout << "footer: nrecs:" << nrecs + 1
            << " global_header_size: " << global_header.size() << " all size:"
            << (global_header_offset + nbytes + local_header.size()) << std::dec
            << std::endl;
  // write everything
  //  fwrite(&local_header[0],sizeof(char),local_header.size(),fp);
  //  fwrite(&npy_header[0],sizeof(char),npy_header.size(),fp);
  //  fwrite(data,sizeof(T),nels,fp);
  //  fwrite(&global_header[0],sizeof(char),global_header.size(),fp);
  //  fwrite(&footer[0],sizeof(char),footer.size(),fp);
  //  fclose(fp);

  fs.write(&local_header[0], local_header.size());
  fs.write(&npy_header[0], npy_header.size());
  fs.write((const char *)data, nels * sizeof(T));
  fs.write(&global_header[0], global_header.size());
  fs.write(&footer[0], footer.size());

  fs.close();
}

template <typename T>
void npy_save(std::string fname, const std::vector<T> data,
              std::string mode = "w") {
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

} // namespace cnpy

#endif
