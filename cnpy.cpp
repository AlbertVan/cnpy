// Copyright (C) 2011  Carl Rogers
// Released under MIT License
// license available in LICENSE file, or at http://www.opensource.org/licenses/mit-license.php

#include "cnpy.h"
#include <stdint.h>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <regex>
#include <stdexcept>

namespace cnpy {

char big_endian_test() {
    int x = 1;
    return (((char*)&x)[0]) ? '<' : '>';
}

char map_type(const std::type_info& t) {
    if (t == typeid(float)) {
        return 'f';
    }
    if (t == typeid(double)) {
        return 'f';
    }
    if (t == typeid(long double)) {
        return 'f';
    }

    if (t == typeid(int)) {
        return 'i';
    }
    if (t == typeid(char)) {
        return 'i';
    }
    if (t == typeid(short)) {
        return 'i';
    }
    if (t == typeid(long)) {
        return 'i';
    }
    if (t == typeid(long long)) {
        return 'i';
    }

    if (t == typeid(unsigned char)) {
        return 'u';
    }
    if (t == typeid(unsigned short)) {
        return 'u';
    }
    if (t == typeid(unsigned long)) {
        return 'u';
    }
    if (t == typeid(unsigned long long)) {
        return 'u';
    }
    if (t == typeid(unsigned int)) {
        return 'u';
    }

    if (t == typeid(bool)) {
        return 'b';
    }

    if (t == typeid(std::complex<float>)) {
        return 'c';
    }
    if (t == typeid(std::complex<double>)) {
        return 'c';
    }
    if (t == typeid(std::complex<long double>)) {
        return 'c';
    }

    else
        return '?';
}

uint32_t crc32(uint32_t crc, const void* data, size_t length) {
    static const uint32_t tinf_crc32tab[16] = {0x00000000, 0x1db71064, 0x3b6e20c8, 0x26d930ac,
                                               0x76dc4190, 0x6b6b51f4, 0x4db26158, 0x5005713c,
                                               0xedb88320, 0xf00f9344, 0xd6d6a3e8, 0xcb61b38c,
                                               0x9b64c2b0, 0x86d3d2d4, 0xa00ae278, 0xbdbdf21c};
    const unsigned char* buf = (const unsigned char*)data;
    crc = ~crc;
    for (size_t i = 0; i < length; ++i) {
        crc ^= buf[i];
        crc = tinf_crc32tab[crc & 0x0f] ^ (crc >> 4);
        crc = tinf_crc32tab[crc & 0x0f] ^ (crc >> 4);
    }
    // return value suitable for passing in next time, for final value invert it
    return crc ^ 0xffffffff;
}

void parse_npy_header(std::istream& is, size_t& word_size, std::vector<size_t>& shape,
                      char& type_class, bool& fortran_order) {
    // get line first;
    std::string header = "";
    std::getline(is, header, '\n');

    size_t loc1 = 0, loc2 = 0;
    // fortran order
    loc1 = header.find("fortran_order");
    if (loc1 == std::string::npos) {
        throw std::runtime_error(
            "parse_npy_header: failed to find header keyword: 'fortran_order'");
    }
    loc1 += 16;
    fortran_order = (header.substr(loc1, 4) == "True" ? true : false);

    // shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    if (loc1 == std::string::npos || loc2 == std::string::npos) {
        throw std::runtime_error("parse_npy_header: failed to find header keyword: '(' or ')'");
    }

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1 + 1, loc2 - loc1 - 1);
    while (std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    // endian, word size, data type
    // byte order code | stands for not applicable.
    // not sure when this applies except for byte array
    loc1 = header.find("descr");
    if (loc1 == std::string::npos) {
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'descr'");
    }
    loc1 += 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    // char type = header[loc1+1];
    // assert(type == map_type(T));

    std::string str_ws = header.substr(loc1 + 1);
    type_class = str_ws.at(0);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(1, loc2).c_str());
}

void parse_zip_footer(std::istream& is, uint16_t& nrecs, size_t& global_header_size,
                      size_t& global_header_offset) {
    std::vector<char> footer(22);
    is.seekg(-22, std::ios::end);
    is.read(&footer[0], 22);
    uint16_t disk_no, disk_start, nrecs_on_disk, comment_len;
    disk_no = *(uint16_t*)&footer[4];
    disk_start = *(uint16_t*)&footer[6];
    nrecs_on_disk = *(uint16_t*)&footer[8];
    nrecs = *(uint16_t*)&footer[10];
    global_header_size = *(uint32_t*)&footer[12];
    global_header_offset = *(uint32_t*)&footer[16];
    comment_len = *(uint16_t*)&footer[20];

    assert(disk_no == 0);
    assert(disk_start == 0);
    assert(nrecs_on_disk == nrecs);
    assert(comment_len == 0);
}

std::string parse_local_header(std::istream& is, uint16_t& compr_method, uint32_t& compr_bytes,
                               uint32_t& uncompr_bytes) {
    std::vector<char> local_header(30);
    is.read(&local_header[0], 30);
    // if we've reached the global header, stop reading
    if (local_header[2] != 0x03 || local_header[3] != 0x04) {
        return "";
    }
    // read in the variable name
    uint16_t name_len = *(uint16_t*)&local_header[26];
    std::string var_name(name_len, ' ');
    is.read(&var_name[0], name_len);
    var_name.erase(var_name.end() - 4, var_name.end());  // erase ".npy"
    // read in the extra field
    uint16_t extra_field_len = *(uint16_t*)&local_header[28];
    if (extra_field_len > 0) {
        std::vector<char> buff(extra_field_len);
        is.read(&buff[0], extra_field_len * sizeof(char));
    }

    compr_method = *reinterpret_cast<uint16_t*>(&local_header[8]);
    compr_bytes = *reinterpret_cast<uint32_t*>(&local_header[18]);
    uncompr_bytes = *reinterpret_cast<uint32_t*>(&local_header[22]);
    if (compr_method != 0) {
        throw std::runtime_error("npz_load: unsupported gzip compress's npz.");
        // ss.read(buffer_compr.data(), compr_bytes);
        // arrays[varname] = load_the_npz_array((unsigned
        // char*)buffer_compr.data(),compr_bytes,uncompr_bytes);
    }
    return var_name;
}

NpyArray load_the_npy_stream(std::istream& is) {
    std::vector<size_t> shape{};
    size_t word_size;
    bool fortran_order;
    char type_class;
    parse_npy_header(is, word_size, shape, type_class, fortran_order);
    NpyArray arr(shape, word_size, fortran_order, type_class);
    is.read(arr.data<char>(), arr.num_bytes());
    return arr;
}

npz_t npz_load_buffer(std::string& serialize_data) {
    npz_t arrays;
    std::stringstream ss(serialize_data);
    while (1) {
        uint16_t compr_method = 0;
        uint32_t compr_bytes = 0, uncompr_bytes = 0;
        std::string var_name = parse_local_header(ss, compr_method, compr_bytes, uncompr_bytes);
        // if var_name == "", stop reading.
        if (var_name == "") {
            break;
        }
        arrays[var_name] = load_the_npy_stream(ss);
    }
    return arrays;
}

npz_t npz_load(const std::string& fname) {
    std::ifstream ifs;
    ifs.open(fname, std::ios::binary | std::ios::in);
    if (!ifs.is_open()) {
        throw std::runtime_error("npz_load: Unable to open file " + fname);
    }
    npz_t arrays;
    while (1) {
        uint16_t compr_method = 0;
        uint32_t compr_bytes = 0, uncompr_bytes = 0;
        std::string var_name = parse_local_header(ifs, compr_method, compr_bytes, uncompr_bytes);
        // if var_name == "", stop reading.
        if (var_name == "") {
            break;
        }
        arrays[var_name] = load_the_npy_stream(ifs);
    }
    ifs.close();
    return arrays;
}

NpyArray npz_load(const std::string& fname, const std::string& varname) {
    std::ifstream ifs;
    ifs.open(fname, std::ios::binary | std::ios::in);
    if (!ifs.is_open()) {
        throw std::runtime_error("npz_load: Unable to open file " + fname);
    }
    while (1) {
        uint16_t compr_method = 0;
        uint32_t compr_bytes = 0, uncompr_bytes = 0;
        std::string var_name = parse_local_header(ifs, compr_method, compr_bytes, uncompr_bytes);
        // if var_name == "", stop reading.
        if (var_name == "") {
            break;
        }
        if (var_name != varname) {
            ifs.seekg(uncompr_bytes, std::ios::cur);
            continue;
        }
        NpyArray array = load_the_npy_stream(ifs);
        ifs.close();
        return array;
    }
    ifs.close();
    // if we get here, we haven't found the variable in the file
    throw std::runtime_error("npz_load: Variable name " + varname + " not found in " + fname);
}

NpyArray npy_load(const std::string& fname) {
    std::ifstream ifs;
    ifs.open(fname, std::ios::binary | std::ios::in);
    if (!ifs.is_open()) {
        throw std::runtime_error("npy_load: Unable to open file " + fname);
    }
    NpyArray arr = load_the_npy_stream(ifs);
    ifs.close();
    return arr;
}

std::string create_npy_header(const std::vector<size_t>& shape, char type_class, size_t word_size) {
    std::string dict;
    dict += "{'descr': '";
    dict += big_endian_test();
    dict += type_class;
    dict += std::to_string(word_size);
    dict += "', 'fortran_order': False, 'shape': (";
    dict += std::to_string(shape[0]);
    for (size_t i = 1; i < shape.size(); ++i) {
        dict += ", ";
        dict += std::to_string(shape[i]);
    }
    if (shape.size() == 1) {
        dict += ",";
    }
    dict += "), }";
    int remainder = 16 - (10 + dict.size()) % 16;
    std::string remainder_str(remainder - 1, ' ');
    dict += remainder_str;
    dict += "\n";

    std::stringstream header;
    header << (char)0x93 << "NUMPY" << (char)0x01 << (char)0x00;
    uint16_t dict_size = (uint16_t)dict.size();
    header.write((char*)&dict_size, 2);
    header << dict;
    return header.str();
}

std::string create_local_header(std::string& varname, uint32_t crc, uint32_t nbytes) {
    const uint16_t local_sig = 0x0403;   // second part of sig
    const uint16_t min_version = 20;     // min version to extract
    const uint16_t bit_flag = 0;         // general purpose bit flag
    const uint16_t compress_method = 0;  // compression method
    const uint16_t last_mod_time = 0;    // file last mod time
    const uint16_t last_mod_date = 0;    // file last mod date
    const uint16_t extra_length = 0;     // exta field length
    uint16_t var_name_size = (uint16_t)varname.size();
    std::stringstream ss;
    ss << "PK";
    ss.write((char*)&local_sig, 2);
    ss.write((char*)&min_version, 2);
    ss.write((char*)&bit_flag, 2);
    ss.write((char*)&compress_method, 2);
    ss.write((char*)&last_mod_time, 2);
    ss.write((char*)&last_mod_date, 2);
    ss.write((char*)&crc, 4);
    ss.write((char*)&nbytes, 4);  // compr bytes
    ss.write((char*)&nbytes, 4);  // uncompr bytes
    ss.write((char*)&var_name_size, 2);
    ss.write((char*)&extra_length, 2);
    ss << varname;
    return ss.str();
}

std::string create_global_header(std::string& varname, std::string& local, uint32_t offset) {
    const uint16_t global_sig = 0x0201;  // second part of sig
    const uint16_t min_version = 20;     // version made by
    const uint16_t comment_len = 0;      // file comment length
    const uint16_t disk_num = 0;         // disk number where file starts
    const uint16_t inter_file_attr = 0;  // internal file attributes
    const uint32_t exter_file_attr = 0;  // external file attributes
    std::stringstream ss;
    ss << "PK";
    ss.write((char*)&global_sig, 2);
    ss.write((char*)&min_version, 2);
    ss.write(&local[4], 26);
    ss.write((char*)&comment_len, 2);
    ss.write((char*)&disk_num, 2);
    ss.write((char*)&inter_file_attr, 2);
    ss.write((char*)&exter_file_attr, 4);
    ss.write((char*)&offset, 4);
    ss << varname;
    return ss.str();
}

std::string create_footer(uint16_t nrecs, uint32_t gh_size, uint32_t gh_offset) {
    const uint16_t footer_sig = 0x0605;  // second part of sig
    const uint16_t num_disk = 0;         // number of this disk
    const uint16_t start_disk = 0;       // disk where footer starts
    const uint16_t zip_comment_len = 0;  // zip file comment length
    std::stringstream ss;
    ss << "PK";
    ss.write((char*)&footer_sig, 2);
    ss.write((char*)&num_disk, 2);
    ss.write((char*)&start_disk, 2);
    ss.write((char*)&nrecs, 2);
    ss.write((char*)&nrecs, 2);
    ss.write((char*)&gh_size, 4);
    ss.write((char*)&gh_offset, 4);
    ss.write((char*)&zip_comment_len, 2);
    return ss.str();
}

std::string npz_save_buffer(const npz_t& arrays) {
    std::stringstream ss_buffer;
    size_t global_header_offset = 0;
    std::string global_header;
    std::string var_name = "";
    uint32_t footer_offset = 0;
    for (auto iter = arrays.begin(); iter != arrays.end(); ++iter) {
        var_name = iter->first + ".npy";
        NpyArray var_data = iter->second;

        std::vector<size_t> shape = var_data.shape;
        char type_class = var_data.type_class;
        size_t word_size = var_data.word_size;
        char* data = var_data.data<char>();

        std::string npy_header = create_npy_header(shape, type_class, word_size);
        size_t nels = std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
        size_t nbytes = nels * word_size + npy_header.size();

        // get the CRC of the data to be added
        uint32_t crc = crc32(0L, &npy_header[0], npy_header.size());
        crc = crc32(crc, (uint8_t*)data, nels * word_size);

        // build the local header
        auto local_header = create_local_header(var_name, crc, nbytes);
        // build global header
        auto cur_global_header = create_global_header(var_name, local_header, global_header_offset);
        global_header += cur_global_header;
        global_header_offset += (nbytes + local_header.size());
        // write everything
        ss_buffer.write(&local_header[0], sizeof(char) * local_header.size());
        ss_buffer.write(&npy_header[0], sizeof(char) * npy_header.size());
        ss_buffer.write(data, word_size * nels);
    }
    // build footer
    auto footer = create_footer(arrays.size(), global_header.size(), global_header_offset);
    ss_buffer.write(&global_header[0], sizeof(char) * global_header.size());
    ss_buffer.write(&footer[0], sizeof(char) * footer.size());
    return ss_buffer.str();
}

}  // namespace cnpy