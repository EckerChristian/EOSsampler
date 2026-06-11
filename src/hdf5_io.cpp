#include "hdf5_io.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace eos::hdf5 {
namespace {

void check_id(hid_t id, const std::string& what) {
    if (id < 0) throw std::runtime_error("HDF5 error while creating/opening " + what);
}

void check_status(herr_t status, const std::string& what) {
    if (status < 0) throw std::runtime_error("HDF5 error while " + what);
}

std::string rank_file_path(const std::string& output_prefix, int rank, int file_start_index) {
    const int file_index = file_start_index + rank;
    return output_prefix + std::to_string(file_index) + ".h5";
}

void create_dataset(hid_t group_id, const std::string& name, hid_t type, hsize_t length) {
    hsize_t dims[1] = {length};
    hid_t space_id = H5Screate_simple(1, dims, nullptr);
    check_id(space_id, "dataspace " + name);

    hid_t dataset_id = H5Dcreate2(
        group_id,
        name.c_str(),
        type,
        space_id,
        H5P_DEFAULT,
        H5P_DEFAULT,
        H5P_DEFAULT
    );

    H5Sclose(space_id);
    check_id(dataset_id, "dataset " + name);
    H5Dclose(dataset_id);
}

} // namespace

bool link_exists(hid_t group_id, const std::string& name) {
    return H5Lexists(group_id, name.c_str(), H5P_DEFAULT) > 0;
}

hid_t create_file(const std::string& output_prefix, int rank, int file_start_index) {
    const std::string path = rank_file_path(output_prefix, rank, file_start_index);
    hid_t id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    check_id(id, "file " + path);

    return id;
}

hid_t open_file_rw(const std::string& output_prefix, int rank, int file_start_index) {
    const std::string path = rank_file_path(output_prefix, rank, file_start_index);
    hid_t id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    check_id(id, "file " + path);

    return id;
}

hid_t create_group(hid_t file_id, int index) {
    hid_t gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
    check_id(gcpl_id, "group property list");

    check_status(
        H5Pset_link_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED),
        "setting group creation order"
    );

    const std::string name = "/" + std::to_string(index);
    hid_t group_id = H5Gcreate2(file_id, name.c_str(), H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    H5Pclose(gcpl_id);
    check_id(group_id, "group " + name);
    return group_id;
}

hid_t open_group(hid_t file_id, int index) {
    const std::string name = "/" + std::to_string(index);
    hid_t id = H5Gopen2(file_id, name.c_str(), H5P_DEFAULT);
    check_id(id, "group " + name);
    return id;
}

hid_t create_subgroup(hid_t group_id, const std::string& name) {
    hid_t id = H5Gcreate2(group_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check_id(id, "subgroup " + name);
    return id;
}

hid_t open_subgroup(hid_t group_id, const std::string& name) {
    hid_t id = H5Gopen2(group_id, name.c_str(), H5P_DEFAULT);
    check_id(id, "subgroup " + name);
    return id;
}

void create_dataset_double(hid_t group_id, const std::string& name, hsize_t length) {
    create_dataset(group_id, name, H5T_NATIVE_DOUBLE, length);
}

void create_dataset_int(hid_t group_id, const std::string& name, hsize_t length) {
    create_dataset(group_id, name, H5T_NATIVE_INT, length);
}

void recreate_dataset_double(hid_t group_id, const std::string& name, hsize_t length) {
    if (link_exists(group_id, name)) {
        check_status(H5Ldelete(group_id, name.c_str(), H5P_DEFAULT), "deleting dataset " + name);
    }
    create_dataset_double(group_id, name, length);
}

void recreate_dataset_int(hid_t group_id, const std::string& name, hsize_t length) {
    if (link_exists(group_id, name)) {
        check_status(H5Ldelete(group_id, name.c_str(), H5P_DEFAULT), "deleting dataset " + name);
    }
    create_dataset_int(group_id, name, length);
}

std::vector<double> read_double_array(hid_t group_id, const std::string& name) {
    hid_t dataset_id = H5Dopen2(group_id, name.c_str(), H5P_DEFAULT);
    check_id(dataset_id, "dataset " + name);

    hid_t space_id = H5Dget_space(dataset_id);
    check_id(space_id, "dataspace " + name);

    const int ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims != 1) {
        H5Sclose(space_id);
        H5Dclose(dataset_id);
        throw std::runtime_error("Expected one-dimensional dataset: " + name);
    }

    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space_id, dims, nullptr);

    std::vector<double> data(static_cast<std::size_t>(dims[0]));
    check_status(
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()),
        "reading dataset " + name
    );

    H5Sclose(space_id);
    H5Dclose(dataset_id);
    return data;
}

double read_double_scalar(hid_t group_id, const std::string& name) {
    auto data = read_double_array(group_id, name);
    if (data.empty()) throw std::runtime_error("Empty scalar dataset: " + name);
    return data.front();
}

int read_int_scalar(hid_t group_id, const std::string& name) {
    hid_t dataset_id = H5Dopen2(group_id, name.c_str(), H5P_DEFAULT);
    check_id(dataset_id, "dataset " + name);

    int value = 0;
    check_status(
        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value),
        "reading dataset " + name
    );

    H5Dclose(dataset_id);
    return value;
}

void write_double_array(hid_t group_id, const std::string& name, const std::vector<double>& data) {
    hid_t dataset_id = H5Dopen2(group_id, name.c_str(), H5P_DEFAULT);
    check_id(dataset_id, "dataset " + name);

    check_status(
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()),
        "writing dataset " + name
    );

    H5Dclose(dataset_id);
}

void write_double_scalar(hid_t group_id, const std::string& name, double value) {
    std::vector<double> data{value};
    write_double_array(group_id, name, data);
}

void write_int_scalar(hid_t group_id, const std::string& name, int value) {
    hid_t dataset_id = H5Dopen2(group_id, name.c_str(), H5P_DEFAULT);
    check_id(dataset_id, "dataset " + name);

    check_status(
        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value),
        "writing dataset " + name
    );

    H5Dclose(dataset_id);
}

} // namespace eos::hdf5
