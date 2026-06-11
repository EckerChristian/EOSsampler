#pragma once

#include <hdf5.h>

#include <string>
#include <vector>

namespace eos::hdf5 {

hid_t create_file(const std::string& output_prefix, int rank, int file_start_index = 0);
hid_t open_file_rw(const std::string& output_prefix, int rank, int file_start_index = 0);

hid_t create_group(hid_t file_id, int index);
hid_t open_group(hid_t file_id, int index);

hid_t create_subgroup(hid_t group_id, const std::string& name);
hid_t open_subgroup(hid_t group_id, const std::string& name);

bool link_exists(hid_t group_id, const std::string& name);

void create_dataset_double(hid_t group_id, const std::string& name, hsize_t length);
void create_dataset_int(hid_t group_id, const std::string& name, hsize_t length);

void recreate_dataset_double(hid_t group_id, const std::string& name, hsize_t length);
void recreate_dataset_int(hid_t group_id, const std::string& name, hsize_t length);

std::vector<double> read_double_array(hid_t group_id, const std::string& name);
double read_double_scalar(hid_t group_id, const std::string& name);
int read_int_scalar(hid_t group_id, const std::string& name);

void write_double_array(hid_t group_id, const std::string& name, const std::vector<double>& data);
void write_double_scalar(hid_t group_id, const std::string& name, double value);
void write_int_scalar(hid_t group_id, const std::string& name, int value);

} // namespace eos::hdf5
