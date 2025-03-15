/*  tsv.hh -- really simplistic TSV reader

    Copyright (c) 2019-2025, Mark Diekhans, University of California, Santa Cruz
**/
#ifndef tsv_hh
#define tsv_hh
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <stdexcept>
#include "type_ops.hh"

using namespace std;


/* read TSV one row at a time, indexed by column name */
class Tsv {
    private:
    typedef map<string, int> ColMap;
    
    const string tsv_file_;
    ifstream tsv_fh_;
    int num_columns_;
    ColMap col_map_;  // name to column index
    vector<string> row_;   // current row

    vector<string> split_line(const string& line) {
        vector<string> words;
        size_t iprevious = 0;
        size_t icurrent = line.find('\t');
        while (icurrent != string::npos) {
            words.push_back(line.substr(iprevious, icurrent - iprevious));
            iprevious = icurrent + 1;
            icurrent = line.find('\t', iprevious);
        }
        words.push_back(line.substr(iprevious, icurrent - iprevious));
        return words;
    }
    
    void read_header() {
        string line;
        if (not getline(tsv_fh_, line)) {
            throw runtime_error("empty TSV file: " + tsv_file_);
        }
        vector<string> words = split_line(line);
        num_columns_ = words.size();
        for (int i = 0; i < num_columns_; i++) {
            col_map_[words[i]] = i;
        }
    }

    public:
    // open file and read headers
    Tsv(const string& tsv_file):
        tsv_file_(tsv_file),
        tsv_fh_(tsv_file.c_str()),
        num_columns_(0) {
        if (tsv_fh_.bad() | tsv_fh_.fail()) {
            throw runtime_error("can't open " + tsv_file);
        }
        read_header();
    }

    // read the next row
    bool next_row() {
        string line;
        if (not getline(tsv_fh_, line)) {
            return false;
        }
        row_ = split_line(line);
        if (row_.size() != num_columns_) {
            throw runtime_error("expected " + int_to_string(num_columns_) +  " columns, got " + int_to_string(num_columns_) + ":" + tsv_file_);
        }
        return true;
    }

    // does column exist?
    bool have_col(const string& col_name) const {
        ColMap::const_iterator it = col_map_.find(col_name);
        return (it != col_map_.end());
    }
    
    // get the specified column
    const string& get_col(const string& col_name) const {
        ColMap::const_iterator it = col_map_.find(col_name);
        if (it == col_map_.end()) {
            throw runtime_error("column '" + col_name + "' not found: " + tsv_file_);
        }
        return row_[it->second];
    }
    
    // get the specified column as an integer
    uint32_t get_col_int(const string& col_name) const {
        return string_to_uint(get_col(col_name));
    }

    // get the specified column as an 640bit
    uint64_t get_col_int64(const string& col_name) const {
        return string_to_uint64(get_col(col_name));
    }

    // get the specified column as a float
    float get_col_float(const string& col_name) const {
        return string_to_float(get_col(col_name));
    }

};

#endif
