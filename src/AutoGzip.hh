#ifndef autogzip_hh
#define autogzip_hh

#include <fstream>
#include <memory>
#include <string>
#include "zfstream.hh"  // defines gzifstream, gzofstream

static bool isGzipName(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

class AutoGzipInput : public std::istream {
public:
    explicit AutoGzipInput(const std::string& filename)
        : std::istream(nullptr), filename_(filename), is_gz_(isGzipName(filename)) {
        if (is_gz_) {
            gz_ = std::make_unique<gzifstream>(filename_);
            if (!gz_->is_open())
                throw std::ios_base::failure("Failed to open gzip file: " + filename_);
            rdbuf(gz_->rdbuf());
        } else {
            file_ = std::make_unique<std::ifstream>(filename_);
            if (!file_->is_open())
                throw std::ios_base::failure("Failed to open file: " + filename_);
            rdbuf(file_->rdbuf());
        }
    }

    bool is_compressed() const { return is_gz_; }
    const std::string& filename() const { return filename_; }

private:
    const std::string filename_;
    bool is_gz_;
    std::unique_ptr<gzifstream> gz_;
    std::unique_ptr<std::ifstream> file_;
};

class AutoGzipOutput : public std::ostream {
public:
    explicit AutoGzipOutput(const std::string& filename)
        : std::ostream(nullptr), filename_(filename), is_gz_(isGzipName(filename)) {
        if (is_gz_) {
            gz_ = std::make_unique<gzofstream>(filename_);
            if (!gz_->is_open())
                throw std::ios_base::failure("Failed to open gzip file: " + filename_);
            rdbuf(gz_->rdbuf());
        } else {
            file_ = std::make_unique<std::ofstream>(filename_);
            if (!file_->is_open())
                throw std::ios_base::failure("Failed to open file: " + filename_);
            rdbuf(file_->rdbuf());
        }
    }

    bool is_compressed() const { return is_gz_; }
    const std::string& filename() const { return filename_; }

private:
    const std::string filename_;
    bool is_gz_;
    std::unique_ptr<gzofstream> gz_;
    std::unique_ptr<std::ofstream> file_;
};

#endif
