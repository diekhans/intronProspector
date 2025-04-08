#ifndef GZIP_STREAM_HH
#define GZIP_STREAM_HH

#include "pipe_process.hh"
#include <streambuf>
#include <istream>
#include <ostream>
#include <string>

class GzipStreamBuf : public std::streambuf {
public:
    enum Mode { Read, Write };
    GzipStreamBuf(const std::string& filename, Mode mode);
    ~GzipStreamBuf() override;

protected:
    int_type underflow() override;
    int_type overflow(int_type ch) override;
    int sync() override;

private:
    static const size_t bufsize = 4096;
    char in_buf_[bufsize], out_buf_[bufsize];
    Mode mode_;
    int fd_ = -1;
    bool is_pipe_ = false;
    PipeProcess* proc_ = nullptr;
    std::string filename_;

    bool flush_out();
    void open_file(const std::string& filename, Mode mode);
};

class GzipInputStream : public std::istream {
    GzipStreamBuf buf_;
public:
    explicit GzipInputStream(const std::string& filename)
        : std::istream(nullptr), buf_(filename, GzipStreamBuf::Read) {
        rdbuf(&buf_);
    }
};

class GzipOutputStream : public std::ostream {
    GzipStreamBuf buf_;
public:
    explicit GzipOutputStream(const std::string& filename)
        : std::ostream(nullptr), buf_(filename, GzipStreamBuf::Write) {
        rdbuf(&buf_);
    }
};

#endif
