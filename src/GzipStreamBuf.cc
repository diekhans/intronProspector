#include "gzip_stream.hh"
#include <unistd.h>
#include <fcntl.h>
#include <stdexcept>
#include <cstring>
#include <cerrno>

void GzipStreamBuf::open_file(const std::string& filename, Mode mode) {
    filename_ = filename;
    bool use_gzip = filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz";
    mode_ = mode;

    if (use_gzip) {
        is_pipe_ = true;
        std::vector<std::string> cmd = (mode == Read)
            ? std::vector<std::string>{"gzip", "-dc", filename}
            : std::vector<std::string>{"gzip", "-c"};
        proc_ = new PipeProcess(cmd, mode == Write);
        fd_ = proc_->get_fd();
    } else {
        int flags = (mode == Read) ? O_RDONLY : (O_WRONLY | O_CREAT | O_TRUNC);
        fd_ = ::open(filename.c_str(), flags, 0666);
        if (fd_ < 0) {
            throw std::runtime_error("Failed to open file: " + filename_ + ": " + std::strerror(errno));
        }
    }

    if (mode == Read) {
        setg(in_buf_, in_buf_, in_buf_);
    } else {
        setp(out_buf_, out_buf_ + bufsize);
    }
}

GzipStreamBuf::GzipStreamBuf(const std::string& filename, Mode mode) {
    open_file(filename, mode);
}

GzipStreamBuf::~GzipStreamBuf() {
    try {
        sync();
        if (is_pipe_ && proc_) {
            proc_->close_and_wait();
            delete proc_;
        } else if (fd_ != -1) {
            if (::close(fd_) != 0) {
                // best effort
            }
        }
    } catch (...) {}
}

std::streambuf::int_type GzipStreamBuf::underflow() {
    if (mode_ != Read) {
        throw std::runtime_error("underflow called in write mode: " + filename_);
    }
    ssize_t n = ::read(fd_, in_buf_, bufsize);
    if (n < 0) {
        throw std::runtime_error("Read error on file: " + filename_ + ": " + std::strerror(errno));
    }
    if (n == 0) {
        return traits_type::eof();
    }
    setg(in_buf_, in_buf_, in_buf_ + n);
    return traits_type::to_int_type(*gptr());
}

std::streambuf::int_type GzipStreamBuf::overflow(int_type ch) {
    if (mode_ != Write) {
        throw std::runtime_error("overflow called in read mode: " + filename_);
    }
    if (ch != traits_type::eof()) {
        *pptr() = ch;
        pbump(1);
    }
    return flush_out() ? traits_type::not_eof(ch) : traits_type::eof();
}

int GzipStreamBuf::sync() {
    if (mode_ == Write) {
        flush_out();
    }
    return 0;
}

bool GzipStreamBuf::flush_out() {
    size_t n = pptr() - pbase();
    if (n > 0) {
        ssize_t written = ::write(fd_, out_buf_, n);
        if (written != static_cast<ssize_t>(n)) {
            throw std::runtime_error("Write error on file: " + filename_ + ": " + std::strerror(errno));
        }
        pbump(-static_cast<int>(n));
    }
    return true;
}
