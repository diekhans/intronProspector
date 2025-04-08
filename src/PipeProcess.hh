#ifndef PIPE_PROCESS_HH
#define PIPE_PROCESS_HH

#include <string>
#include <vector>

class PipeProcess {
public:
    PipeProcess(const std::vector<std::string>& argv, bool write_mode);
    int get_fd() const;
    void close_and_wait();
    ~PipeProcess();

private:
    int fd_;
    pid_t pid_;
    bool write_mode_;
    std::string command_;
};

#endif
