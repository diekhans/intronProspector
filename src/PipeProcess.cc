#include "pipe_process.hh"
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include <sstream>

static std::string join_command(const std::vector<std::string>& argv) {
    std::ostringstream oss;
    for (const auto& arg : argv) {
        oss << arg << " ";
    }
    return oss.str();
}

PipeProcess::PipeProcess(const std::vector<std::string>& argv, bool write_mode)
    : fd_(-1), pid_(-1), write_mode_(write_mode), command_(join_command(argv)) {
    int pipefd[2];
    if (pipe(pipefd) != 0) {
        throw std::runtime_error("pipe() failed for command [" + command_ + "]: " + std::strerror(errno));
    }

    pid_ = fork();
    if (pid_ < 0) {
        close(pipefd[0]);
        close(pipefd[1]);
        throw std::runtime_error("fork() failed for command [" + command_ + "]: " + std::strerror(errno));
    }

    if (pid_ == 0) {
        if (write_mode_) {
            close(pipefd[1]);
            dup2(pipefd[0], STDIN_FILENO);
        } else {
            close(pipefd[0]);
            dup2(pipefd[1], STDOUT_FILENO);
        }

        for (int fd = 3; fd < 1024; ++fd) close(fd);

        std::vector<char*> args;
        for (const auto& s : argv) {
            args.push_back(const_cast<char*>(s.c_str()));
        }
        args.push_back(nullptr);
        execvp(args[0], args.data());
        _exit(127); // exec failed
    }

    if (write_mode_) {
        close(pipefd[0]);
        fd_ = pipefd[1];
    } else {
        close(pipefd[1]);
        fd_ = pipefd[0];
    }
}

int PipeProcess::get_fd() const {
    return fd_;
}

void PipeProcess::close_and_wait() {
    if (fd_ != -1) {
        ::close(fd_);
        fd_ = -1;
    }
    if (pid_ != -1) {
        int status;
        if (waitpid(pid_, &status, 0) < 0) {
            throw std::runtime_error("waitpid() failed for command [" + command_ + "]: " + std::strerror(errno));
        }
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            throw std::runtime_error("child process failed: [" + command_ + "]");
        }
        pid_ = -1;
    }
}

PipeProcess::~PipeProcess() {
    try {
        close_and_wait();
    } catch (...) {
        // suppress destructor exceptions
    }
}
