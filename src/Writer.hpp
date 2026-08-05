#ifndef WRITER_HPP
#define WRITER_HPP

#include <cstdio>
#include <ios>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <zlib.h>

namespace sanity
{

enum class OutputFormat
{
    Auto,
    PlainText,
    Gzip
};

class WriterBackend
{
public:
    virtual ~WriterBackend() {}

    virtual void write(const char *data, std::size_t size) = 0;
    virtual void flush() = 0;
    virtual void close() = 0;
    virtual bool good() const = 0;
};

class PlainTextWriterBackend : public WriterBackend
{
public:
    PlainTextWriterBackend() : file_(nullptr) {}

    explicit PlainTextWriterBackend(const std::string &filename) : file_(nullptr)
    {
        open(filename);
    }

    ~PlainTextWriterBackend()
    {
        close();
    }

    void open(const std::string &filename)
    {
        close();
        file_ = std::fopen(filename.c_str(), "wb");
        if (file_ == nullptr)
        {
            throw std::runtime_error("Failed to open output file: " + filename);
        }
    }

    virtual void write(const char *data, std::size_t size)
    {
        if (file_ == nullptr)
        {
            throw std::runtime_error("Attempted to write to a closed plain-text file");
        }
        if (size == 0)
        {
            return;
        }

        const std::size_t written = std::fwrite(data, 1, size, file_);
        if (written != size)
        {
            throw std::runtime_error("Failed while writing plain-text output");
        }
    }

    virtual void flush()
    {
        if (file_ != nullptr && std::fflush(file_) != 0)
        {
            throw std::runtime_error("Failed to flush plain-text output");
        }
    }

    virtual void close()
    {
        if (file_ != nullptr)
        {
            std::fclose(file_);
            file_ = nullptr;
        }
    }

    virtual bool good() const
    {
        return file_ != nullptr && std::ferror(file_) == 0;
    }

private:
    std::FILE *file_;
};

class GzipWriterBackend : public WriterBackend
{
public:
    GzipWriterBackend() : file_(nullptr) {}

    explicit GzipWriterBackend(const std::string &filename) : file_(nullptr)
    {
        open(filename);
    }

    ~GzipWriterBackend()
    {
        close();
    }

    void open(const std::string &filename)
    {
        close();
        file_ = gzopen(filename.c_str(), "wb");
        if (file_ == nullptr)
        {
            throw std::runtime_error("Failed to open gzip output file: " + filename);
        }
    }

    virtual void write(const char *data, std::size_t size)
    {
        if (file_ == nullptr)
        {
            throw std::runtime_error("Attempted to write to a closed gzip file");
        }
        if (size == 0)
        {
            return;
        }

        std::size_t offset = 0;
        while (offset < size)
        {
            const std::size_t chunk_size = size - offset > static_cast<std::size_t>(std::numeric_limits<int>::max())
                ? static_cast<std::size_t>(std::numeric_limits<int>::max())
                : size - offset;
            const int written = gzwrite(file_, data + offset, static_cast<unsigned int>(chunk_size));
            if (written <= 0)
            {
                throw std::runtime_error("Failed while writing gzip output");
            }
            offset += static_cast<std::size_t>(written);
        }
    }

    virtual void flush()
    {
        if (file_ != nullptr && gzflush(file_, Z_SYNC_FLUSH) != Z_OK)
        {
            throw std::runtime_error("Failed to flush gzip output");
        }
    }

    virtual void close()
    {
        if (file_ != nullptr)
        {
            gzclose(file_);
            file_ = nullptr;
        }
    }

    virtual bool good() const
    {
        return file_ != nullptr;
    }

private:
    gzFile file_;
};

class Writer
{
public:
    Writer() : format_(OutputFormat::PlainText) {}

    explicit Writer(const std::string &filename, OutputFormat format = OutputFormat::Auto)
        : format_(OutputFormat::PlainText)
    {
        open(filename, format);
    }

    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;

    Writer(Writer &&other) noexcept
        : backend_(std::move(other.backend_)),
          format_state_(other.format_state_.str()),
          filename_(std::move(other.filename_)),
          format_(other.format_)
    {
        format_state_.copyfmt(other.format_state_);
    }

    Writer &operator=(Writer &&other) noexcept
    {
        if (this != &other)
        {
            close();
            backend_ = std::move(other.backend_);
            format_state_.str(other.format_state_.str());
            format_state_.copyfmt(other.format_state_);
            filename_ = std::move(other.filename_);
            format_ = other.format_;
        }
        return *this;
    }

    ~Writer()
    {
        close();
    }

    void open(const std::string &filename, OutputFormat format = OutputFormat::Auto)
    {
        close();

        filename_ = filename;
        format_ = (format == OutputFormat::Auto) ? detect_format(filename) : format;

        switch (format_)
        {
        case OutputFormat::PlainText:
            backend_.reset(new PlainTextWriterBackend(filename));
            break;
        case OutputFormat::Gzip:
            backend_.reset(new GzipWriterBackend(filename));
            break;
        default:
            throw std::runtime_error("Unsupported output format");
        }
    }

    void close()
    {
        if (backend_)
        {
            backend_->close();
            backend_.reset();
        }
        filename_.clear();
    }

    void flush()
    {
        require_open();
        backend_->flush();
    }

    bool is_open() const
    {
        return backend_.get() != nullptr;
    }

    bool good() const
    {
        return backend_ && backend_->good();
    }

    const std::string &filename() const
    {
        return filename_;
    }

    OutputFormat format() const
    {
        return format_;
    }

    void write_raw(const char *data, std::size_t size)
    {
        require_open();
        backend_->write(data, size);
    }

    void write_raw(const std::string &data)
    {
        write_raw(data.data(), data.size());
    }

    template <typename T>
    Writer &operator<<(const T &value)
    {
        require_open();

        std::ostringstream chunk;
        chunk.copyfmt(format_state_);
        chunk << value;

        const std::string rendered = chunk.str();
        if (!rendered.empty())
        {
            backend_->write(rendered.data(), rendered.size());
        }

        format_state_.width(0);
        return *this;
    }

    Writer &operator<<(std::ostream &(*manip)(std::ostream &))
    {
        require_open();

        std::ostringstream chunk;
        chunk.copyfmt(format_state_);
        manip(chunk);

        const std::string rendered = chunk.str();
        if (!rendered.empty())
        {
            backend_->write(rendered.data(), rendered.size());
        }
        else
        {
            backend_->flush();
        }

        format_state_.width(0);
        return *this;
    }

    Writer &operator<<(std::ios &(*manip)(std::ios &))
    {
        manip(format_state_);
        return *this;
    }

    Writer &operator<<(std::ios_base &(*manip)(std::ios_base &))
    {
        manip(format_state_);
        return *this;
    }

    static OutputFormat detect_format(const std::string &filename)
    {
        if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz")
        {
            return OutputFormat::Gzip;
        }
        return OutputFormat::PlainText;
    }

private:
    void require_open() const
    {
        if (!backend_)
        {
            throw std::runtime_error("Writer is not open");
        }
    }

    std::unique_ptr<WriterBackend> backend_;
    std::ostringstream format_state_;
    std::string filename_;
    OutputFormat format_;
};

} // namespace sanity

#endif