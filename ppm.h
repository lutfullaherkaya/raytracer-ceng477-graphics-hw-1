#ifndef __ppm_h__
#define __ppm_h__
#include <stdexcept>
#include "parser.h"
void write_ppm(const char* filename, unsigned char* data, int width, int height);
class PPMWriter {
public:
    const char *filename;
    int width, height;
    FILE *outfile;
    unsigned char* data;
    size_t idx = 0;
    int lastWrittenRow = -1;
    PPMWriter(const char* filename, unsigned char* data, int width, int height) : filename(filename), width(width), height(height), data(data) {
        if ((outfile = fopen(filename, "w")) == nullptr) {
            throw std::runtime_error("Error: The ppm file cannot be opened for writing.");
        }

        (void) fprintf(outfile, "P3\n%d %d\n255\n", width, height);
    }

    void writeNextRow() {
        unsigned char color;

        for (size_t i = 0; i < width; ++i) {
            for (size_t c = 0; c < 3; ++c, ++idx) {
                color = data[idx];

                if (i == width - 1 && c == 2) {
                    (void) fprintf(outfile, "%d", color);
                } else {
                    (void) fprintf(outfile, "%d ", color);
                }
            }
        }

        (void) fprintf(outfile, "\n");
        lastWrittenRow++;

    }

    ~PPMWriter() {
        (void) fclose(outfile);
    }
};

#endif // __ppm_h__
