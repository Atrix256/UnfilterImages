#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include <direct.h>
#include <vector>
#include <chrono>

#define SAVE_HDR() true
#define OUTPUT_FILTERED_IMAGE() false

static const char* c_inputFileDir = "../";
static const char* c_inputFileBaseName = "clear64";

float c_kernel[] =
{
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f
};

static const int c_kernelWidth = 2;
static const int c_kernelHeight = 1;

static const int c_kernelRadiusX = int(c_kernelWidth / 2);
static const int c_kernelRadiusY = int(c_kernelHeight / 2);

float Clamp(float value, float themin, float themax)
{
    if (value <= themin)
        return themin;
    else if (value >= themax)
        return themax;
    else
        return value;
}

struct AugmentedMatrix
{
    AugmentedMatrix(uint64_t numUnknowns)
        : m_numUnknowns(numUnknowns)
    {
        m_rowsPerPage = c_pageSizeBytes / ((numUnknowns + 1) * sizeof(float));
        uint64_t pageSizeFloats = m_rowsPerPage * (numUnknowns + 1);
        uint64_t numPages = (numUnknowns + m_rowsPerPage - uint64_t(1)) / m_rowsPerPage;
        m_pages.resize(numPages);
        {
            for (uint64_t i = 0; i < numPages; ++i)
            {
                if (i < numPages - 1)
                {
                    m_pages[i] = new float[pageSizeFloats];
                    memset(m_pages[i], 0, pageSizeFloats * sizeof(float));
                }
                else
                {
                    uint64_t remainderRows = numUnknowns - (numPages - 1) * m_rowsPerPage;
                    uint64_t remainderRowsFloats = remainderRows * (numUnknowns + 1);
                    m_pages[i] = new float[remainderRowsFloats];
                    memset(m_pages[i], 0, remainderRowsFloats * sizeof(float));
                }
            }
        }
        m_temporaryRow.resize(numUnknowns + 1);
    }

    float* GetRow(uint64_t rowIndex)
    {
        uint64_t pageIndex = rowIndex / m_rowsPerPage;
        uint64_t pageOffsetFloats = (rowIndex % m_rowsPerPage) * (m_numUnknowns + 1);
        return &m_pages[pageIndex][pageOffsetFloats];
    }

    ~AugmentedMatrix()
    {
        for (float* p : m_pages)
            delete[] p;
        m_pages.clear();
    }

    std::vector<float> Solve()
    {
        print("out/start");

        // solve each column
        {
            int lastPercent = -1;
            for (uint64_t index = 0; index < m_numUnknowns; ++index)
            {
                int percent = int(10000.0f * float(index) / float(m_numUnknowns - 1));
                if (lastPercent != percent)
                {
                    printf("\r%0.2f%%", float(percent) / 100.0f);
                    lastPercent = percent;
                }
                SolveColumn(index);
            }
            printf("\r100.00%%\n");
        }

        print("out/end");

        // now gather the results
        std::vector<float> results(m_numUnknowns);
        for (uint64_t index = 0; index < m_numUnknowns; ++index)
            results[index] = GetRow(index)[m_numUnknowns];
        return results;
    }

private:

    void print(const char* fileNameBase)
    {
        static int count = -1;
        count++;
        char fileName[1024];
        sprintf_s(fileName, "%s.%i.csv", fileNameBase, count);

        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+b");
        if (!file)
        {
            printf("could not open %s", fileName);
            ((int*)0)[0] = 0;
        }

        for (uint64_t i = 0; i < m_numUnknowns; ++i)
        {
            float* row = GetRow(i);
            for (uint64_t j = 0; j < m_numUnknowns + 1; ++j)
                fprintf(file, "\"%f\",", row[j]);
            fprintf(file, "\n");
        }
        fclose(file);
    }

    void SolveColumn(uint64_t columnIndex)
    {
        //1) Find the first row that has non zero in this column
        uint64_t rowIndex = columnIndex;
        float* row = GetRow(rowIndex);
        while (rowIndex < m_numUnknowns && row[columnIndex] == 0.0f)
        {
            rowIndex++;
            row = GetRow(rowIndex);
        }
        if (rowIndex == m_numUnknowns)
        {
            // TODO: what to do about this??
            printf("ERROR!");
            return;
            //((int*)0)[0] = 0;
        }

        //2) Swap this row with the row it should be in
        SwapRows(rowIndex, columnIndex);
        rowIndex = columnIndex;
        row = GetRow(rowIndex);

        //3) Make this row have a 1 in the column
        {
            float rowValue = row[columnIndex];
            for (uint64_t index = 0; index < m_numUnknowns + 1; ++index)
                row[index] /= rowValue;
            row[columnIndex] = 1.0f; // avoid numerical issues
        }

        //3) Make other rows have a zero in this column by subtracting a multiple of this row from them
        for (uint64_t index = 0; index < m_numUnknowns; ++index)
        {
            if (index == rowIndex)
                continue;

            float* zeroRow = GetRow(index);
            if (zeroRow[columnIndex] == 0.0f)
                continue;

            RowSubtract(zeroRow, row, zeroRow[columnIndex]);
            zeroRow[columnIndex] = 0.0f; // avoid numerical issues
        }
    }

    void RowSubtract(float* row, float* srcRow, float multiplier)
    {
        for (uint64_t index = 0; index < m_numUnknowns + 1; ++index)
            row[index] -= srcRow[index] * multiplier;
    }

    void SwapRows(uint64_t rowIndexA, uint64_t rowIndexB)
    {
        // Nothing to do if the row is in the right place
        if (rowIndexA == rowIndexB)
            return;

        float* rowA = GetRow(rowIndexA);
        float* rowB = GetRow(rowIndexB);

        memcpy(m_temporaryRow.data(), rowA, (m_numUnknowns + 1) * sizeof(float));
        memcpy(rowA, rowB, (m_numUnknowns + 1) * sizeof(float));
        memcpy(rowB, m_temporaryRow.data(), (m_numUnknowns + 1) * sizeof(float));
    }

    static const uint64_t c_pageSizeBytes = uint64_t(2) * uint64_t(1024) * uint64_t(1024) * uint64_t(1024);

    std::vector<float*> m_pages;
    std::vector<float> m_temporaryRow;
    uint64_t m_numUnknowns;
    uint64_t m_rowsPerPage;
};

std::vector<float> Convolve(const float* image, int imageWidth, int imageHeight, const float* kernel, int kernelWidth, int kernelHeight)
{
    int startx = -(int)kernelWidth / 2;
    int starty = -(int)kernelHeight / 2;

    std::vector<float> ret(imageWidth * imageHeight, 0.0f);

    for (int iy = 0; iy < (int)imageHeight; ++iy)
    {
        for (int ix = 0; ix < (int)imageWidth; ++ix)
        {
            for (int ky = 0; ky < kernelHeight; ++ky)
            {
                int sy = (ky + iy + starty + (int)imageHeight) % (int)imageHeight;

                for (int kx = 0; kx < kernelWidth; ++kx)
                {
                    int sx = (kx + ix + startx + (int)imageWidth) % (int)imageWidth;

                    ret[iy * imageWidth + ix] += image[sy * imageWidth + sx] * kernel[ky * kernelWidth + kx];
                }
            }
        }
    }

    return ret;
}

int main(int argc, char** argv)
{
    _mkdir("out");
    char fileName[1024];

    // normalize the kernel
    {
        float sum = 0.0f;
        for (int i = 0; i < c_kernelWidth * c_kernelHeight; ++i)
            sum += c_kernel[i];
        for (int i = 0; i < c_kernelWidth * c_kernelHeight; ++i)
            c_kernel[i] /= sum;
    }

    // load the image into memory
    std::vector<std::vector<float>> pixels;
    int width, height, components;
    {
        sprintf(fileName, "%s%s.png", c_inputFileDir, c_inputFileBaseName);
        unsigned char* pixels_ = stbi_load(fileName, &width, &height, &components, 0);
        if (!pixels_)
        {
            printf("Could not load %s\n", fileName);
            return 1;
        }

        pixels.resize(components);
        for (int i = 0; i < components; ++i)
            pixels[i].resize(width * height);

        pixels.resize(components);
        for (int i = 0; i < width * height * components; ++i)
            pixels[i % components][i / components] = float(pixels_[i]) / 255.0f;

        stbi_image_free(pixels_);

        printf("Loaded %s\n%ix%i with %i components\n\n", fileName, width, height, components);
    }

#if OUTPUT_FILTERED_IMAGE()
    {
        // Convolve
        std::vector<std::vector<float>> pixelsOut(components);
        for (int ic = 0; ic < components; ++ic)
            pixelsOut[ic] = Convolve(pixels[ic].data(), width, height, c_kernel, c_kernelWidth, c_kernelHeight);

        // put the image back together
        std::vector<unsigned char> output(width * height * components);
        for (size_t i = 0; i < output.size(); ++i)
            output[i] = (unsigned char)Clamp(pixelsOut[i % components][i / components] * 256.0f, 0.0f, 255.0f);

        // save it out
        sprintf(fileName, "out/%s.filtered.png", c_inputFileBaseName);
        stbi_write_png(fileName, width, height, components, output.data(), 0);
    }
#endif

    // handle each component individually
    for (int ic = 0; ic < components; ++ic)
    {
        printf("Solving component %i / %i\n", ic + 1, components);
        std::vector<float>& src = pixels[ic];

        // make the augmented matrix
        AugmentedMatrix augmentedMatrix(size_t(width) * size_t(height));
        for (int i = 0; i < width * height; i++)
        {
            float* row = augmentedMatrix.GetRow(i);

            // set the value of the new pixel
            row[width * height] = src[i];

            // put the kernel values in
            int pixelX = i % width;
            int pixelY = i / width;
            for (int ky = 0; ky < c_kernelHeight; ++ky)
            {
                int kernelOffsetY = ky - c_kernelRadiusY;
                int destY = (pixelY + kernelOffsetY + height) % height;

                for (int kx = 0; kx < c_kernelWidth; ++kx)
                {
                    int kernelOffsetX = kx - c_kernelRadiusX;
                    int destX = (pixelX + kernelOffsetX + width) % width;

                    row[destY * width + destX] = c_kernel[ky * c_kernelWidth + kx];
                }
            }
        }

        // Solve it
        pixels[ic] = augmentedMatrix.Solve();
    }

    // put the image back together
    std::vector<unsigned char> output(width * height * components);
    std::vector<float> outputF(width * height * components);
    for (size_t i = 0; i < output.size(); ++i)
    {
        outputF[i] = pixels[i % components][i / components];
        output[i] = (unsigned char)Clamp(pixels[i % components][i / components] * 256.0f, 0.0f, 255.0f);
    }

    // save the image to disk
    sprintf(fileName, "out/%s.unfilter.png", c_inputFileBaseName);
    stbi_write_png(fileName, width, height, components, output.data(), 0);

#if SAVE_HDR()
    sprintf(fileName, "out/%s.unfilter.hdr", c_inputFileBaseName);
    stbi_write_hdr(fileName, width, height, components, outputF.data());
#endif

    // refilter
    {
        // Convolve
        std::vector<std::vector<float>> pixelsOut(components);
        for (int ic = 0; ic < components; ++ic)
            pixelsOut[ic] = Convolve(pixels[ic].data(), width, height, c_kernel, c_kernelWidth, c_kernelHeight);

        // put the image back together
        std::vector<unsigned char> output(width * height * components);
        std::vector<float> outputF(width * height * components);
        for (size_t i = 0; i < output.size(); ++i)
        {
            output[i] = (unsigned char)Clamp(pixelsOut[i % components][i / components] * 256.0f, 0.0f, 255.0f);
            outputF[i] = pixelsOut[i % components][i / components];
        }

        // save it out
        sprintf(fileName, "out/%s.refilter.png", c_inputFileBaseName);
        stbi_write_png(fileName, width, height, components, output.data(), 0);

#if SAVE_HDR()
        sprintf(fileName, "out/%s.refilter.hdr", c_inputFileBaseName);
        stbi_write_hdr(fileName, width, height, components, outputF.data());
#endif
    }

    return 0;
}

// NOTE: 512x512 was too much to fit into memory so went to 256x256 and then 64x64. numPixels is squared, so it's dimension^4.
// could try using https://petsc.org/release/overview/
// or something else from https://twitter.com/Atrix256/status/1518338372829204480?s=20&t=MQXKkpGKHmqQMr38VkmBiw
