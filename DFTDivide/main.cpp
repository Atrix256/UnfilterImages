#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "simple_fft/fft_settings.h"
#include "simple_fft/fft.h"

#include <vector>

static const char* c_inputFile = "../blurry.png";
static const char* c_outFile = "blurry.div.png";

const float c_kernel[] =
{
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
};
static const int c_kernelWidth = 3;
static const int c_kernelHeight = 3;

float Clamp(float value, float themin, float themax)
{
    if (value <= themin)
        return themin;
    else if (value >= themax)
        return themax;
    else
        return value;
}

void AntiConvolve(std::vector<float>& pixels, std::vector<float>& kernel, int width, int height)
{
    
}

int main(int argc, char** argv)
{
    // Load the image
    std::vector<std::vector<float>> pixels;
    int width, height, components;
    {
        unsigned char* pixels_ = stbi_load(c_inputFile, &width, &height, &components, 0);
        if (!pixels_)
        {
            printf("Could not load %s\n", c_inputFile);
            return 1;
        }

        pixels.resize(components);
        for (int c = 0; c < components; ++c)
            pixels[c].resize(width * height, 0.0f);

        for (size_t i = 0; i < width * height * components; ++i)
            pixels[i % components][i / components] = float(pixels_[i]) / 255.0f;

        stbi_image_free(pixels_);
    }

    // zero pad the kernel
    std::vector<float> kernel(width * height, 0.0f);
    {
        int startx = (width - c_kernelWidth) / 2;
        int starty = (height - c_kernelHeight) / 2;
        for (int iy = 0; iy < c_kernelHeight; ++iy)
        {
            for (int ix = 0; ix < c_kernelWidth; ++ix)
            {
                kernel[(starty + iy) * width + ix] = c_kernel[iy * c_kernelWidth + ix];
            }
        }
    }

    // Anti convolve each channel
    for (int c = 0; c < components; ++c)
        AntiConvolve(pixels[c], kernel, width, height);

    // put the image back together
    std::vector<unsigned char> outputPixels(width * height * components);
    for (size_t i = 0; i < outputPixels.size(); ++i)
        outputPixels[i] = (unsigned char)Clamp(pixels[i % components][i / components] * 256.0f, 0.0f, 255.0f);

    // write the output file
    stbi_write_png(c_outFile, width, height, components, outputPixels.data(), 0);

    return 0;
}

// TODO: should work in float, and should handle sRGB
