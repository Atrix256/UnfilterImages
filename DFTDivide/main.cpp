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

float c_kernel[] =
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

inline size_t NextPowerOf2(size_t n)
{
    size_t ret = 1;
    while (ret < n)
        ret = ret << 1;
    return ret;
}

struct ComplexImage2D
{
    ComplexImage2D(size_t w = 0, size_t h = 0)
    {
        Resize(w, h);
    }

    size_t m_width;
    size_t m_height;
    std::vector<complex_type> pixels;

    void Resize(size_t w, size_t h)
    {
        m_width = w;
        m_height = h;
        pixels.resize(w * h, real_type(0.0f));
    }

    complex_type& operator()(size_t x, size_t y)
    {
        return pixels[y * m_width + x];
    }

    const complex_type& operator()(size_t x, size_t y) const
    {
        return pixels[y * m_width + x];
    }

    void IndexToCoordinates(size_t index, size_t& x, size_t& y) const
    {
        y = index / m_width;
        x = index % m_width;
    }
};

ComplexImage2D ZeroPad(const ComplexImage2D& image, size_t width, size_t height)
{
    size_t offsetx = (width - image.m_width) / 2;
    size_t offsety = (height - image.m_height) / 2;

    ComplexImage2D ret(width, height);

    for (size_t iy = 0; iy < image.m_height; ++iy)
    {
        const complex_type* src = &image.pixels[iy * image.m_width];
        complex_type* dest = &ret.pixels[(iy + offsety) * ret.m_width + offsetx];
        memcpy(dest, src, sizeof(complex_type) * image.m_width);
    }

    return ret;
}

ComplexImage2D UnZeroPad(const ComplexImage2D& image, size_t width, size_t height)
{
    size_t offsetx = (image.m_width - width) / 2;
    size_t offsety = (image.m_height - height) / 2;

    ComplexImage2D ret(width, height);

    for (size_t iy = 0; iy < height; ++iy)
    {
        const complex_type* src = &image.pixels[(iy + offsety) * image.m_width + offsetx];
        complex_type* dest = &ret.pixels[iy * ret.m_width];
        memcpy(dest, src, sizeof(complex_type) * width);
    }

    return ret;
}

ComplexImage2D ShiftImage(const ComplexImage2D& image, size_t offsetx, size_t offsety)
{
    ComplexImage2D ret = image;

    for (size_t iy = 0; iy < image.m_height; ++iy)
    {
        size_t srciy = (iy + offsety) % image.m_height;
        for (size_t ix = 0; ix < image.m_width; ++ix)
        {
            size_t srcix = (ix + offsetx) % image.m_width;
            ret(ix, iy) = image(srcix, srciy);
        }
    }

    return ret;
}

void AntiConvolve(ComplexImage2D& pixels, ComplexImage2D& kernel)
{
    // DFT the image
    ComplexImage2D pixelsFT;
    const char* error = nullptr;
    pixelsFT.Resize(pixels.m_width, pixels.m_height);
    simple_fft::FFT(pixels, pixelsFT, pixels.m_width, pixels.m_height, error);
    if (error)
        printf(error);

    // DFT the kernel
    ComplexImage2D kernelFT;
    error = nullptr;
    kernelFT.Resize(kernel.m_width, kernel.m_height);
    simple_fft::FFT(kernel, kernelFT, kernel.m_width, kernel.m_height, error);
    if (error)
        printf(error);

    // Anti convolve
    for (size_t i = 0; i < pixelsFT.pixels.size(); ++i)
        pixelsFT.pixels[i] *= kernelFT.pixels[i];

    // Inverse DFT the image
    error = nullptr;
    simple_fft::IFFT(pixelsFT, pixels, pixelsFT.m_width, pixelsFT.m_height, error);
    if (error)
        printf(error);
}

int main(int argc, char** argv)
{
    // Load the image
    std::vector<ComplexImage2D> pixels;
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
            pixels[c].Resize(width, height);

        for (size_t i = 0; i < width * height * components; ++i)
            pixels[i % components].pixels[i / components] = float(pixels_[i]) / 255.0f;

        stbi_image_free(pixels_);
    }

    // make a normalized kernel
    ComplexImage2D kernel(c_kernelWidth, c_kernelHeight);
    {
        float sum = 0.0f;
        for (int iy = 0; iy < c_kernelHeight; ++iy)
            for (int ix = 0; ix < c_kernelWidth; ++ix)
                sum += c_kernel[iy * c_kernelWidth + ix];

        for (int iy = 0; iy < c_kernelHeight; ++iy)
            for (int ix = 0; ix < c_kernelWidth; ++ix)
                kernel.pixels[iy * c_kernelWidth + ix] /= sum;
    }

    // calculate the size that the images need to be, to be multiplied in frequency space
    size_t desiredWidth = NextPowerOf2(pixels[0].m_width + kernel.m_width + 1);
    size_t desiredHeight = NextPowerOf2(pixels[0].m_height + kernel.m_height + 1);

    // zero pad the images to be the right size
    for (int c = 0; c < components; ++c)
        pixels[c] = ZeroPad(pixels[c], desiredWidth, desiredHeight);
    kernel = ZeroPad(kernel, desiredWidth, desiredHeight);

    // shift the kernel because it needs to be centered at (0,0)
    kernel = ShiftImage(kernel, kernel.m_width / 2, kernel.m_height / 2);

    // Anti convolve each channel and remove the extra bordering we put on
    for (int c = 0; c < components; ++c)
    {
        AntiConvolve(pixels[c], kernel);
        pixels[c] = UnZeroPad(pixels[c], width, height);
    }

    // put the image back together
    std::vector<unsigned char> outputPixels(width * height * components);
    for (size_t i = 0; i < outputPixels.size(); ++i)
        outputPixels[i] = (unsigned char)Clamp((float)pixels[i % components].pixels[i / components].real() * 256.0f, 0.0f, 255.0f);

    // write the output file
    stbi_write_png(c_outFile, width, height, components, outputPixels.data(), 0);

    return 0;
}

// TODO: should work in float, and should handle sRGB
// TODO: do i need to center the kernel at (0,0)? check that blog post
