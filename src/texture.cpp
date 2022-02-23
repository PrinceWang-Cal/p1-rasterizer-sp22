#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {

    float level = get_level(sp);
    if (level < 0) {
              level = 0;
    } else if (level >= mipmap.size()) {
        level = float(mipmap.size() - 1);
    }
    if (sp.lsm == L_NEAREST) {
        int level_rounded = (int) round(level);
        if (sp.psm == P_NEAREST) {
            return sample_nearest(sp.p_uv, level_rounded);
        } else if (sp.psm == P_LINEAR) {
            return sample_bilinear(sp.p_uv, level_rounded);
        }
    } else if (sp.lsm == L_LINEAR) {
        int level1 = (int) floor(level);
        int level2 = (int) ceil(level);
        float level2_diff = (float) level2 - level;
        float level1_diff = (float) 1 - level2_diff;

        if (sp.psm == P_NEAREST) {
            Color level1_prime = sample_nearest(sp.p_uv, level1);
            Color level2_prime = sample_nearest(sp.p_uv, level2);
            return ((level1_prime * level1_diff) + (level2_prime * level2_diff));
        } else if (sp.psm == P_LINEAR) {
            Color level1_prime = sample_bilinear(sp.p_uv, level1);
            Color level2_prime = sample_bilinear(sp.p_uv, level2);
            return ((level1_prime * level1_diff) + (level2_prime * level2_diff));
        }
    }
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
      Vector2D distance_vec1 = (sp.p_dx_uv - sp.p_uv);
      Vector2D distance_vec2 = (sp.p_dy_uv - sp.p_uv);
      
      //scale
      distance_vec2[0] = distance_vec2[0]*width;
      distance_vec2[1] = distance_vec2[1]*height;
      distance_vec1[0] = distance_vec1[0]*width;
      distance_vec1[1] = distance_vec1[1]*height;
      
      float L = max(distance_vec1.norm(), distance_vec2.norm());
      float D = log2(L);
      
      if (D<0) {
          return 0;
      } else if (D >= mipmap.size()) {
          return (float) (mipmap.size() - 1);
      } else {
          return D;
      }
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    //if layer is invalid
      if (level<0 || level >= mipmap.size()) {
          return Color(1,0,1);
      }

    auto& mip = mipmap[level];
    int u = (int) (uv[0] * mip.width);
    int v = (int) (uv[1] * mip.height);
    return mip.get_texel((int) round(u), (int) round(v));
  }


//  Color Texture::lerp(float x, Color c0, Color c1) {
//      return (1-x) * c0 + x * c1;
//  }

Color Texture::lerp(float x, Color c0, Color c1) {
    return (1-x) * c0 + x * c1;
}

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
      if (level<0 || level >= mipmap.size()) {
          return Color(1,0,1);
      }
      auto& mip = mipmap[level];

      // first two intermediate interpolation
      int u0 = (int) floor((uv[0] * (mip.width - 1)) / (pow(2, level)));
      int v0 = (int) floor((uv[1] * (mip.height - 1)) / (pow(2, level)));

      int u1 = (int) ceil((uv[0] * (mip.width - 1)) / (pow(2, level)));
      int v1 = (int) ceil((uv[1] * (mip.height - 1)) / (pow(2, level)));

      Color u00 = mip.get_texel(u0, v0);
      Color u01 = mip.get_texel(u0, v1);
      Color u10 = mip.get_texel(u1, v0);
      Color u11 = mip.get_texel(u1, v1);

      float a = (float)((uv[0] * (width - 1)) / (pow(2, level)));
      float b = (float)((uv[1] * (height - 1)) / (pow(2, level)));

      Color u0_lerp = lerp(b-(float)v0, u00, u01);
      Color u1_lerp = lerp(b-(float)v0, u10, u11);

    return lerp(a-(float)u0, u0_lerp, u1_lerp);
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
