#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  // allocate spaces
  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);
  }
        
    
    // fill int sub levels
    int blockSize = 2;
    for(size_t level = 1; level < tex.mipmap.size(); ++level){
        // scale down 2^level
        int subMip_w = tex.mipmap[level].width;
        int subMip_h = tex.mipmap[level].height;

        // fill in mipmap in the level
        for(int w = 0; w < subMip_w; w++){
            for(int h = 0; h < subMip_h; h++){
                // get color values from upper level
                Color avgC(0,0,0,0);
                for(int ww = 0; ww < blockSize; ww++){
                    for(int hh = 0; hh < blockSize; hh++){
                        // get color values from upper level
                        Color sampleC((float)(tex.mipmap[level - 1].texels[4 * ((w*blockSize+ww) + (h*blockSize+hh) * subMip_w * 2)    ] / 255.f),
                                      (float)(tex.mipmap[level - 1].texels[4 * ((w*blockSize+ww) + (h*blockSize+hh) * subMip_w * 2) + 1] / 255.f),
                                      (float)(tex.mipmap[level - 1].texels[4 * ((w*blockSize+ww) + (h*blockSize+hh) * subMip_w * 2) + 2] / 255.f),
                                      (float)(tex.mipmap[level - 1].texels[4 * ((w*blockSize+ww) + (h*blockSize+hh) * subMip_w * 2) + 3] / 255.f));
                        
                        avgC += sampleC;
                    }
                }
                
                avgC = avgC * 0.25;
                
                tex.mipmap[level].texels[4 * (w + h * subMip_w)    ] = (uint8_t) (255.f * max( 0.0f, min( 1.0f, avgC.r)));
                tex.mipmap[level].texels[4 * (w + h * subMip_w) + 1] = (uint8_t) (255.f * max( 0.0f, min( 1.0f, avgC.g)));
                tex.mipmap[level].texels[4 * (w + h * subMip_w) + 2] = (uint8_t) (255.f * max( 0.0f, min( 1.0f, avgC.b)));
                tex.mipmap[level].texels[4 * (w + h * subMip_w) + 3] = (uint8_t) (255.f * max( 0.0f, min( 1.0f, avgC.a)));
            }
        }

    }
    
    
    
  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
//  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
//  for(size_t i = 1; i < tex.mipmap.size(); ++i) {
//
//    Color c = colors[i % 3];
//    MipLevel& mip = tex.mipmap[i];
//
//    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
//      float_to_uint8( &mip.texels[i], &c.r );
//    }
//  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
    
    // scale according to texture map size
    size_t w = tex.mipmap[level].width;
    size_t h = tex.mipmap[level].height;
    // u,v in [0,1]
    int i = floor(u * w - 0.5);
    int j = floor(v * h - 0.5);
    float s = u - (i + 0.5);
    float t = v - (j + 0.5);
    
    if(s > 0.5 && i + 1 < w){
        i++;
    }
    
    if(t > 0.5 && j + 1 < h){
        j++;
    }
    
    Color color((float)(tex.mipmap[level].texels[4 * (i + j * w)    ] / 255.f),
                (float)(tex.mipmap[level].texels[4 * (i + j * w) + 1] / 255.f),
                (float)(tex.mipmap[level].texels[4 * (i + j * w) + 2] / 255.f),
                (float)(tex.mipmap[level].texels[4 * (i + j * w) + 3] / 255.f));
        
    return color;

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
    
    // scale according to texture map size
    size_t w = tex.mipmap[level].width;
    size_t h = tex.mipmap[level].height;
    // u,v in [0,1]
    int i = floor(u * w - 0.5);
    int j = floor(v * h - 0.5);
    float s = u * w - (i + 0.5);
    float t = v * h - (j + 0.5);
    
    // retrieve color values in neightbourhood
    Color c00((float)(tex.mipmap[level].texels[4 * (i + j * w)    ] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + j * w) + 1] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + j * w) + 2] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + j * w) + 3] / 255.f));
    
    Color c01((float)(tex.mipmap[level].texels[4 * (i + (j + 1) * w)    ] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + (j + 1) * w) + 1] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + (j + 1) * w) + 2] / 255.f),
              (float)(tex.mipmap[level].texels[4 * (i + (j + 1) * w) + 3] / 255.f));
    
    Color c10((float)(tex.mipmap[level].texels[4 * ((i + 1) + j * w)    ] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + j * w) + 1] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + j * w) + 2] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + j * w) + 3] / 255.f));
    
    Color c11((float)(tex.mipmap[level].texels[4 * ((i + 1) + (j + 1) * w)    ] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + (j + 1) * w) + 1] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + (j + 1) * w) + 2] / 255.f),
              (float)(tex.mipmap[level].texels[4 * ((i + 1) + (j + 1) * w) + 3] / 255.f));
    
    Color c_bilinear = (1 - t) * ((1 - s) * c00 + s * c10) + t * ((1 - s) * c01 + s * c11);
    
    return c_bilinear;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

    // scale u,v according to x,y
    float level = max(9 + log2f(sqrt(max(u_scale * u_scale, v_scale * v_scale))), 0.0f);
    //cout << "level: " << level << " total level: " << tex.mipmap.size() << endl;
    // base level
    Color c0 = sample_bilinear(tex, u, v, (int)floor(level));
    // next level
    Color c1 = sample_bilinear(tex, u/2, v/2, (int)ceil(level));
    
    return (1 - level) * c0 + level * c1;
    
//  // return magenta for invalid level
//  return Color(1,0,1,1);

}

} // namespace CMU462
