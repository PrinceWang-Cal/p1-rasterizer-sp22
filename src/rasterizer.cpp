#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


    sample_buffer[y * width*sample_rate + x*sample_rate] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

    bool RasterizerImp::point_in_triangle(float x, float y,
      float x0, float y0,
      float x1, float y1,
      float x2, float y2,
      float dx0, float dy0,
      float dx1, float dy1,
      float dx2, float dy2) {
        float l0 = -(x-x0) * dy0 + (y-y0) * dx0;
        float l1 = -(x-x1) * dy1 + (y-y1) * dx1;
        float l2 = -(x-x2) * dy2 + (y-y2) * dx2;
        
        return (((l0 >= 0) && (l1 >= 0)) && (l2 >= 0)) ||
        (((l0 <= 0) && (l1 <= 0)) && (l2 <= 0));
    }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      
      unsigned int sample_rate_uns = sample_rate;
      float increment = float(1)/float(sqrt(sample_rate_uns));
      float increment_center = increment/float(2);

    //locate the bounding boxes of the triangle
      float minx = floor(std::min({x0,x1,x2}));
      float miny = floor(std::min({y0,y1,y2}));

      float maxx = ceil(std::max({x0,x1,x2}));
      float maxy = ceil(std::max({y0,y1,y2}));

      float dx0 = x1-x0;
      float dx1 = x2-x1;
      float dx2 = x0-x2;

      float dy0 = y1-y0;
      float dy1 = y2-y1;
      float dy2 = y0-y2;
      

      for (float x = minx; x < maxx; x++) {
          for (float y = miny; y < maxy; y++) {
              int sx = (int)floor(x);
              int sy = (int)floor(y);
              
              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;
              

              for (int subx = 0; subx < sqrt(sample_rate); subx++) {
                  for (int suby = 0; suby < sqrt(sample_rate); suby++) {
                      int idx_sub = suby * sqrt(sample_rate) + subx;
                      int idx = y*width*sample_rate + x*sample_rate + idx_sub;

                      //check bound

                      if (idx < 0 || idx >= width*height*sample_rate) {
                          continue; //we skip it
                      } else {
                          
                          float add_to_x = float(subx) * increment + increment_center;
                          float add_to_y = float(suby) * increment + increment_center;
                          
                          float y_prime = float(y) + add_to_y;
                          float x_prime = float(x) + add_to_x;

                          //check bound
                          if (x_prime<0 || y_prime < 0 || x_prime >= width || y_prime >= height) {
                              continue;
                          }

                          if (point_in_triangle(x_prime, y_prime, x0, y0, x1, y1, x2, y2, dx0, dy0, dx1, dy1, dx2, dy2)) {
                              sample_buffer[idx] = color;
                          }

                      }
                  }
              }

          }
      }
      

  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      unsigned int sample_rate_uns = sample_rate;
      float increment = float(1)/float(sqrt(sample_rate_uns));
      float increment_center = increment/float(2);

    //locate the bounding boxes of the triangle
      float minx = floor(std::min({x0,x1,x2}));
      float miny = floor(std::min({y0,y1,y2}));

      float maxx = ceil(std::max({x0,x1,x2}));
      float maxy = ceil(std::max({y0,y1,y2}));

      float dx0 = x1-x0;
      float dx1 = x2-x1;
      float dx2 = x0-x2;

      float dy0 = y1-y0;
      float dy1 = y2-y1;
      float dy2 = y0-y2;
      

      for (float x = minx; x < maxx; x++) {
          for (float y = miny; y < maxy; y++) {
              int sx = (int)floor(x);
              int sy = (int)floor(y);
              
              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;
              

              for (int subx = 0; subx < sqrt(sample_rate); subx++) {
                  for (int suby = 0; suby < sqrt(sample_rate); suby++) {
                      int idx_sub = suby * sqrt(sample_rate) + subx;
                      int idx = y*width*sample_rate + x*sample_rate + idx_sub;

                      //check bound

                      if (idx < 0 || idx >= width*height*sample_rate) {
                          continue; //we skip it
                      } else {
                          
                          float add_to_x = float(subx) * increment + increment_center;
                          float add_to_y = float(suby) * increment + increment_center;
                          
                          float y_prime = float(y) + add_to_y;
                          float x_prime = float(x) + add_to_x;
                          
//                          tuple<float,float,float> bary_coord = convert_to_barycentric(x, y, x0, y0, x1, y1, x2, y2);
            
                          tuple<float, float, float> bary_coord = convert_to_barycentric(x, y, x0, y0, x1, y1, x2, y2);
                          
                          float alpha = get<0>(bary_coord);
                          float beta = get<1>(bary_coord);
                          float gamma = get<2>(bary_coord);
                          
                          Color col = alpha * c0 + beta * c1 + gamma * c2;
                          

                          //check bound
                          if (x_prime < 0 || y_prime < 0 || x_prime >= width || y_prime >= height) {
                              continue;
                          }

                          if (point_in_triangle(x_prime, y_prime, x0, y0, x1, y1, x2, y2, dx0, dy0, dx1, dy1, dx2, dy2)) {
                              sample_buffer[idx] = col;
                          }

                      }
                  }
              }

          }
      }


  }

  tuple<float, float, float> RasterizerImp::convert_to_barycentric(float x, float y,
                               float xa, float ya,
                               float xb, float yb,
                               float xc, float yc) {
      float alpha_numerator = -(x - xb) * (yc - yb) + (y - yb) * (xc - xb);
      float beta_numerator = - (x - xc) * (ya - yc) + (y - yc) * (xa - xc);
      
      float alpha_denomenator = -(xa - xb) * (yc - yb) + (ya - yb) * (xc - xb);
      float beta_denomenator = - (xb - xc) * (ya - yc) + (yb - yc) * (xa - xc);

      float alpha = alpha_numerator / alpha_denomenator;
      float beta = beta_numerator / beta_denomenator;
      float gamma = float(1) - alpha - beta;
      
      tuple<float, float, float> tup;
      
      tup = make_tuple(alpha, beta, gamma);
      return tup;

  }




  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
      unsigned int sample_rate_uns = sample_rate;
      float increment = float(1)/float(sqrt(sample_rate_uns));
      float increment_center = increment/float(2);

    //locate the bounding boxes of the triangle
      float minx = floor(std::min({x0,x1,x2}));
      float miny = floor(std::min({y0,y1,y2}));

      float maxx = ceil(std::max({x0,x1,x2}));
      float maxy = ceil(std::max({y0,y1,y2}));

      float dx0 = x1-x0;
      float dx1 = x2-x1;
      float dx2 = x0-x2;

      float dy0 = y1-y0;
      float dy1 = y2-y1;
      float dy2 = y0-y2;
      SampleParams sample_params;

      for (float x = minx; x < maxx; x++) {
          for (float y = miny; y < maxy; y++) {
              int sx = (int)floor(x);
              int sy = (int)floor(y);

              if (sx < 0 || sx >= width) continue;
              if (sy < 0 || sy >= height) continue;


              for (int subx = 0; subx < sqrt(sample_rate); subx++) {
                  for (int suby = 0; suby < sqrt(sample_rate); suby++) {
                      int idx_sub = suby * sqrt(sample_rate) + subx;
                      int idx = y*width*sample_rate + x*sample_rate + idx_sub;

                      //check bound

                      if (idx < 0 || idx >= width*height*sample_rate) {
                          continue; //we skip it
                      } else {

                          float add_to_x = float(subx) * increment + increment_center;
                          float add_to_y = float(suby) * increment + increment_center;

                          float y_prime = float(y) + add_to_y;
                          float x_prime = float(x) + add_to_x;

//                          tuple<float,float,float> bary_coord = convert_to_barycentric(x, y, x0, y0, x1, y1, x2, y2);
//                          if (!point_in_triangle(x_prime, y_prime, x0, y0, x1, y1, x2, y2, dx0, dy0, dx1, dy1, dx2, dy2)) {
//                              continue;
//                          }


                          float l0 = -(x_prime - x0) * dy0 + (y_prime - y0) * dx0;
                          float l1 = -(x_prime - x1) * dy1 + (y_prime - y1) * dx1;
                          float l2 = -(x_prime - x2) * dy2 + (y_prime - y2) * dx2;

                          if (!((l0 >= 0) && (l1 >= 0) && (l2 >= 0))) {
                              if (!((l0 < 0) && (l1 < 0) && (l2 < 0))) {
                                  continue;
                              }
                          }

                          tuple<float, float, float> bary_coord = convert_to_barycentric(x_prime, y_prime, x0, y0, x1, y1, x2, y2);

                          float alpha = get<0>(bary_coord);
                          float beta = get<1>(bary_coord);
                          float gamma = get<2>(bary_coord);

                          Vector2D uv = Vector2D(u0,v0) * alpha + Vector2D(u1,v1) * beta + Vector2D(u2,v2) * gamma;
                          Color col;

                          if (psm == P_NEAREST) {
                              if (lsm == L_ZERO) {
                                  col = tex.sample_nearest(uv, 0);
                              } else {
                                  sample_params.p_uv = uv;
                                  tuple<float,float,float> bary_coord_dx = convert_to_barycentric(x_prime+float(1), y_prime, x0, y0, x1, y1, x2, y2);
                                  float alpha_dx = get<0>(bary_coord_dx);
                                  float beta_dx = get<1>(bary_coord_dx);
                                  float gamma_dx = get<2>(bary_coord_dx);
                                  sample_params.p_dx_uv = Vector2D(u0, v0) * alpha_dx + Vector2D(u1, v1) * beta_dx + Vector2D(u2, v2) * gamma_dx;
                                  

                                  tuple<float,float,float> bary_coord_dy = convert_to_barycentric(x_prime, y_prime+float(1), x0, y0, x1, y1, x2, y2);
                                  float alpha_dy = get<0>(bary_coord_dy);
                                  float beta_dy = get<1>(bary_coord_dy);
                                  float gamma_dy = get<2>(bary_coord_dy);
                                  sample_params.p_dy_uv = Vector2D(u0, v0) * alpha_dy + Vector2D(u1, v1) * beta_dy + Vector2D(u2, v2) * gamma_dy;

                                  sample_params.psm = psm;
                                  sample_params.lsm = lsm;
                                  col = tex.sample(sample_params);

                              }
                          } else if (psm == P_LINEAR) {
                              if (lsm == L_ZERO) {
                                  col = tex.sample_bilinear(uv, 0);
                              } else {
                                  col = tex.sample_nearest(uv, 0);
                              }
                          }
                          sample_buffer[idx] = col;

                      }
                  }
              }

          }
      }



  }


  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    clear_buffers();
    this->sample_buffer.resize(width * height*sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
      
    //aggregate color

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {

          int sx = (int)floor(x);
          int sy = (int)floor(y);

          // check bounds
          if (sx < 0 || sx >= width) continue;
          if (sy < 0 || sy >= height) continue;

          Color aggregate;
          Color col;
          for (int subx = 0; subx < sqrt(sample_rate); ++subx) {
              for (int suby = 0; suby < sqrt(sample_rate); ++suby) {
                  int idx_sub = suby * sqrt(sample_rate) + subx;
                  int idx = y*width*sample_rate + x*sample_rate + idx_sub;

                  //check bound

                  if (idx < 0 || idx >= width*height*sample_rate) {
                      continue; //we skip it
                  } else {
                      Color local_color = sample_buffer[idx];
                      aggregate.r += local_color.r;
                      aggregate.g += local_color.g;
                      aggregate.b += local_color.b;
                  }
              }
          }

          // now compute the average of the color
          col.r = aggregate.r / float(sample_rate);
          col.g = aggregate.g / float(sample_rate);
          col.b = aggregate.b / float(sample_rate);


        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }


  }

  Rasterizer::~Rasterizer() { }


}// CGL
