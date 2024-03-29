#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  
    
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan;
    
    // store transform
    double mat[] = {1/(2*vspan),           0, -(centerX-vspan)/(2*vspan),
                    0, 1/(2*vspan), -(centerY-vspan)/(2*vspan),
                    0,           0,                          1};
    Matrix3x3 m = Matrix3x3(mat);
    set_svg_2_norm(m);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;

  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
