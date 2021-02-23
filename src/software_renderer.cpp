#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  
    
    // update sampling buffer size
    buff_w = target_w * sample_rate;
    buff_h = target_h * sample_rate;
    sample_buffer = std::vector<unsigned char>(buff_w * buff_h * 4,255);
    
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

    // update smpling buffer size
    sample_buffer.clear();
    buff_w = target_w * sample_rate;
    buff_h = target_h * sample_rate;
    sample_buffer = std::vector<unsigned char>(buff_w * buff_h * 4,255);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x * sample_rate);
  int sy = (int) floor(y * sample_rate);

  // check bounds
  if ( sx < 0 || sx >= buff_w ) return;
  if ( sy < 0 || sy >= buff_h ) return;

  // fill sample - NOT doing alpha blending!
    for(size_t i=0; i<sample_rate; i++){
        for(size_t j=0; j<sample_rate; j++){
            sample_buffer[4 * ((sx+i) + (sy+j) * buff_w)    ] = (uint8_t) (color.r * 255);
            sample_buffer[4 * ((sx+i) + (sy+j) * buff_w) + 1] = (uint8_t) (color.g * 255);
            sample_buffer[4 * ((sx+i) + (sy+j) * buff_w) + 2] = (uint8_t) (color.b * 255);
            sample_buffer[4 * ((sx+i) + (sy+j) * buff_w) + 3] = (uint8_t) (color.a * 255);
        }
    }
    

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
    
    // round to integer
    int sx0 = (int) floor(x0 * sample_rate);
    int sy0 = (int) floor(y0 * sample_rate);
    int sx1 = (int) floor(x1 * sample_rate);
    int sy1 = (int) floor(y1 * sample_rate);
    
    // check bounds
    if ( sx0 < 0 || sx0 >= buff_w ) return;
    if ( sy0 < 0 || sy0 >= buff_h ) return;
    if ( sx1 < 0 || sx1 >= buff_w ) return;
    if ( sy1 < 0 || sy1 >= buff_h ) return;
    
    // Bresenham algorithm
    int dx = sx1 - sx0;
    int dy = sy1 - sy0;
    
    // vertical or horizontal
    if(dx == 0){
        // vertical
        int x = sx0;
        for(int y = sy0; y != sy1; y += dy/abs(dy)){
            sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
            sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
            sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
            sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
        }
        return;
    }else if(dy == 0){
        // horizontal
        int y = sy0;
        for(int x = sx0; x != sx1; x += dx/abs(dx)){
            sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
            sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
            sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
            sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
        }
        return;
    }
    
    // nonzero and finite slope
    if(dy*dx > 0){
        // positive slope
        if(sy1<sy0 && sx1<sx0){
            // make (sx1,sy1) always bigger
            int tempY = sy1;
            sy1 = sy0;
            sy0 = tempY;
            int tempX = sx1;
            sx1 = sx0;
            sx0 = tempX;
        }
        if(abs(dy) <= abs(dx)){
            // slope < 0.5
            int y = sy0;
            int eps = 0;
            for ( int x = sx0; x <= sx1; x++ )  {
                sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
                sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
                sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
                sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
                eps += abs(dy);

                if (2*eps >= abs(dx))  {
                    y++;
                    eps -= abs(dx);
                }
            }
        }else{
            // slope > 0.5
            int x = sx0;
            int eps = 0;
            for ( int y = sy0; y <= sy1; y++ )  {
                sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
                sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
                sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
                sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
                eps += abs(dx);

                if (2*eps >= abs(dy))  {
                    x++;
                    eps -= abs(dy);
                }
            }
        }
    }else{
        // negative slope
        if(sy1>sy0 && sx1<sx0){
            // make sy0<sy1, sx0>sx1
            int tempY = sy1;
            sy1 = sy0;
            sy0 = tempY;
            int tempX = sx1;
            sx1 = sx0;
            sx0 = tempX;
        }
        if(abs(dy) < abs(dx)){
            // slope > -0.5
            int y = sy0;
            int eps = 0;
            for ( int x = sx0; x <= sx1; x++ )  {
                sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
                sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
                sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
                sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
                eps += abs(dy);

                if (2*eps >= abs(dx))  {
                    y--;
                    eps -= abs(dx);
                }
            }
        }else{
            // slope < -0.5
            int x = sx0;
            int eps = 0;
            for ( int y = sy0; y >= sy1; y-- )  {
                sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
                sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
                sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
                sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
                eps += abs(dx);

                if (2*eps >= abs(dy))  {
                    x++;
                    eps -= abs(dy);
                }
            }
        }
        
    }
    return;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization
    
    // scale up by sampling rate
    x0 *= sample_rate; y0 *= sample_rate;
    x1 *= sample_rate; y1 *= sample_rate;
    x2 *= sample_rate; y2 *= sample_rate;
        
    // check bounds
    if ( x0 < 0 || x0 >= buff_w ) return;
    if ( y0 < 0 || y0 >= buff_h ) return;
    if ( x1 < 0 || x1 >= buff_w ) return;
    if ( y1 < 0 || y1 >= buff_h ) return;
    if ( x2 < 0 || x2 >= buff_w ) return;
    if ( y2 < 0 || y2 >= buff_h ) return;
     
    int xMax = max(x0,max(x1,x2));
    int xMin = min(x0,min(x1,x2));
    int yMax = max(y0,max(y1,y2));
    int yMin = min(y0,min(y1,y2));
    
    // edge as vector
    
    Vector2D v01(x1 - x0, y1 - y0); // vector: 0->1
    Vector2D v12(x2 - x1, y2 - y1);
    Vector2D v20(x0 - x2, y0 - y2);
    
    
    for(int y = yMin; y <= yMax; y++){
        for(int x = xMin; x <= xMax; x++){
            // cross product
            
            Vector2D v0p(x - x0, y - y0); // vector 0->p
            Vector2D v1p(x - x1, y - y1);
            Vector2D v2p(x - x2, y - y2);
            
            double cross0 = cross(v0p, v01);
            double cross1 = cross(v1p, v12);
            double cross2 = cross(v2p, v20);
            
            if(cross0*cross1>0 && cross0*cross2>0){
                sample_buffer[4 * (x + y * buff_w)    ] = (uint8_t) (color.r * 255);
                sample_buffer[4 * (x + y * buff_w) + 1] = (uint8_t) (color.g * 255);
                sample_buffer[4 * (x + y * buff_w) + 2] = (uint8_t) (color.b * 255);
                sample_buffer[4 * (x + y * buff_w) + 3] = (uint8_t) (color.a * 255);
            }
        }
    }
    
    
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
      
    int block_size = sample_rate * sample_rate;
    
    for(size_t y=0; y<target_h; y++){
        for(size_t x=0; x<target_w; x++){
            // box filter
            int avgR(0), avgG(0), avgB(0);
            for(int yy=0; yy<sample_rate; yy++){
                for(int xx=0; xx<sample_rate; xx++){
                    avgR += sample_buffer[4 * ((x*sample_rate+xx) + (y*sample_rate+yy) * buff_w)    ];
                    avgG += sample_buffer[4 * ((x*sample_rate+xx) + (y*sample_rate+yy) * buff_w) + 1];
                    avgB += sample_buffer[4 * ((x*sample_rate+xx) + (y*sample_rate+yy) * buff_w) + 2];
                    
                }
            }
            
            int originVal = (render_target[4 * (x + y * target_w)    ]
                            + render_target[4 * (x + y * target_w) + 1]
                            + render_target[4 * (x + y * target_w) + 2]) / 3;
            
                render_target[4 * (x + y * target_w)    ] = avgR / block_size;
                render_target[4 * (x + y * target_w) + 1] = avgG / block_size;
                render_target[4 * (x + y * target_w) + 2] = avgB / block_size;
                render_target[4 * (x + y * target_w) + 3] = 255;
        }
    }
    
    // clean buffer
    sample_buffer = std::vector<unsigned char>(buff_w * buff_h * 4,255);
    
    return;

}


} // namespace CMU462
