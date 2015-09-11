#ifndef PHYSICS_H
#define PHYSICS_H

#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>


#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace ph
{

const double EPS = 1E-9;
const double fEPS = 1E-6;

class intPolygon
{
public:
    intPolygon(std::vector< sf::Vector2f > point);
    intPolygon(sf::Rect < int > rect);
    sf::Vector2f CheckIntersect(sf::Vector2f p0, sf::Vector2f);
    sf::Vector2f ReflexSpeed(sf::Vector2f p0, sf::Vector2f * speedvec);
    float GetArea();
private:
    std::vector< sf::Vector2f > Points;
    float Area;
};

sf::Vector2f GetPointOfIntersection(sf::Vector2f a0, sf::Vector2f a1, sf::Vector2f b0, sf::Vector2f b1);
int sgn(float x);

}

namespace ph
{

inline float det (float a, float b, float c, float d)
{
    return a * d - b * c;
}

bool get_line_intersection(float p0_x, float p0_y, float p1_x, float p1_y,
    float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y)
{
    float s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    float s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != NULL)
            *i_x = p0_x + (t * s1_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s1_y);
        return true;
    }

    return false; // No collision
}


float sqr(float x){ return x*x; }
float length_squared(sf::Vector2f v, sf::Vector2f w){
    return sqr(w.x - v.x) + sqr(w.y - v.y);
}
float distance(sf::Vector2f p, sf::Vector2f v){
    return sqr(p.x - v.x) + sqr(p.y - v.y);
}
float dot(sf::Vector2f a, sf::Vector2f b){
    return (a.x*b.x + a.y*b.y);
}

float minimum_squared_distance(sf::Vector2f v, sf::Vector2f w, sf::Vector2f p) {
  // Return minimum distance between line segment vw and point p
  const float l2 = length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0) return distance(p, v);   // v == w case
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line.
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  const float t = dot(sf::Vector2f(p.x - v.x, p.y - v.y),
                      sf::Vector2f(w.x - v.x, w.y - v.y)) / l2;
  if (t < 0.0) return distance(p, v);       // Beyond the 'v' end of the segment
  else if (t > 1.0) return distance(p, w);  // Beyond the 'w' end of the segment
  const sf::Vector2f projection (v.x + t * (w.x - v.x),
                                 v.y + t * (w.y - v.y));  // Projection falls on the segment
  return distance(p, projection);
}

sf::Vector2f GetPointOfIntersection(sf::Vector2f a0, sf::Vector2f a1, sf::Vector2f b0, sf::Vector2f b1)
{
    sf::Vector2f a (a0.x, a0.y);
    sf::Vector2f b (a1.x, a1.y);
    sf::Vector2f c (b0.x, b0.y);
    sf::Vector2f d (b1.x, b1.y);
    float x,y;
    bool is_intersect =  get_line_intersection(a0.x, a0.y, a1.x, a1.y,
                                               b0.x, b0.y, b1.x, b1.y,
                                               &x, &y); //intersect(a,b,c,d,&x, &y);
    if (is_intersect) return sf::Vector2f((float) x, (float) y);
    return sf::Vector2f(-1.0f, -1.0f);
}
int sgn(float x){
    return (fabs(x) < fEPS)? 0.0f : ( (x>fEPS)? 1.0f : -1.0f );
}




intPolygon::intPolygon(std::vector< sf::Vector2f > point){
    this->Points = point;
    this->Area = 0;
    if (this->Points.size() == 0) return ;
    // find area by the 1-th (with index 0) point
    sf::Vector2f p0 = this->Points[0];
    for (int i=0; i<this->Points.size(); i++){

        int ind_1 = (i+0 % this->Points.size());
        int ind_2 = (i+1 % this->Points.size());
        sf::Vector2f v1 = Points[ind_1] - p0;
        sf::Vector2f v2 = Points[ind_2] - p0;
        Area += ph::det(v1.x, v1.y, v2.x, v2.y);
    }
}
intPolygon::intPolygon(sf::Rect<int> rect){
    this->Points.clear();
    this->Points.push_back(sf::Vector2f(rect.left,              rect.top              ));
    this->Points.push_back(sf::Vector2f(rect.left + rect.width, rect.top              ));
    this->Points.push_back(sf::Vector2f(rect.left + rect.width, rect.top + rect.height));
    this->Points.push_back(sf::Vector2f(rect.left             , rect.top + rect.height));
    this->Area = rect.height * rect.width;
}

float intPolygon::GetArea(){
    return this->Area;
}
sf::Vector2f intPolygon::CheckIntersect(sf::Vector2f p0, sf::Vector2f p1){
    sf::Vector2f iP = sf::Vector2f(-1.0f,-1.0f);

    if (p0 == p1) return iP;

    float distance = 1000000000.0f;
    for (int i=0; i<this->Points.size(); i++){
        int ind_1 = ((i+0) % this->Points.size());
        int ind_2 = ((i+1) % this->Points.size());
        sf::Vector2f v1 = Points[ind_1];
        sf::Vector2f v2 = Points[ind_2];
        sf::Vector2f temp = ph::GetPointOfIntersection(p0, p1, v1, v2);
        if (temp == sf::Vector2f(-1.0f,-1.0f)) continue;

        //sf::Vector2f dis(temp.x - p0.x, temp.y - p0.y);
        float tmpDis = ph::minimum_squared_distance(sf::Vector2f(v1.x, v1.y),
                                                     sf::Vector2f(v2.x, v2.y),
                                                     sf::Vector2f(p0.x, p0.y)); //((dis.x * dis.x) + (dis.y * dis.y));
        if (tmpDis < distance + ph::fEPS){
            distance = tmpDis;
            iP = sf::Vector2f(temp.x, temp.y);
        }
    }
    //if (iP != sf::Vector2i(-1, -1))
    //    std::cout << std::endl << "intersect: " << iP.x << " " << iP.y << std::endl << std::endl;
    return iP;
}
sf::Vector2f intPolygon::ReflexSpeed(sf::Vector2f p0, sf::Vector2f *speedvec){
    sf::Vector2f iP = sf::Vector2f(-1.0f,-1.0f);
    float distance = 1000000000.0f;
    sf::Vector2f pDest (p0.x + speedvec->x,
                        p0.y + speedvec->y);
    sf::Vector2f p1 = sf::Vector2f(-1.0f,-1.0f);
    sf::Vector2f p2 = sf::Vector2f(-1.0f,-1.0f);

    for (int i=0; i<this->Points.size(); i++){
        int ind_1 = ((i+0) % this->Points.size());
        int ind_2 = ((i+1) % this->Points.size());
        sf::Vector2f v1 = Points[ind_1];
        sf::Vector2f v2 = Points[ind_2];
        sf::Vector2f temp = ph::GetPointOfIntersection(p0, pDest, v1, v2);
        if ((temp == sf::Vector2f(-1.0f,-1.0f))) continue;
                //|| (temp == sf::Vector2f(p0.x, p0.y))) continue;

        //sf::Vector2f dis(temp.x - p0.x, temp.y - p0.y);
        float tmpDis = ph::minimum_squared_distance(sf::Vector2f(v1.x, v1.y),
                                                     sf::Vector2f(v2.x, v2.y),
                                                     sf::Vector2f(p0.x, p0.y)); //((dis.x * dis.x) + (dis.y * dis.y)); // ((dis.x * dis.x) + (dis.y * dis.y));
        if (tmpDis < distance + ph::EPS){
            distance = tmpDis;
            iP = sf::Vector2f(temp.x, temp.y);
            p1 = v1;
            p2 = v2;
        }
    }

    if (iP == sf::Vector2f(-1.0f, -1.0f)) return iP; // no intersection point - nothing to do!
    //std::cout << "bounce line: " << p1.x << " " << p1.y << " ; " << p2.x << " " << p2.y << std::endl;(it's new)


    sf::Vector2f f (p2.x - p1.x, p2.y - p1.y);
    sf::Vector2f n (p2.y - p1.y, p1.x - p2.x);
    sf::Vector2f v = *speedvec;

    float len = sqrt(n.x * n.x + n.y * n.y);
    n.x /= len;
    n.y /= len;

    float dot2 = 2*(n.x * v.x + n.y * v.y);
    sf::Vector2f v_ref (v.x - dot2*n.x, v.y - dot2*n.y);
    *speedvec = sf::Vector2f(v_ref.x, v_ref.y);

    return iP;
}



};


#endif // PHYSICS_H

