#include "primitive.hpp"

#include <cassert>
#include <cmath>
#include <numbers>

Point::Point(double x_, double y_) : x(x_), y(y_) {}
double Point::getX() const { return x; }
double Point::getY() const { return y; }
void Point::setX(double x_) { x = x_; }
void Point::setY(double y_) { y = y_; }

Line::Line(std::shared_ptr<Point> first_, std::shared_ptr<Point> second_)
    : first(std::move(first_)), second(std::move(second_)) {}
std::shared_ptr<Point> Line::getFirstPoint() const { return first; }
std::shared_ptr<Point> Line::getSecondPoint() const { return second; }
void Line::setFirstPoint(std::shared_ptr<Point> point) { first = point; }
void Line::setSecondPoint(std::shared_ptr<Point> point) { second = point; }
double Line::getLength() const {
  double dx = second->getX() - first->getX();
  double dy = second->getY() - first->getY();
  return std::sqrt(dx * dx + dy * dy);
}

Edge::Edge(const std::vector<Line> &lines_) : lines(lines_) {}
const std::vector<Line> &Edge::getLines() const { return lines; }
std::vector<Line> &Edge::getLines() { return lines; }
bool Edge::isContinuous() const {
  for (int i = 1; i < lines.size(); i++) {
    if (lines[i].getFirstPoint() != lines[i - 1].getSecondPoint()) {
      return false;
    }
  }
  return true;
}

// PrimitiveRectangle::PrimitiveRectangle(double x0_, double y0_, double width_,
//                                        double height_, double phi)
//     : x{0.0, 0.0, 0.0, 0.0}, y{0.0, 0.0, 0.0, 0.0} {
//   double pi = std::numbers::pi;
//   phi = std::fmod(phi, 2 * pi);
//
//   if (phi < 0.25 * pi || phi >= 3.75 * pi) {
//     x[0] = x0_;
//     y[0] = y0_;
//
//     x[3] = std::cos(phi) * width_ + x[0];
//     y[3] = std::sin(phi) * width_ + y[0];
//
//     x[1] = std::cos(phi + 0.5 * pi) * height_ + x[0];
//     y[1] = std::sin(phi + 0.5 * pi) * height_ + y[0];
//
//     x[2] = std::cos(phi) * width_ + x[1];
//     y[2] = std::sin(phi) * width_ + y[1];
//
//   } else if (phi < 0.75 * pi) {
//     x[3] = x0_;
//     y[3] = y0_;
//
//     x[2] = std::cos(phi) * width_ + x[3];
//     y[2] = std::sin(phi) * width_ + y[3];
//
//     x[0] = std::cos(phi + 0.5 * pi) * height_ + x[3];
//     y[0] = std::sin(phi + 0.5 * pi) * height_ + y[3];
//
//     x[1] = std::cos(phi) * width_ + x[0];
//     y[1] = std::sin(phi) * width_ + y[0];
//   } else if (phi < 3.25 * pi) {
//     x[2] = x0_;
//     y[2] = y0_;
//
//     x[1] = std::cos(phi) * width_ + x[2];
//     y[1] = std::sin(phi) * width_ + y[2];
//
//     x[3] = std::cos(phi + 0.5 * pi) * height_ + x[2];
//     y[3] = std::sin(phi + 0.5 * pi) * height_ + y[2];
//
//     x[0] = std::cos(phi) * width_ + x[3];
//     y[0] = std::sin(phi) * width_ + y[3];
//   } else if (phi < 3.75 * pi) {
//     x[1] = x0_;
//     y[1] = y0_;
//
//     x[0] = std::cos(phi) * width_ + x[1];
//     y[0] = std::sin(phi) * width_ + y[1];
//
//     x[2] = std::cos(phi + 0.5 * pi) * height_ + x[1];
//     y[2] = std::sin(phi + 0.5 * pi) * height_ + y[1];
//
//     x[3] = std::cos(phi) * width_ + x[2];
//     y[3] = std::sin(phi) * width_ + y[2];
//   }
// }
//
// bool PrimitiveRectangle::isInside(double x_, double y_) const {
//   double dot3 = (x_ - x[0]) * (x[3] - x[0]) + (y_ - y[0]) * (y[3] - y[0]);
//   if (dot3 < 0 ||
//       dot3 > (x[3] - x[0]) * (x[3] - x[0]) + (y[3] - y[0]) * (y[3] - y[0])) {
//     return false;
//   }
//
//   double dot2 = (x_ - x[0]) * (x[1] - x[0]) + (y_ - y[0]) * (y[1] - y[0]);
//   if (dot2 < 0 ||
//       dot2 > (x[1] - x[0]) * (x[1] - x[0]) + (y[1] - y[0]) * (y[1] - y[0])) {
//     return false;
//   }
//
//   return true;
// }
//
// void PrimitiveRectangle::getVertices(std::array<double, 4> &x,
//                                      std::array<double, 4> &y) const {
//   x = this->x;
//   y = this->y;
// }
//
// void PrimitiveRectangle::setBoundaries(const std::vector<BoundaryType> &bnds)
// {
//   assert(bnds.size() == 4);
//   for (int i = 0; i < bnds.size(); i++) {
//     boundaries[i] = bnds[i];
//   }
// }
//
// std::vector<BoundaryType> PrimitiveRectangle::getBoundaries() const {
//   return std::vector<BoundaryType>(boundaries.cbegin(), boundaries.cend());
// }
//
// PrimitiveCircle::PrimitiveCircle(double x0_, double y0_, double r_)
//     : x0(x0_), y0(y0_), r(r_) {}
//
// bool PrimitiveCircle::isInside(double x_, double y_) const {
//   return (x_ - x0) * (x_ - x0) + (y_ - y0) * (y_ - y0) <= r * r;
// }
