#include "primitive.hpp"

#include <cassert>
#include <cmath>
#include <numbers>

Point::Point(double x_, double y_) : x(x_), y(y_) {}
double Point::getX() const { return x; }
double Point::getY() const { return y; }
void Point::setX(double x_) { x = x_; }
void Point::setY(double y_) { y = y_; }

Line::Line(Point first_, Point second_) : first(first_), second(second_) {}
Point Line::getFirstPoint() const { return first; }
Point Line::getSecondPoint() const { return second; }
void Line::setFirstPoint(Point point) { first = point; }
void Line::setSecondPoint(Point point) { second = point; }
double Line::getLength() const {
  double dx = second.getX() - first.getX();
  double dy = second.getY() - first.getY();
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

Polygon::Polygon(std::shared_ptr<Edge> leftEdge_,
                 std::shared_ptr<Edge> topEdge_,
                 std::shared_ptr<Edge> rightEdge_,
                 std::shared_ptr<Edge> bottomEdge_)
    : leftEdge(leftEdge_),
      topEdge(topEdge_),
      rightEdge(rightEdge_),
      bottomEdge(bottomEdge_) {}
bool Polygon::isEnclosed() const {
  if (!(leftEdge->isContinuous() && topEdge->isContinuous() &&
        rightEdge->isContinuous() && bottomEdge->isContinuous())) {
    return false;
  }

  // Get first and last points of the edges
  auto leftEdgeFirst = leftEdge->getLines()[0].getFirstPoint();
  auto leftEdgeLast =
      leftEdge->getLines()[leftEdge->getLines().size() - 1].getSecondPoint();

  auto topEdgeFirst = topEdge->getLines()[0].getFirstPoint();
  auto topEdgeLast =
      topEdge->getLines()[topEdge->getLines().size() - 1].getSecondPoint();

  auto rightEdgeFirst = rightEdge->getLines()[0].getFirstPoint();
  auto rightEdgeLast =
      rightEdge->getLines()[rightEdge->getLines().size() - 1].getSecondPoint();

  auto bottomEdgeFirst = bottomEdge->getLines()[0].getFirstPoint();
  auto bottomEdgeLast =
      bottomEdge->getLines()[bottomEdge->getLines().size() - 1]
          .getSecondPoint();

  // Check if geometry constructed clockwise or counterclockwise
  bool bClockwise =
      leftEdgeLast == topEdgeFirst && topEdgeLast == rightEdgeFirst &&
      rightEdgeLast == bottomEdgeFirst && bottomEdgeLast == leftEdgeFirst;
  bool bCounterClockwise =
      leftEdgeLast == bottomEdgeFirst && bottomEdgeLast == rightEdgeFirst &&
      rightEdgeLast == topEdgeFirst && topEdgeLast == leftEdgeFirst;

  return bClockwise || bCounterClockwise;
}

/* Checks if the point with coordinates (x_, y_) is inside the polygon */
bool Polygon::isInside(double x_, double y_) const {
  // Cast a ray from point to the right.
  // If the ray intersects with odd number of lines - the point is inside

  size_t counter = 0;

  for (const auto &line : leftEdge->getLines()) {
    if (isLineIntersecting(x_, y_, line)) {
      counter++;
    }
  }
  for (const auto &line : topEdge->getLines()) {
    if (isLineIntersecting(x_, y_, line)) {
      counter++;
    }
  }
  for (const auto &line : rightEdge->getLines()) {
    if (isLineIntersecting(x_, y_, line)) {
      counter++;
    }
  }
  for (const auto &line : bottomEdge->getLines()) {
    if (isLineIntersecting(x_, y_, line)) {
      counter++;
    }
  }

  return counter % 2 == 1;
}
/* Check if line is intersecting with the ray coming from (x_, y_) in Ox
 * direction */
bool Polygon::isLineIntersecting(double x_, double y_, const Line &line) const {
  double x1 = line.getFirstPoint().getX();
  double y1 = line.getFirstPoint().getY();
  double x2 = line.getSecondPoint().getX();
  double y2 = line.getSecondPoint().getY();

  if (!((y_ < y1) != (y_ < y2))) {
    return false;
  }

  if (y2 == y1 && y_ == y1) {
    return true;
  }

  if (x_ < x1 + ((y_ - y1) / (y2 - y1)) * (x2 - x1)) {
    return true;
  }
  return false;
}

Rectangle::Rectangle(double x0, double y0, double width, double height)
    : Polygon(nullptr, nullptr, nullptr, nullptr) {
  Point bottomLeft{x0, y0};
  Point topLeft{x0, y0 + height};
  Point topRight{x0 + width, y0 + height};
  Point bottomRight{x0 + width, y0};
  Line left{bottomLeft, topLeft};
  Line top{topLeft, topRight};
  Line right{topRight, bottomRight};
  Line bottom{bottomRight, bottomLeft};
  leftEdge = std::make_shared<Edge>(std::vector<Line>{left});
  topEdge = std::make_shared<Edge>(std::vector<Line>{top});
  rightEdge = std::make_shared<Edge>(std::vector<Line>{right});
  bottomEdge = std::make_shared<Edge>(std::vector<Line>{bottom});
}
/* Checks if point is inside the rectangle by calculating two scalar products */
bool Rectangle::isInside(double x_, double y_) const {
  double xbl =
      leftEdge->getLines()[0].getFirstPoint().getX();  // bottom left vert
  double ybl = leftEdge->getLines()[0].getFirstPoint().getY();
  double xbr =
      bottomEdge->getLines()[0].getFirstPoint().getX();  // bottom right vert
  double ybr = bottomEdge->getLines()[0].getFirstPoint().getY();
  double xtl = topEdge->getLines()[0].getFirstPoint().getX();  // top left vert
  double ytl = leftEdge->getLines()[0].getFirstPoint().getY();
  double dot1 = (x_ - xbl) * (xbr - xbl) + (y_ - ybl) * (ybr - ybl);

  if (dot1 < 0 ||
      dot1 > (xbr - xbl) * (xbr - xbl) + (ybr - ybl) * (ybr - ybl)) {
    return false;
  }
  double dot2 = (x_ - xbl) * (xtl - xbl) + (y_ - ybl) * (ytl - ybl);
  if (dot2 < 0 ||
      dot2 > (xtl - xbl) * (xtl - xbl) + (ytl - ybl) * (ytl - ybl)) {
    return false;
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
