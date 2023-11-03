#include "primitive.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>

Point::Point(double x_, double y_) : x(x_), y(y_) {}
double Point::getX() const { return x; }
double Point::getY() const { return y; }
void Point::setX(double x_) { x = x_; }
void Point::setY(double y_) { y = y_; }

Line::Line(Point first_, Point second_) : first(first_), second(second_) {}
void Line::reverse() { std::swap(first, second); }
double Line::getLength() const {
  double dx = second.getX() - first.getX();
  double dy = second.getY() - first.getY();
  return std::sqrt(dx * dx + dy * dy);
}
Point Line::getFirstPoint() const { return first; }
Point Line::getSecondPoint() const { return second; }
void Line::setFirstPoint(Point point) { first = point; }
void Line::setSecondPoint(Point point) { second = point; }

Edge::Edge(const std::vector<Line> &lines_) : lines(lines_) {}
void Edge::divideLines(size_t newSize) {
  auto currentSize = lines.size();
  if (newSize > currentSize) {
    auto divisions = newSize - currentSize;

    std::vector<int> divPerLine(currentSize, 0);
    for (int i = 0; i < divisions; i++) {
      double maxLength = lines[0].getLength() / (divPerLine[0] + 1);
      int maxIndex = 0;
      for (int j = 0; j < lines.size(); j++) {
        double length = lines[j].getLength() / (divPerLine[j] + 1);
        if (length > maxLength) {
          maxLength = length;
          maxIndex = j;
        }
      }
      divPerLine[maxIndex]++;
    }

    lines.reserve(newSize);
    int icr = 0;
    for (int i = 0; i < divPerLine.size(); i++) {
      if (divPerLine[i] == 0) {
        continue;
      }
      auto &divLine = lines[i + icr];
      Point firstPoint = divLine.getFirstPoint();
      Point secondPoint = divLine.getSecondPoint();
      double dx =
          (secondPoint.getX() - firstPoint.getX()) / (divPerLine[i] + 1);
      double dy =
          (secondPoint.getY() - firstPoint.getY()) / (divPerLine[i] + 1);

      divLine = Line(Point(firstPoint.getX() + divPerLine[i] * dx,
                           firstPoint.getY() + divPerLine[i] * dy),
                     secondPoint);

      std::vector<Line> newLines{};
      newLines.reserve(divPerLine[i]);
      for (int j = 0; j < divPerLine[i]; j++) {
        newLines.emplace_back(
            Point(firstPoint.getX() + j * dx, firstPoint.getY() + j * dy),
            Point(firstPoint.getX() + (j + 1) * dx,
                  firstPoint.getY() + (j + 1) * dy));
      }
      lines.insert(lines.begin() + i + icr, newLines.begin(), newLines.end());
      icr += divPerLine[i];
    }
  }
}
void Edge::reverse() {
  for (auto &line : lines) {
    line.reverse();
  }
  std::reverse(lines.begin(), lines.end());
}
bool Edge::isContinuous() const {
  for (int i = 1; i < lines.size(); i++) {
    if (lines[i].getFirstPoint() != lines[i - 1].getSecondPoint()) {
      return false;
    }
  }
  return true;
}
const std::vector<Line> &Edge::getLines() const { return lines; }
std::vector<Line> &Edge::getLines() { return lines; }

Polygon::Polygon(std::shared_ptr<Edge> leftEdge_,
                 std::shared_ptr<Edge> topEdge_,
                 std::shared_ptr<Edge> rightEdge_,
                 std::shared_ptr<Edge> bottomEdge_)
    : leftEdge(leftEdge_),
      topEdge(topEdge_),
      rightEdge(rightEdge_),
      bottomEdge(bottomEdge_) {}
bool Polygon::isEnclosed() const {
  return getClockwiseDirection() != ClockwiseDirection::eNone;
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
bool Polygon::isInside(const Point &point) const {
  return isInside(point.getX(), point.getY());
}
bool Polygon::isInside(const Line &line) const {
  return isInside(line.getFirstPoint()) && isInside(line.getSecondPoint());
}
bool Polygon::isInside(std::shared_ptr<Edge> edge) const {
  for (const auto &line : edge->getLines()) {
    if (!(isInside(line.getSecondPoint()))) {
      return false;
    }
  }
  return isInside(edge->getLines()[0].getFirstPoint());
}
bool Polygon::isIntersecting(const Line &line) const {
  return isInside(line.getFirstPoint()) != isInside(line.getSecondPoint());
}
bool Polygon::isIntersecting(const std::shared_ptr<Edge> edge) const {
  for (const auto &line : edge->getLines()) {
    if (isIntersecting(line)) {
      return true;
    }
  }
  return false;
}
ClockwiseDirection Polygon::getClockwiseDirection() const {
  if (!(leftEdge->isContinuous() && topEdge->isContinuous() &&
        rightEdge->isContinuous() && bottomEdge->isContinuous())) {
    return ClockwiseDirection::eNone;
  }

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
  if (bClockwise) {
    return ClockwiseDirection::eClockwise;
  }
  bool bCounterClockwise =
      leftEdgeLast == bottomEdgeFirst && bottomEdgeLast == rightEdgeFirst &&
      rightEdgeLast == topEdgeFirst && topEdgeLast == leftEdgeFirst;
  if (bCounterClockwise) {
    return ClockwiseDirection::eCounterclockwise;
  }

  return ClockwiseDirection::eNone;
}
void Polygon::reverseClockwiseDirection() {
  leftEdge->reverse();
  topEdge->reverse();
  rightEdge->reverse();
  bottomEdge->reverse();
}
std::shared_ptr<Edge> Polygon::getLeftEdge() const { return leftEdge; }
std::shared_ptr<Edge> Polygon::getTopEdge() const { return topEdge; }
std::shared_ptr<Edge> Polygon::getRightEdge() const { return rightEdge; }
std::shared_ptr<Edge> Polygon::getBottomEdge() const { return bottomEdge; }
PrimitiveType Polygon::getType() const { return type; }
/* Check if line is intersecting with the ray coming from (x_, y_) in Ox
 * direction */
bool Polygon::isLineIntersecting(double x_, double y_, const Line &line) const {
  double x1 = line.getFirstPoint().getX();
  double y1 = line.getFirstPoint().getY();
  double x2 = line.getSecondPoint().getX();
  double y2 = line.getSecondPoint().getY();

  if (y2 == y1 && y_ == y1) {
    return true;
  }

  if (!((y_ < y1) != (y_ < y2))) {
    return false;
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
Rectangle::Rectangle(double x0, double y0, double width, double height,
                     double phi)
    : Polygon(nullptr, nullptr, nullptr, nullptr) {
  using namespace std::numbers;
  phi = std::fmod(phi, 2.0 * pi);
  double xbl;
  double ybl;
  double xbr;
  double ybr;
  double xtl;
  double ytl;
  double xtr;
  double ytr;

  if (phi < 0.25 * pi || phi >= 1.75 * pi) {
    xbl = x0;
    ybl = y0;
    xbr = std::cos(phi) * width + xbl;
    ybr = std::sin(phi) * width + ybl;
    xtl = std::cos(phi + 0.5 * pi) * height + xbl;
    ytl = std::sin(phi + 0.5 * pi) * height + ybl;
    xtr = std::cos(phi) * width + xtl;
    ytr = std::sin(phi) * width + ytl;
  } else if (phi < 0.75 * pi) {
    xbr = x0;
    ybr = y0;
    xtr = std::cos(phi) * width + xbr;
    ytr = std::sin(phi) * width + ybr;
    xbl = std::cos(phi + 0.5 * pi) * height + xbr;
    ybl = std::sin(phi + 0.5 * pi) * height + ybr;
    xtl = std::cos(phi) * width + xbl;
    ytl = std::sin(phi) * width + ybl;
  } else if (phi < 1.25 * pi) {
    xtr = x0;
    ytr = y0;
    xtl = std::cos(phi) * width + xtr;
    ytl = std::sin(phi) * width + ytr;
    xbr = std::cos(phi + 0.5 * pi) * height + xtr;
    ybr = std::sin(phi + 0.5 * pi) * height + ytr;
    xbl = std::cos(phi) * width + xbr;
    ybl = std::sin(phi) * width + ybr;
  } else if (phi < 1.75 * pi) {
    xtl = x0;
    ytl = y0;
    xbl = std::cos(phi) * width + xtl;
    ybl = std::sin(phi) * width + ytl;
    xtr = std::cos(phi + 0.5 * pi) * height + xtl;
    ytr = std::sin(phi + 0.5 * pi) * height + ytl;
    xbr = std::cos(phi) * width + xtr;
    ybr = std::sin(phi) * width + ytr;
  }

  Point bottomLeft(xbl, ybl);
  Point topLeft(xtl, ytl);
  Point topRight(xtr, ytr);
  Point bottomRight(xbr, ybr);
  Line left{bottomLeft, topLeft};
  Line top{topLeft, topRight};
  Line right{topRight, bottomRight};
  Line bottom{bottomRight, bottomLeft};
  leftEdge = std::make_shared<Edge>(std::vector<Line>{left});
  topEdge = std::make_shared<Edge>(std::vector<Line>{top});
  rightEdge = std::make_shared<Edge>(std::vector<Line>{right});
  bottomEdge = std::make_shared<Edge>(std::vector<Line>{bottom});
}
void Rectangle::generateMesh(std::vector<std::vector<double>> &x,
                             std::vector<std::vector<double>> &y) const {
  // TEMPORARY VALUES TODO remove later
  size_t sizeX = 100;
  size_t sizeY = 100;

  double xbl = leftEdge->getLines()[0].getFirstPoint().getX();
  double ybl = leftEdge->getLines()[0].getFirstPoint().getY();
  double xtl = topEdge->getLines()[0].getFirstPoint().getX();
  double ytl = topEdge->getLines()[0].getFirstPoint().getY();
  double xtr = rightEdge->getLines()[0].getFirstPoint().getX();
  double ytr = rightEdge->getLines()[0].getFirstPoint().getY();

  x.clear();
  y.clear();
  x.reserve(sizeX + 1);
  y.reserve(sizeX + 1);
  double dxi = (xtr - xbl) / sizeX;
  double dyi = (ytr - ybl) / sizeY;
  double dxj = (xtl - xbl) / sizeX;
  double dyj = (ytl - ybl) / sizeY;
  for (int i = 0; i < sizeX + 1; i++) {
    std::vector<double> tmpX;
    std::vector<double> tmpY;
    tmpX.reserve(sizeY + 1);
    tmpY.reserve(sizeY + 1);
    for (int j = 0; j < sizeY + 1; j++) {
      tmpX.push_back(xbl + dxi * i + dxj * j);
      tmpY.push_back(ybl + dyi * i + dyj * j);
    }
    x.push_back(tmpX);
    y.push_back(tmpY);
  }
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

Region2D::Region2D(std::unique_ptr<Polygon> polygon) {
  if (!(polygon->isEnclosed())) {
    logger.Log("class Region2D: polygon is not enclosed", LogLevel::eWarning);
  }
  std::vector<std::unique_ptr<Polygon>> polygonsTmp{};
  polygonsTmp.push_back(std::move(polygon));
  polygons.push_back(std::move(polygonsTmp));
}
bool Region2D::isInside(double x_, double y_) const {
  for (const auto &polygonLine : polygons) {
    for (const auto &polygon : polygonLine) {
      if (polygon->isInside(x_, y_)) {
        return true;
      }
    }
  }
  return false;
}
const std::vector<std::vector<std::unique_ptr<Polygon>>> &
Region2D::getPolygons() const {
  return polygons;
}
std::vector<std::vector<std::unique_ptr<Polygon>>> &Region2D::getPolygons() {
  return polygons;
}
