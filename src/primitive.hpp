#ifndef ALE_SOLVER_SRC_PRIMITIVE_HPP_
#define ALE_SOLVER_SRC_PRIMITIVE_HPP_

#include <memory>
#include <vector>

#include "logger.hpp"

enum class PrimitiveType { ePolygon, eRectangle, eCircle };
enum class ClockwiseDirection { eClockwise, eCounterclockwise, eNone };

/* Point is a dot with coordinates (x, y) */
class Point {
 public:
  Point(double x_, double y_);
  Point(const Point &rhs) = default;
  Point(Point &&rhs) = default;

  Point &operator=(const Point &rhs) = default;
  Point &operator=(Point &&rhs) = default;

  ~Point() = default;

  bool operator==(const Point &) const = default;

 public:
  double getX() const;
  double getY() const;

  void setX(double x_);
  void setY(double y_);

 private:
  double x;
  double y;
};

/* Line is a directional vector with start at first Point and end at the second
 */
class Line {
 public:
  Line(Point first_, Point second_);
  Line(const Line &rhs) = default;
  Line(Line &&rhs) = default;

  Line &operator=(const Line &rhs) = default;

 public:
  void reverse();  // Swap first and second point

  double getLength() const;  // Compute the length of the line

  Point getFirstPoint() const;
  Point getSecondPoint() const;

  void setFirstPoint(Point point);
  void setSecondPoint(Point point);

 private:
  Point first;
  Point second;
};

/* Edge is an array of continuous lines which will correspond to a particular
 * side of a quadrilateral grid */
class Edge {
 public:
  Edge(const std::vector<Line> &lines_);
  Edge(const Edge &rhs) = default;
  Edge(Edge &&rhs) = default;

  Edge &operator=(const Edge &rhs) = default;
  Edge &operator=(Edge &&rhs) = default;

 public:
  void divideLines(
      size_t newSize);  // if newSize is bigger than current number of lines
                        // divide them optimally to get newSize number of lines

  void reverse();  // Reverses edge

  bool isContinuous() const;  // returns if lines of the edge are connected

  const std::vector<Line> &getLines() const;
  std::vector<Line> &getLines();

 private:
  std::vector<Line> lines;
};

/* Polygon is an enclosed space which contains edges */
class Polygon {
 public:
  Polygon(std::shared_ptr<Edge> leftEdge, std::shared_ptr<Edge> topEdge,
          std::shared_ptr<Edge> rightEdge, std::shared_ptr<Edge> bottomEdge);

 public:
  virtual void generateMesh(std::vector<std::vector<double>> &x,
                            std::vector<std::vector<double>> &y)
      const;  // Generate a quadrilateral mesh

  bool isEnclosed()
      const;  // Check if edges are continuous and their ends are connected
  virtual bool isInside(
      double x_, double y_) const;  // Check if a point is inside the polygon
  bool isInside(const Point &point) const;
  bool isInside(const Line &line) const;
  bool isInside(std::shared_ptr<Edge> edge) const;

  bool isIntersecting(const Line &line) const;
  bool isIntersecting(const std::shared_ptr<Edge> edge) const;

  ClockwiseDirection getClockwiseDirection()
      const;  // returns if the edge structure is clockwise,
              // counterClockwise or the polygon is not enclosed
  void reverseClockwiseDirection();  // if the polygon is enclosed changes
                                     // clockwise direction

  std::shared_ptr<Edge> getLeftEdge() const;
  std::shared_ptr<Edge> getTopEdge() const;
  std::shared_ptr<Edge> getRightEdge() const;
  std::shared_ptr<Edge> getBottomEdge() const;
  PrimitiveType getType() const;

 private:
  bool isLineIntersecting(double x_, double y_, const Line &line) const;

 protected:
  std::shared_ptr<Edge> leftEdge;
  std::shared_ptr<Edge> topEdge;
  std::shared_ptr<Edge> rightEdge;
  std::shared_ptr<Edge> bottomEdge;

 private:
  PrimitiveType type;

  Logger logger;
};

/* rectangle is a polygon with four straight perpendicular edges */
class Rectangle : public Polygon {
 public:
  Rectangle(double x0, double y0, double width, double height);
  Rectangle(double x0, double y0, double width, double height, double phi);

 public:
  virtual void generateMesh(std::vector<std::vector<double>> &x,
                            std::vector<std::vector<double>> &y) const override;

  virtual bool isInside(double x_, double y_) const override;
};

/* Region2D is an array of polygons to support arbitrary shape */
class Region2D {
 public:
  Region2D(std::unique_ptr<Polygon> polygon);

 public:
  bool isInside(double x_,
                double y_) const;  // Check if a point is inside the region

  const std::vector<std::vector<std::unique_ptr<Polygon>>> &getPolygons() const;
  std::vector<std::vector<std::unique_ptr<Polygon>>> &getPolygons();

 protected:
  Logger logger;

 private:
  std::vector<std::vector<std::unique_ptr<Polygon>>> polygons;
};

#endif
