#ifndef ALE_SOLVER_METHOD_HPP_
#define ALE_SOLVER_METHOD_HPP_

enum class MethodType { eALE };

class Method {
 public:
  Method();
  Method(const Method &method);
  Method(Method &&method);

  Method &operator=(const Method &method);
  Method &operator=(Method &&method);

  ~Method();

 private:
  MethodType type;
};

#endif
