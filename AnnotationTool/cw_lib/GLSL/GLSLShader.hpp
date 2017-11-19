#ifndef GLSL_SHADER
#define GLSL_SHADER

#include <GL/glew.h>
#include <string>

class GLSLShader
{
public:
  GLSLShader(const std::string &filename, GLenum shaderType = GL_VERTEX_SHADER);
  GLSLShader(GLenum shaderType = GL_VERTEX_SHADER );
  ~GLSLShader();
  void compile();
  bool isCompiled() const; 
  void getShaderLog(std::string &log) const;
  void getShaderSource(std::string &shader) const;
  void setShaderSource(std::string &code);
  
  GLuint getHandle() const;
  void getParameter(GLenum param, GLint *data) const;

private:
  char *readShader(const std::string &filename, GLint * length);
  bool compiled_;
  GLint handle_;
};

#endif //GLSL_SHADER
