#ifndef guiOutput_h
#define guiOutput_h

struct Tracer;
struct Vortex;

void genFName(char *strBuffer, int frameNum);

void drawToConsole(struct Vortex *vorts, int numVorts, struct Tracer *tracers);
void drawToScreen(struct Vortex *vorts, int numVorts, struct Tracer *tracers);
void drawToFile(struct Vortex *vorts, int numVorts, struct Tracer *tracers, char filename[]);

#endif
