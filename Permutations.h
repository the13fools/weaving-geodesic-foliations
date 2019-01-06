#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H
#include <vector>

class Weave;

int reassignAllPermutations(Weave &weave);
int reassignCutPermutations(Weave &weave);
void findSingularVertices(const Weave &weave, std::vector<std::pair<int, int> > &topologicalSingularVerts, std::vector<std::pair<int, int> > &geometricSingularVerts);

#endif

