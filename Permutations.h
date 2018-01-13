#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H
#include <vector>

class Weave;

int reassignPermutations(Weave &weave);
void findSingularVertices(const Weave &weave, std::vector<int> &topologicalSingularVerts, std::vector<std::pair<int, int> > &geometricSingularVerts);

#endif

