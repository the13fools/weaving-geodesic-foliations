#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H
#include <vector>

class Weave;

int reassignPermutations(Weave &weave);
void findSingularVertices(const Weave &weave, std::vector<int> &singularVerts);

#endif

