#include <ctime>
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Functions.h"

using namespace std;
int main() {
    clock_t start = clock();
    FEMGrid grid;
    simulation(grid);
    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nTime elapsed: %.5f\n", elapsed);
    return 0;
}
